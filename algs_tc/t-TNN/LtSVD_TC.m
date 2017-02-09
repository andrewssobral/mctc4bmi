% ===================================================================
%                 WTNN-Based Tensor Completion using ADMM
% ===================================================================
%
% We are trying to solve the following problem with acceleration.
%       
%          min ||X||_{WTNN}
%        s.t.  P_\Omega(X) = P_\Omega(M)
%
% The Lagrange:
% L(X,Z,Q) = \|Z\|_{TNN} + 1_{Y = P_\Omega(X)} + ...
%                    Q(:)^T (X(:)-Z(:)) +  rho/2 ||X-Z||_F^2
%
% ADMM algorithm is as follows:
%
% X^{k+1} = argmin_{X:Y = P_\Omega(X)} ||X-(Z^k - 1/rho*Q^k)||_F^2
%
% Z^{k+1} = argmin_Z  ||blkdiag(\hat{Z})||_* + ...
%           rho/(2n_3) ||\hat{Z} - (\hat{X} + 1/rho \hat{Q})||_F^2
%
% Q^{k+1} = Q^{k} + rho*( X^{k+1}-Z^{k+1})
%
% 
%  INPUTS:
%   
%   A - Matrix operator on unobserved X.
% 
%   y - Measurement vector A*X = y
%
%   alpha - over-relaxation parameter (typical values for alpha are 
%           between 1.0 and 1.8).
%
%   rho - augmented Lagrangian parameter. 
%
%   sX - size X object, i.e [n1 n2 n3 .... nP]
%
%   maxItr - maximum number of iterations

%   OUTPUTS:
%
%   x_hat - estimated solution.
%

function [x,RSE,History_RSE] = LtSVD_TC(A,y,rho,alpha,sX,maxItr,xref,QUIET,ImgShow)

%% Input checking.
if ~exist('QUIET','var')
    %显示迭代过程
    QUIET = false;
end

if ~exist('alpha','var')
    %显示迭代过程
    alpha = 1.05;
end

if ~exist('ImgShow','var')
    %显示补全中间过程
    ImgShow = true;
    sframe = round(sX(3)/2);  
end
%% set Default                          ;
    ABSTOL        =    1e-6                                ;
    RELTOL        =    1e-4                                ;


%% -----
t_start       =    tic                                 ;
if ImgShow
    figure('Name','Tensor Completion Process');
    nframe1 = sprintf('Original: %3dth frame', sframe);
    subplot(231); imshow(getframe(xref,sX,sframe)); title(nframe1);
    nframe2 = sprintf('Observed: %3dth frame', sframe);
    subplot(232); imshow(getframe(y,sX,sframe)); title(nframe2);
%     nframe3 = sprintf('Diff_init: %3dth frame', sframe);
%     subplot(233); imshow(getframe(abs(y-xref),sX,sframe)); title(nframe3);
end

%% ADMM solver

[~,n]         =    size(A)                             ;
x             =    zeros(n,1)                          ;
z             =    zeros(n,1)                          ;

u             =    zeros(n,1)                          ;

zhat          =    zeros(n,1)                          ;
uhat          =    zeros(n,1)                          ;

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective', 'RSE');
end

% precompute static variables for x-update (projection on to Ax=b)
% only diagonal form allowed
P             =    double(~logical(diag(A)))           ;
P             =    diag(sparse(P))                     ;
q             =    sparse(y)                           ;


for k         =     1 : maxItr
    
    rho = rho* alpha;
    zold      =     z                                  ;
    uold      =     u                                  ;

    x         =     P*(zhat - 1/rho*uhat) + q          ;     % x-update
    % z-update with relaxation
    [z, objV] =     wshrinkObj(x + 1/rho*uhat,1/rho,sX,0,1)   ;
    %
    u         =     uhat + rho*(x - z)                 ;
    
    r_norm    =     norm(x-z)                          ;
    s_norm    =     norm(-rho*(z - zold))              ;
    %% apply acceleration.
%     Ek = -1;
%     if k > 1
%        Ek = max(history.r_norm(k-1),...
%                history.s_norm(k-1))-max(r_norm,s_norm) ; 
%     end
   
%     if Ek > 0
%         alpha_old = alpha                              ;
%         alpha =  (1+sqrt(1+4*alpha^2))/2               ;
%         zhat  =  z + (alpha_old-1)/alpha*(z - zold)    ;
%         uhat  =  u + (alpha_old-1)/alpha*(u - uold)    ;       
%     else
%         alpha =  1                                     ;
        uhat  =  u                                     ;
        zhat  =  z                                     ; 
%     end

    if ~QUIET 
        %disp(k);
    end 
    %% diagnostics, reporting, termination checks
    history.objval(k)   =  objV                        ;

    history.r_norm(k)   =  r_norm                      ;
    history.s_norm(k)   =  s_norm                      ;
    
    history.eps_pri(k)  = (sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z)))*5;
    history.eps_dual(k) = (sqrt(n)*ABSTOL +RELTOL*norm(rho*u))*5           ;
%     history.eps_pri(k)  = (sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z)));
%     history.eps_dual(k) = (sqrt(n)*ABSTOL +RELTOL*norm(rho*u))          ;
       
    history.rse(k)      = rse(x,xref);
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%5d\t%3d\n', k, ...
           history.r_norm(k), history.eps_pri(k), history.s_norm(k),...
            history.eps_dual(k),history.objval(k),history.rse(k))      ;
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
        history.s_norm(k) < history.eps_dual(k))
        break;
    end
    
    if ImgShow
        name_x = sprintf('Recons:%sth iter', num2str(k));
        subplot(234), imshow(getframe(x,sX,sframe)), title(name_x);
        name_z = sprintf('Z:%sth iter', num2str(k));
        subplot(235), imshow(getframe(z,sX,sframe)), title(name_z);
        name_dif = sprintf('Difference:%sth iter', num2str(k));
        subplot(236), imshow(getframe(abs(x-xref),sX,sframe)), title(name_dif);
        name_dif = sprintf('-RSE Curve', num2str(k));
        subplot(233), plot(history.rse,'r','linewidth',1.5), title(name_dif);
        ylabel('-RSE');xlabel('iter');axis auto;
        pause(0.1);
    end
end

if ~QUIET
    toc(t_start);
end
History_RSE = history.rse;
RSE = history.rse(k);
end