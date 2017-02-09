% ===================================================================
%                  Tensor Completion using ADMM
% ===================================================================
%
% We are trying to solve the following problem with acceleration.
%       
%          min ||X||_{TNN}
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
% W^{k+1} = W^{k} + rho*( L^{k+1}+S^{k+1}-M )
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
%
%   myNorm - string of one of the supported norms.  see shrinkObj.
%  
%
%   parOp - (optional bool) run the algorithm in parallel
%
%   OUTPUTS:
%
%   x_hat - estimated solution.
%
%   Authors: G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015


function [x,history] = tensor_cpl_admm(A,y,rho,alpha,sX,maxItr,myNorm,QUIET,varargin)

%% Input checking.
if ~exist('QUIET','var')
    QUIET = false;
end

%set Default
parOP         =    false                               ;
ABSTOL        =    1e-6                                ;
RELTOL        =    1e-4                                ;


if ~exist('vararg','var')
    nARG      =    length(varargin)                    ;
    if nARG   >    0
        parOP =    varargin{1}                         ;
    end
    
    if nARG   >    1
        lip   =    varargin{2}                         ;
        
    end
    
    
end
t_start       =    tic                                 ;

%% ADMM solver

[~,n]         =    size(A)                             ;
x             =    zeros(n,1)                          ;
z             =    zeros(n,1)                          ;

u             =    zeros(n,1)                          ;

zhat          =    zeros(n,1)                          ;
uhat          =    zeros(n,1)                          ;

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

% precompute static variables for x-update (projection on to Ax=b)
% only diagonal form allowed
P             =    double(~logical(diag(A)))           ;
P             =    diag(sparse(P))                     ;
q             =    sparse(y)                           ;


for k         =     1 : maxItr
    
    zold      =     z                                  ;
    uold      =     u                                  ;

    x         =     P*(zhat - 1/rho*uhat) + q          ;     % x-update
    % z-update with relaxation
    [z, objV] =     shrinkObj(x + 1/rho*uhat,...
                          1/rho,myNorm,sX,parOP)   ;
    %
    u         =     uhat + rho*(x - z)                 ;
    
    r_norm    =     norm(x-z)                          ;
    s_norm    =     norm(-rho*(z - zold))              ;
    %% apply acceleration.
    Ek = -1;
    if k > 1
       Ek = max(history.r_norm(k-1),...
               history.s_norm(k-1))-max(r_norm,s_norm) ; 
    end
   
    if Ek > 0
        alpha_old = alpha                              ;
        alpha =  (1+sqrt(1+4*alpha^2))/2               ;
        zhat  =  z + (alpha_old-1)/alpha*(z - zold)    ;
        uhat  =  u + (alpha_old-1)/alpha*(u - uold)    ;
        
    else
        alpha =  1                                     ;
        uhat  =  u                                     ;
        zhat  =  z                                     ; 
    end

    if ~QUIET 
        disp(k);
    end 
    %% diagnostics, reporting, termination checks
    history.objval(k)   =  objV                        ;

    history.r_norm(k)   =  r_norm                      ;
    history.s_norm(k)   =  s_norm                      ;
    
    history.eps_pri(k)  = sqrt(n)*ABSTOL + ...
                          RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k) = sqrt(n)*ABSTOL + ...
                          RELTOL*norm(rho*u)           ;
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%5d\n', k, ...
           history.r_norm(k), history.eps_pri(k), ...
           history.s_norm(k), history.eps_dual(k),...
                              history.objval(k) )      ;
    end

    if (history.r_norm(k) < history.eps_pri(k) && ...
        history.s_norm(k) < history.eps_dual(k))
        break;
    end
end

if ~QUIET
    toc(t_start);
end

end