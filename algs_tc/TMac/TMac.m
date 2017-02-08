function [X,Y,Out] = TMac(data,known,Nway,coreNway,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tensor completion by parallel matrix factorization
% ==============================================================================
% solve
%  min_{X,Y} sum_{n=1}^N 0.5*alpha(n)*||Proj_{known}(Xn*Yn-Unfold_n(M))||_F^2
% ===============================================================================
%
% Input:
%       data: observed entries of the underlying tensor
%       known: indices of observed entries
%       Nway: the dimension of the underlying tensor
%       coreNway: estimated ranks of all mode matricizations
%       opts.
%           maxit: maximum number of iterations (default: 500)
%           tol: stopping tolerance (default: 1e-4)
%           maxT: maximum running time (sec) (default: 1e6)
%           alpha: weights in the model (default: alpha(n) = 1/N, any n)
%           alpha_adj: determine whether dynamically update alpha
%                       (default: 1 (yes))
%           rank_adj: determine rank-adjusting strategy
%                       (1: increase; -1: decrease; 0: fix)
%           rank_inc: rank increment if rank-increasing strategy is used
%           rank_min: minimum rank estimation
%           rank_max: maximum rank estimation
%
% Output:
%       X,Y: cell structs
%       Out.
%           iter: number of iterations
%           relerr1: relative change array of total fitting
%           relerr2: array of total fitting
%           alpha: final weights alpha (may be different from input)

N = length(Nway); coNway = zeros(1,N);
for n = 1:N
    coNway(n) = prod(Nway)/Nway(n);
end
nrmb = norm(data);
%% Parameters and defaults
if isfield(opts,'maxit')      maxit = opts.maxit;     else maxit = 500;             end
if isfield(opts,'tol')        tol = opts.tol;         else tol = 1e-4;              end
if isfield(opts,'maxT')       maxT = opts.maxT;       else maxT = 1e6;              end
if isfield(opts,'alpha')      alpha = opts.alpha;     else alpha = ones(1,N)/N;     end
if isfield(opts,'alpha_adj') 
    alpha_adj = opts.alpha_adj;
else
    alpha_adj = 1;
end

if isfield(opts,'rank_adj') rank_adj = opts.rank_adj; else rank_adj = zeros(1,N);   end
if isfield(opts,'rank_inc') rank_inc = opts.rank_inc; else rank_inc = ones(1,N);    end
if isfield(opts,'rank_min') rank_min = opts.rank_min; else rank_min = ones(1,N);    end
if isfield(opts,'rank_max') rank_max = opts.rank_max; else rank_max = 50*ones(1,N); end

%% Data preprocessing and initialization
if isfield(opts,'X0')
    X = opts.X0;
else
    X = cell(1,N);
    for n = 1:N
        X{n} = randn(Nway(n),coreNway(n));
    end
end

if isfield(opts,'Y0')
    Y = opts.Y0;
else
    Y = cell(1,N);
    for n = 1:N
        Y{n} = randn(coreNway(n),coNway(n));
    end
end

% rescale the initial point based on number of elements
estMnrm = sqrt(nrmb^2*(prod(Nway)/length(known)));

for n = 1:N
    X{n} = X{n}/norm(X{n},'fro')*estMnrm^(Nway(n)/(Nway(n)+coNway(n)));
    Y{n} = Y{n}/norm(Y{n},'fro')*estMnrm^(coNway(n)/(Nway(n)+coNway(n)));
end

X0 = X; Y0 = Y;

[known,id] = sort(known); data = data(id);

M = zeros(Nway); M(known) = data;


sx = cell(1,N);
reschg = ones(1,N); 
reschg_tol = max(1e-2,10*tol);
rank_inc_num = sum(rank_adj==1);

res0 = zeros(1,N); TotalRes = 0;
res = res0;
for n = 1:N
    Mn = Fold(X{n}*Y{n},Nway,n);
    res0(n) = norm(Mn(known)-data);
    TotalRes = TotalRes+res0(n);
end
solX = ones(1,N);
Xsq = cell(1,N); Yt = cell(1,N); spI = cell(1,N);
for n = 1:N
    Yt{n} = Y{n}';
end

Out.rank = coreNway;

start_time = tic;

alpha = alpha/sum(alpha);

fprintf('Iteration:     ');

for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k);
    
    % update (X,Y)
    for n = 1:N
        if alpha(n) > 0
            Mn = Unfold(M,Nway,n);
            if solX(n)   
                X{n} = Mn*Yt{n};
            end
            solX(n) = 1;

            Xsq{n} = X{n}'*X{n};
            Y{n} = pinv(Xsq{n})*X{n}'*Mn;
            Yt{n} = Y{n}';

            if rank_adj(n) == -1
                rank_dec_adaptive();
            end
        end
    end
    
    % update M
    Mn = Fold(X{1}*Y{1},Nway,1);
    res(1) = norm(Mn(known)-data);
    M = alpha(1)*Mn;
    for n = 2:N
        if alpha(n) > 0
        Mn = Fold(X{n}*Y{n},Nway,n);
        res(n) = norm(Mn(known)-data);
        M = M+alpha(n)*Mn;
        end
    end
    
    % pass the true tensor M for evaluation
    if isfield(opts,'Mtr')
        Out.truerel(k) = norm(M(:)-opts.Mtr(:))/norm(opts.Mtr(:));
    end
    
    M(known) = data;
    
    TotalRes0 = TotalRes;
    TotalRes = 0;
    for n = 1:N
        if alpha(n) > 0
            TotalRes = TotalRes+res(n)^2;
        end
    end
    ratio = res./res0;   reschg = abs(1-ratio);    
     
    if rank_inc_num > 0
        for n = 1:N
            if alpha(n) > 0
            if coreNway(n) < rank_max(n) && reschg(n) < reschg_tol
                rank_inc_adaptive();
            end
            end
        end
    end
    
    % adaptively update weight
    
    if alpha_adj ~= 0
        alpha = 1./(res.^2); alpha = alpha/sum(alpha);
    end
    
    % record how the rank estimates are updated
    Out.rank = [Out.rank; coreNway];
    
    % --- diagnostics, reporting, stopping checks ---
    relerr1 = abs(TotalRes-TotalRes0)/(TotalRes0+1); 
    relerr2 = sum(alpha.*res)/nrmb;
    
    % reporting
    Out.hist_rel(1,k) = relerr1;
    Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion
    crit = relerr1<tol;
    if crit; nstall = nstall+1; else nstall = 0; end
    if nstall>=3 || relerr2<tol break; end
    if toc(start_time)>maxT; break; end;
    
    X0 = X; Y0 = Y; M0 = M;
    res0 = res;
end % main
fprintf('\n'); Out.iter = k;
Out.alpha = alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rank_dec_adaptive()
        % decrease the estimated rank
        sx{n} = svd(Xsq{n});
        dR = sx{n}(rank_min(n):end);
        drops = dR(1:end-1)./dR(2:end);
        [dmx,imx] = max(drops);
        rel_drp = (coreNway(n)-rank_min(n))*dmx/(sum(drops)-dmx);
        %If a large drop is found, adjust the rank estimate to imx
        if rel_drp>10
            coreNway(n) = imx+rank_min(n)-1;
            % set rank_adj(n) to 0, so only decrease the rank once
            rank_adj(n) = 0;
            [Qx,Rx] = qr(X{n},0);
            [Qy,Ry] = qr(Y{n}',0);
            [U,S,V] = svd(Rx*Ry');
            sigv = diag(S);
            X{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),...
                0,coreNway(n),coreNway(n)); 
            X0{n} = X{n};
            
            Yt{n} = Qy*V(:,1:coreNway(n));
            Y{n} = Yt{n}';
            Y0{n} = Y{n};
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rank_inc_adaptive()
        % increase the estimated rank
                   
        [Q,R] = qr(Y{n}',0);
        for ii = 1:rank_inc(n)
            rdnx = randn(coNway(n),1);
            rdnx = rdnx-Q*(Q'*rdnx);
            rdnx = rdnx/norm(rdnx);
            Q = [Q,rdnx];
        end
        Y{n} = Q'; Yt{n} = Q; Y0{n} = Q';
        
        coreNway(n) = coreNway(n)+rank_inc(n);
        if coreNway(n) >= rank_max(n)
            rank_inc_num = rank_inc_num - 1;
        end

        if rank_inc_num == 0
            nstall = 0;
        end

        Mn = Unfold(M,Nway,n);
        X{n} = Mn*Y{n}';
        X0{n} = X{n};
        solX(n) = 0;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function W = Fold(W, dim, i)
        dim = circshift(dim, [1-i, 1-i]);
        W = shiftdim(reshape(W, dim), length(dim)+1-i);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function W = Unfold(W, dim, i)
        W = reshape(shiftdim(W,i-1), dim(i), []);
    end

end    

