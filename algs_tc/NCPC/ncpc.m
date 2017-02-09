function [A,Out] = ncpc(data,known,ndim,r,opts)
% ncpi: nonnegative tensor decomposition (CANDECOMP/PARAFAC) from partial data by block-coordinate
%  min 0.5*||P_known(M - A_1 \circ ... \circ A_N)||_F^2  
%  subject to A_1>=0, ..., A_N>=0
%
%input:
%       data: observed data of a nonnegative tensor
%       known: observed index set
%       ndim: size of the underlying tensor
%       r: estimated rank; require exact or moderate overestimates
%       opts.
%           tol: tolerance for relative change of function value, default:
%           1e-4
%           maxit: max number of iterations, default: 500
%           maxT: max running time, default: 1e3
%           rw: control the extrapolation weight, default: 1
%           A0: initial point in cell struct, default: Gaussian random
%           matrices
% output:
%       A: nonnegative ktensor
%       Out.
%           iter: number of iterations
%           hist_obj: history of objective values
%           hist_rel: history of rrelative objective changes (row 1) and relative residuals (row 2)
%
% require MATLAB Tensor Toolbox from
% http://www.sandia.gov/~tgkolda/TensorToolbox/
%
% More information can be found at:
% http://www.caam.rice.edu/~optimization/bcu/

%% Parameters and defaults
if isfield(opts,'tol');   tol = opts.tol;     else tol = 1e-4;   end % stopping tolerance
if isfield(opts,'maxit'); maxit = opts.maxit; else maxit = 500;  end % max # of iterations
if isfield(opts,'maxT');  maxT = opts.maxT;   else maxT = 1e3;   end % max time in seconds 
if isfield(opts,'rw');    rw = opts.rw;       else rw = 1;       end % initial extrapolation weight

%% Data preprocessing and initialization
[known,I] = sort(known);
data = data(I); % data contains the observed entries
N = length(ndim); % the underlying tensor is N-way, and its dimensions are given in ndim
M = zeros(ndim); M(known) = data; M0 = M; % the initial estimated tensor M
nrmb = norm(data); % norm of observed entries
obj0 = 0.5*nrmb^2; % initial objective value

sizeN = zeros(1,N);
for n = 1:N; sizeN(n) = prod(ndim)/ndim(n); end

% initial tensor factors
if isfield(opts,'A0') 
    A0 = opts.A0; 
else
    A0 = cell(1,N);
    for n = 1:N
        A0{n} = max(0,randn(ndim(n),r)); % randomly generate each factor
    end
end

% normalize A0 and cache its square
Asq = cell(1,N);
for n = 1:N
    A0{n} = A0{n}/norm(A0{n},'fro')*nrmb^(1/N);
    Asq{n} = A0{n}'*A0{n};
end
Am = A0; A = A0; % see below for the definitions of A, A0, Am

nstall = 0; % # of stalled iterations
t0 = 1; % used for extrapolation weight update
wA = ones(N,1); % extrapolation weight array
L0 = ones(N,1); L = ones(N,1); % Lipschitz constant array

%% Store data to save computing time if it is not too large
if N*prod(ndim)<4000^2
    storedata = true;   kroA = cell(1,N);
    for n = 1:N
        kroA{n} = zeros(sizeN(n),r);
    end
else
    storedata = false;
    kroA1 = zeros(prod(ndim)/ndim(1),r);
end

% kroA{1} = A{N}\odot...\odot A{n+1}\odot A{n-1}...\odot A{2}
% "odot" denotes Khatri-Rao-product
matorder = N:-1:2;
for j = 1:r
    ab = A{matorder(1)}(:,j);
    for i = matorder(2:N-1)
        ab = A{i}(:,j) * ab(:).';
    end
    if storedata
        kroA{1}(:,j) = ab(:);
    else
        kroA1(:,j) = ab(:);
    end
end

%% Iterations of block-coordinate update 
%
%  iteratively updated variables:
%       M: estimated tensor
%       Gn: gradients with respect to A{n}
%       A: current factors
%       A0: previous factors
%       Am: extrapolated factors
%       L, L0: current and previous Lipschitz bounds
%       obj, obj0: current and previous objective values

start_time = tic;
fprintf('Iteration:     ');

for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k); 
    
    %--- update each factor matrix ----%
    % For derivation of the update, see Section 3.2 of this report:
    % Yangyang Xu and Wotao Yin. A block coordinate descent method for 
    % regularized multi-convex optimization with application to nonnegative
    % tensor factorization and completion. 
    % Rice University CAAM Technical Report TR12-15
    
    %--- update the first factor matrix ----%
    Bsq = Asq{2};
    for i = 3:N            
        Bsq = Bsq.*Asq{i};
    end
    L0(1) = L(1);     L(1) = norm(Bsq);
    % form the 1-mode matricization of M
    MB = reshape(M, ndim(1), sizeN(1));
    if storedata        
        MB = MB*kroA{1};
    else
        MB = MB*kroA1;
    end
    % compute the gradient
    Gn = Am{1}*Bsq-MB; 
    A{1} = max(0,Am{1}-Gn/L(1));
    Asq{1} = A{1}'*A{1};
    
    %--- update each other factor matrix ----%
    for n = 2:N
        matorder = [2:n-1,n+1:N];
        % Bsq = Asq{1}.*Asq{2}...Asq{n-1}.*Asq{n+1}....*Asq{N}
        Bsq = Asq{1};
        for i = matorder
            Bsq = Bsq.*Asq{i};            
        end
        L0(n) = L(n);     L(n) = norm(Bsq);
        if storedata
            % do Khatri-Rao-product
            % kroA{n} = A{N}\odot...\odot A{n+1}\odot A{n-1}...\odot A{1}
            matorder = [N:-1:n+1,n-1:-1:1];
            for j = 1:r
                ab = A{matorder(1)}(:,j);
                for i = matorder(2:N-1)
                    ab = A{i}(:,j) * ab(:).';
                end
                kroA{n}(:,j) = ab(:);
            end
            % do matricized-tensor-times-Khatri-Rao-product
            % which equals n-th mode matricization of M times kroA
            MB = permute(M,[n,1:n-1,n+1:N]);
            MB = reshape(MB, ndim(n), sizeN(n));
            MB = MB*kroA{n};
        else
            % do matricized-tensor-times-Khatri-Rao-product
            MB = mttkrp(tensor(M),A,n);
        end
        % compute the gradient
        Gn = Am{n}*Bsq-MB;
        A{n} = max(0,Am{n}-Gn/L(n));
        Asq{n} = A{n}'*A{n};
    end
    
    % --- update the estimate of the underlying tensor ----%
    matorder = N:-1:2;
    for j = 1:r
        ab = A{matorder(1)}(:,j);
        for i = matorder(2:N-1)
            ab = A{i}(:,j) * ab(:).';
        end
        if storedata
            % kroA{1} = A{N}\odot...\odot A{n+1}\odot A{n-1}...\odot A{2}
            kroA{1}(:,j) = ab(:);
        else
            kroA1(:,j) = ab(:);
        end
    end
    if storedata; M = A{1}*kroA{1}'; else M = A{1}*kroA1'; end
    
    % --- diagnostics, reporting, stopping checks ---
    obj = 0.5*norm(M(known)-data)^2;
    M(known) = data; M = reshape(M,ndim);
    relerr1 = abs(obj-obj0)/(obj0+1);    relerr2 = (2*obj)^.5/nrmb;
    
    % reporting
    Out.hist_obj(k) = obj;
    Out.hist_rel(1,k) = relerr1;      
    Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion  
    crit = relerr1<tol; 
    if crit; nstall = nstall+1; else nstall = 0; end
    if nstall>=3 || relerr2<tol; break; end
    if toc(start_time)>maxT; break; end;
    
    % --- correction and extrapolation ---
    t = (1+sqrt(1+4*t0^2))/2;
    if obj>=obj0 
        % restore A,M to make the objective nonincreasing
        Am = A0; M = M0;
    else
        % apply extrapolation
        w = (t0-1)/t; % extrapolation weight
        for n = 1:N
            wA(n) =  min([w,rw*sqrt(L0(n)/L(n))]); % choose smaller weight for convergence
            Am{n} = A{n}+wA(n)*(A{n}-A0{n}); % extrapolation
        end
        A0 = A; M0 = M; t0 = t; obj0 = obj;
    end    
end
fprintf('\n');
A = ktensor(A);
Out.iter = k; % report # of iterations

%% <http://www.caam.rice.edu/~optimization/bcu/ncpc/ncpc.m Download this m-file>