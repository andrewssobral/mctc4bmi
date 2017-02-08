function quicktest
%% Generate synthetic 3-order tensor
clear; close all;
rand('seed',2013); randn('seed',2013);
Nway = [50,50,50]; % dimension of tensor
r = 10;
coreNway = [r,r,r];

% randomly generate core tensor
G = tensor(randn(coreNway));
A = cell(1,ndims(G));
% randomly generate factor matrices
for i = 1:ndims(G)
    A{i} = randn(Nway(i),coreNway(i));
end
% generate tensor
M = full(ttensor(G,A));
% M = tensor(M.data/max(M.data(:)));
N = ndims(M);
M = M.data;

sr = 0.5;
p = round(sr*prod(Nway));
known = randsample(prod(Nway),p); data = M(known);

opts = [];
opts.maxit = 100; 
opts.tol = -1e-5; % run to maxit by using negative tolerance
opts.Mtr = M; % pass the true tensor to calculate the fitting

%% rank_dec strategy
opts.alpha_adj = 0;
opts.rank_adj = -1*ones(1,3);
opts.rank_min = 5*ones(1,3);
opts.rank_max = 20*ones(1,3);

EstCoreNway = round(1.25*coreNway);
coNway = zeros(1,N);
for n = 1:N
    coNway(n) = prod(Nway)/Nway(n);
end

% use random generated starting point
for i = 1:3
    X0{i} = randn(Nway(i),EstCoreNway(i));
    Y0{i} = randn(EstCoreNway(i),coNway(i));
end

opts.X0 = X0; opts.Y0 = Y0;
t0 = tic;
[X_dec,Y_dec,Out_dec] = TMac(data,known,Nway,EstCoreNway,opts);
time = toc(t0);

% use the weighted sum of all mode matrix factorizations as estimated tensor
Mrec = zeros(Nway);
for i = 1:N
    Mrec = Mrec+Out_dec.alpha(i)*Fold(X_dec{i}*Y_dec{i},Nway,i);
end

% Reporting
relerr = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('rank_dec:\n')
fprintf('SR-%3.2f: time = %4.2e, ',sr,time);
fprintf('relative error = %4.2e\n\n',relerr);

%% rank_inc strategy
opts.alpha_adj = 0;
opts.rank_adj = 1*ones(1,3);
opts.rank_min = 5*ones(1,3);
opts.rank_max = round(1.5*coreNway);

EstCoreNway = round(0.75*coreNway);
coNway = zeros(1,N);
for n = 1:N
    coNway(n) = prod(Nway)/Nway(n);
end

for i = 1:3
    X0{i} = randn(Nway(i),EstCoreNway(i));
    Y0{i} = randn(EstCoreNway(i),coNway(i));
end

opts.X0 = X0; opts.Y0 = Y0;
t0 = tic;
[X_inc,Y_inc,Out_inc] = TMac(data,known,Nway,EstCoreNway,opts);
time = toc(t0);

% use the weighted sum of all mode matrix factorizations as estimated tensor
Mrec = zeros(Nway);
for i = 1:N
    Mrec = Mrec+Out_inc.alpha(i)*Fold(X_inc{i}*Y_inc{i},Nway,i);
end

% Reporting
relerr = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('rank_inc:\n')
fprintf('SR-%3.2f: time = %4.2e, ',sr,time);
fprintf('relative error = %4.2e\n\n',relerr);

%%%
% plot the fitting with respect to iterations
fig = figure('papersize',[5,4],'paperposition',[0,0,5,4]);
semilogy(Out_dec.truerel,'r-.','linewidth',4)
hold on;
semilogy(Out_inc.truerel,'k-','linewidth',4)
legend('TMac-dec','TMac-inc','location','best')
set(gca,'fontsize',12)
xlabel('iteration','fontsize',16)
ylabel('relative error','fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function W = Fold(W, dim, i)
        dim = circshift(dim, [1-i, 1-i]);
        W = shiftdim(reshape(W, dim), length(dim)+1-i);
    end
end