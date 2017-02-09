% An example of nonnegative CP decomposition from partial observations
%% Generate synthetic 4-order tensor
rand('seed',0); randn('seed',0);
N1 = 30; N2 = 30; N3 = 30; N4 = 30; % tensor dimension
R = 10; % tensor rank

% randomly generate factor matrices
A1 = max(0,randn(N1,R));
A2 = max(0,randn(N2,R));
A3 = max(0,randn(N3,R));
A4 = max(0,randn(N4,R));
% generate tensor M using the above factor matrices
M = ktensor({A1,A2,A3,A4});
M = full((arrange(M)));

sr = 0.4; % percentage of samples
% randomly choose samples
known = randsample(N1*N2*N3*N4,round(sr*N1*N2*N3*N4));
data = M.data(known);

sn = 60; % signal to noise
% -- add noise --
Noise = max(0,randn(size(data)));
data = data + 10^(-sn/20)*norm(data)/norm(Noise)*Noise;
ndim = [N1,N2,N3,N4];

%% Solve problem
% exact rank works but a rough rank overestimate is also fine
esr = round(1.25*R); % overestimate of rank
opts.tol = 1e-4; opts.maxit = 1000;
t0 = tic;
[A,Out] = ncpc(data,known,ndim,esr,opts);
time = toc(t0);

%% Reporting
relerr = norm(full(A)-M)/norm(M);
fprintf('time = %4.2e, ',time);
fprintf('solution relative error = %4.2e\n\n',relerr);

figure;
semilogy(1:Out.iter, Out.hist_obj,'k-','linewidth',2);
xlabel('iteration','fontsize',12);
ylabel('objective value','fontsize',12)

figure;
semilogy(1:Out.iter, Out.hist_rel(2,:),'k-','linewidth',2);
xlabel('iteration','fontsize',12);
ylabel('relative residual','fontsize',12)

%% <http://www.caam.rice.edu/~optimization/bcu/ncpc/example_ncpc.m Download this m-file>