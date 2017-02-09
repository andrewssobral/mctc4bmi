% Tensor completion with noisy entries.
%
% This example reproduces Figure 5 in 
%       D. Kressner, M. Steinlechner, B. Vandereycken:
%       Low-Rank Tensor Completion by Riemannian Optimization.
%       MATHICSE Technical Report 20.2013, June 2013. Submitted to BIT Numerical Mathematics.

% GeomCG Tensor Completion. Copyright 2013 by
% Michael Steinlechner
% Questions and contact: michael.steinlechner@epfl.ch
% BSD 2-clause license, see ../LICENSE.txt

n = [100 100 100];
k1 = 6;
k2 = 6;
k3 = 6;
k = [k1, k2, k3];

% sparsity: p = nnz / (n^3)
p = 0.1;

opts = struct( 'maxiter', 30, 'tol', 1e-15, 'testtol', 1e-15, 'verbose', true);

% set the seed of the mersenne twister to 11 for reproducible results
rng( 11 );

A = makeRandTensor( n, k );
subs = makeOmegaSet( n, round(p*prod(n)) );
vals = getValsAtIndex(A, subs);

% random initial guess:
X_init = makeRandTensor( n, k );

% setup plot colors
noise = [0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4];
color = [0,0,0;lines(5)];
style = {'-k.','-ob','-xr','-*g','-cs','-md'};
style = {'-k.','-ob','-xr','-*g','-cs','-md'};
style2 = {'--k.','--ob','--xr','--*g','--cs','--md'};
levels = {'--k','--b','--r','--g','--c','--m'};

set(gca,'fontsize',14)
norm_A_Omega = norm(vals);
N = tensor(randn(n));
norm_N_Omega = norm(N(subs));

for i=1:length(noise)
    e = noise(i);
    % create noisy observations
    A_Omega = sptensor( subs, vals + e*norm_A_Omega/norm_N_Omega*N(subs), n, 0);

    % complete the tensor
    [~, err_Omega, ~] = geomCG( A_Omega, X_init, [], opts)

    semilogy( err_Omega(1,:), style{i},'Linewidth',1.5,'color',color(i,:));
    hold on
    semilogy( [1, opts.maxiter], [noise(i), noise(i)], levels{i},'color',color(i,:));
    drawnow;
end

xlabel('Iteration')
ylabel('Rel. residual and noise levels')
title('Noisy completion, n = 100, r = [6, 6, 6]')
set(gca,'Ytick',[1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0])
ylim([1e-15,10])
set(figure(1), 'Position', [0 0 600 500])
