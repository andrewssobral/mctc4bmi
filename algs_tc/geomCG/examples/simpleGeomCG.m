% Simple example of the geomCG tensor completion Code.
%
% GeomCG Tensor Completion. Copyright 2013 by
% Michael Steinlechner
% Questions and contact: michael.steinlechner@epfl.ch
% BSD 2-clause license, see ../LICENSE.txt

n = [100 150 200];
r1 = 4;
r2 = 5;
r3 = 6;
r = [r1, r2, r3];

% sparsity: p = nnz / (n^3)
p = 0.01;

opts = struct( 'maxiter', 100, 'tol', 1e-9 );

% set the seed of the mersenne twister to 11 for reproducible results
rng( 11 );

% create the data to complete:
% create random rank-r tensor ...
A = makeRandTensor( n, r );
% create the sampling set ...
subs = makeOmegaSet( n, round(p*prod(n)) );
% get the values of A at the sampling points ...
vals = getValsAtIndex(A, subs);
% save indices and values together in a sparse tensor
A_Omega = sptensor( subs, vals, n, 0);

% create the test set to compare:
subs_Test = makeOmegaSet( n, round(p*prod(n)) );
% get the values of A at the test points ...
vals_Test = getValsAtIndex(A, subs_Test);
% save indices and values together in a sparse tensor
A_Test = sptensor( subs_Test, vals_Test, n, 0);

% random initial guess:
X_init = makeRandTensor( n, r );

% RUN THE TENSOR COMPLETION ALGORITHM
% -----------------------------------
[~, err, ~] = geomCG( A_Omega, X_init, A_Test, opts);

% Plot the results.
set(gca,'fontsize',14)
semilogy( err(1,:),'-or','Linewidth',1.5);
hold on
semilogy( err(2,:),'-xb','Linewidth',1.5);
hold off
xlabel('Iteration')
ylabel('Rel. residual and rel. error on test set')
legend('Rel. error on \Omega', 'Rel. error on \Gamma')
title(sprintf('Tensor completion with %i%% sampling. n = 100, r = [4, 5, 6]', round(100*p)))
set(figure(1), 'Position', [0 0 600 500])
