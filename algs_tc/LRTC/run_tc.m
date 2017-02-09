function T_hat = run_tc(params)
  T = params.T;
  Idx = params.Idx;
  
  T(T == 0) = 1e-3;
  %T = tensor(T);
  T = T.*Idx;
  %T = sptensor(T);
  Omega = logical(Idx);
  
  alpha = [1, 1, 1e-3];
  alpha = alpha / sum(alpha);
  maxIter = 100;
  epsilon = 1e-5;

  %% Fast LRTC (solve the relaxed formulation, FaLRTC algorithm in the paper)
  mu = 5 * alpha ./ sqrt(size(T));
  C =  0.6;
  L0 = 1e-5;
  [X_F, errList_F] = FaLRTC(...
     T,...                    % a tensor whose elements in Omega are used for estimating missing value
     Omega,...          % the index set indicating the obeserved elements
     alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     mu,...                 % the relaxation parameters, mu >= 0. See the function for definitions
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon...          % the tolerance of the relative difference of outputs of two neighbor iterations 
     );

   T_hat = X_F;
end
