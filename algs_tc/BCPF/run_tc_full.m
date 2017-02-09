function T_hat = run_tc_full(params)
  T = params.T;
  Idx = params.Idx;
  
  DIM = size(T);
  R = 1;
  T(T == 0) = 1e-3;
  T = tensor(T);
  T = T.*Idx;
  %T = sptensor(T);
  
  % Bayes CP algorithm for incomplete tensor and tensor completion    
  [model] = BCPF_TC(T, 'obs', Idx, 'init', 'ml', ...
    'maxRank', max([DIM 2*R]), 'dimRed', 1, 'tol', 1e-6, 'maxiters', 100, 'verbose', 0);
  
  T_hat = double(model.X);
end
