function T_hat = run_tc_full(params)
  T = params.T;
  Idx = params.Idx;
  
  T(T == 0) = 1e-3;
  T = tensor(T.*Idx);
  T = sptensor(T);
  %Idx = sptensor(Idx);
  
  n = size(T);
  r1 = n(1);
  r2 = n(2);
  r3 = 1;
  r = [r1, r2, r3];
  T_init = makeRandTensor(n, r);
  
  opts = struct('maxiter', 100, 'tol', 1e-6, 'verbose', false);
  [resX, err, conv] = geomCG(T, T_init, [], opts);
  
  T_hat = double(resX);
end
