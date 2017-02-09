function T_hat = run_tc_full(params)
  T = params.T;
  Idx = params.Idx;
  
  ndim = size(T);
  known = find(Idx);
  data = T(known);
  R = 2; % tensor rank
  esr = round(1.25*R); % overestimate of rank
  opts.tol = 1e-4; opts.maxit = 1000;
  [T_hat] = ncpc(data,known,ndim,esr,opts);
  
  T_hat = double(T_hat);
end
