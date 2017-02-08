function T_hat = run_tc(params)
  T = params.T;
  Idx = params.Idx;
  
  Omega = find(Idx);
  Ak = T(Omega);
  
  N = ndims(T);
  Nway = size(T); % dimension of tensor
  coreNway = [Nway(1),Nway(2),1];
  
  % rank_dec strategy
  opts = [];
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
  [X_dec,Y_dec,Out_dec] = TMac(Ak,Omega,Nway,EstCoreNway,opts);
  % use the weighted sum of all mode matrix factorizations as estimated tensor
  T_hat = zeros(Nway);
  for i = 1:N
    T_hat = T_hat+Out_dec.alpha(i)*Fold(X_dec{i}*Y_dec{i},Nway,i);
  end
end
