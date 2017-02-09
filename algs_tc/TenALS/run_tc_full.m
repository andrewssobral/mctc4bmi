function T_hat = run_tc_full(params)
  T = params.T;
  Idx = params.Idx;
  T_hat = zeros(size(T));
  [U1,U2,U3,S,~] = TenALS(T, Idx, 2, [], [], []);
  for i = 1:size(T,3)
    T_hat(:,:,i) = U1*diag(U3(i,:).*S(:)')*U2';
  end
end
