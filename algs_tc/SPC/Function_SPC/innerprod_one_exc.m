function u = innerprod_one_exc(E,U,k,n)

  N = length(U);
  for m = 1:N
    V{m} = U{m}(:,k);
  end
  u = tensor_allprod_exc(E,V,1,n);
