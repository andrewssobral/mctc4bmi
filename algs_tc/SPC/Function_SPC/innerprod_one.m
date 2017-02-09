function u = innerprod_one(E,U,k)

  N = length(U);
  for m = 1:N
    V{m} = U{m}(:,k);
  end
  u = tensor_allprod(E,V,1);
