function Z = tensor_allprod(G,U,tr)

  N = ndims(G);
  N = length(U);

  Z = G;
  for n = 1:N
    if ~isempty(U{n})
      if tr == 0
        Z = tmult(Z,U{n},n);
      else
        Z = tmult(Z,U{n}',n);
      end
    end
  end
