function H = Gaussian_window(L,sig)

  t1 = 1:L;
  t2 = 1:L;

  mu = [mean(t1) mean(t2)]';

  for ii = t1
  for jj = t2
  
    T = [ii; jj] - mu;
    H(ii,jj) = exp(-0.5*T'*T/sig/sig);

  end
  end
  H = H/sum(H(:));
