function T_hat = run_tc(params)
  T = params.T;
  Idx = params.Idx;
  
  T(T == 0) = 1e-3;
  %T = tensor(T);
  T = T.*Idx;
  %T = sptensor(T);
  Q = logical(Idx);
  
  %% hyperparameters and run SPC-TV

  TVQV    = 'tv';        % 'tv' or 'qv' ;
  rho     = [0.01 0.01 0.01]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
  K       = 10;          % Number of components which are updated in one iteration.
  SNR     = 50;          % error bound
  nu      = 0.01;        % threshold for R <-- R + 1.
  maxiter = 100;       % maximum number of iteration
  tol     = 1e-7;        % tolerance
  out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.

  [Xtv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);
  T_hat = Xtv;

  %% hyperparameters and run SPC-QV

%   TVQV    = 'qv';        % 'tv' or 'qv' ;
%   rho     = [1.0 1.0 1.0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
%   K       = 10;          % Number of components which are updated in one iteration.
%   SNR     = 50;          % error bound
%   nu      = 0.01;        % threshold for R <-- R + 1.
%   maxiter = 100;       % maximum number of iteration
%   tol     = 1e-7;        % tolerance
%   out_im  = 0;           % you can monitor the process of 'image' completion if out == 1.
% 
%   [Xqv Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);
%   T_hat = Xqv;

end
