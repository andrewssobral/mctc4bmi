function T_hat = run_tc_full(params)
  T = params.T;
  Idx = params.Idx;
  
  %T(T == 0) = 1e-6;
  %T = T.*Idx;
  %T = sptensor(T);
  %T = tensor(T);
  
  Omega = find(Idx);
  data.linInd = Omega;
  data.b = T(Omega);
  data.T = T;
  data.X = T;
  N = ndims(data.T);
  r = 1/sqrt(max(size(data.T))); 
  params.E0 = tenzeros(size(data.T));
  params.X0 = tenzeros(size(data.T));
  params.V0 = cell(1, N);
  for i = 1:N 
    params.V0{i} = tenzeros(size(data.T));
  end
  params.mu0 = 1/(N+1);
  params.mode = N;
  params.IsTC = true; % is tensor completion
  params.rRatio = 1/4;
  params.opt_tol = 1e-3;
  params.eta = 1/(N+1);
  params.max_iter = 1000;
  params.mu1fac = 10;
  params.mu1 = params.mu1fac*std(T(:));
  params.mu2 = params.mu1;
  params.mu_min = 1e-4;
  params.mu_max = 1e2;
  params.lambdaS = 1;
  params.lambda = params.lambdaS*r*params.rRatio; 
  params.verbose = 0;
  params.use_cont = true;
  params.k = [size(T,1) size(T,2) 1];
  %%%%%%%%%% for PROPACK %%%%%%%%%%%%
  % declare global var 'sv'
  global sv;
  global tmode;
  global use_propack;
  global curr_mu;
  sv = ceil(min(size(data.T)) * 0.1) * ones(1, N);
  use_propack = true;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  results = rpca_for_tensor(data, params);
  %results = tensor_rpca_tc_adal2(data, params);
  %if(strcmp(algorithm_id,'HoRPCA-IALM')) results = rpca_for_tensor(data, params); end
  %if(strcmp(algorithm_id,'HoRPCA-S')) results = tensor_rpca_adal2(data, params); end % tensor_rpca_adal
  %if(strcmp(algorithm_id,'HoRPCA-S-NCX')) results = tensor_rpca_adal_ncx(data, params); end
  %if(strcmp(algorithm_id,'Tucker-ADAL')) results = tensor_tucker_adal_ncx(data, params); end
  L = double(results.X);
  S = double(results.E);
  clear sv tmode use_propack curr_mu;
  T_hat = L;
end
