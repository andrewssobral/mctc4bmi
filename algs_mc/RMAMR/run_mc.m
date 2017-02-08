function M_hat = run_mc(params)
  D = params.M;
  CM = params.Idx;
  algorithm_id = 'RMAMR';
  
  L = zeros(size(D));
  S = zeros(size(D));
  E = zeros(size(D));

  %CM = ones(size(D));
  %CM = zeros(size(D));
  [numr,numc] = size(D);
  %CM = randi([0 1],numr,numc); % simulated confidence map (binary matrix)

  Omega = find(CM ~= 0);
  [I,J] = ind2sub([numr numc],Omega);  

  % Motion-Assisted Matrix Restoration (Ye et al. 2015)
  if(strcmp(algorithm_id,'MAMR'))
    lambda = 1; 
    [L,S,iter1] = core_MAMR(D,lambda,I,J);
  end

  % Robust Motion-Assisted Matrix Restoration (Ye et al. 2015)
  if(strcmp(algorithm_id,'RMAMR'))
    lambda = 10; 
    [L,S,E,iter1] = core_RMAMR(D,lambda,I,J);
  end

  M_hat = L;
  
  %M_hat = L + S + E;
  %error = norm(M_hat(:)-M(:))/norm(M(:));
  %disp(['Error: ' num2str(error)]);
end
