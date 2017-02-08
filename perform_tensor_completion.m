function perform_tensor_completion(params)
  %% Load pre-processed sequence
  current_seq_name = params.sequences_name(params.current_sequence).name;
  current_seq_path = fullfile(params.sequences_mat,strcat(current_seq_name,'.mat'));
  displog(['Loading sequence: ' current_seq_name]);
  load(current_seq_path,'V2','O2');
  % show_4dvideo(V2);
  % show_3dvideo(O2); slice3(O2);
  
  %%% Load algorithm
  current_alg_name = params.algs_name(params.current_algorithm).name;
  current_alg_path = fullfile(pwd,params.algs_path,current_alg_name);
  displog(['Loading algorithm: ' current_alg_name]);
  addpath(genpath(current_alg_path));
  
  %%% Block size for subtensors
  params_tc = [];
  switch(current_seq_name)
    case{'Board'}
      params_tc.blksize = [4 8];
    case {'CaVignal','Foliage','Snellen'}
      params_tc.blksize = [8 8];
    case{'Toscana'}
      params_tc.blksize = [25 25];
    otherwise
      params_tc.blksize = [16 16];
  end
  
  %% Tensor completion
  clc;
  displog('Performing tensor completion');
  T = permute(V2,[1 2 4 3]); % width*height*time*channels
  Idx = double(~O2); % show_3dvideo(Idx)
  timerVal = tic;
  BG_hat = [];
  for k = 1:size(T,4)
    displog(['Processing channel: ' num2str(k)]);
    % 4D > 3D
    Tk = T(:,:,:,k); % show_3dvideo(Tk); slice3(Tk); colormap('gray');
    % imagesc(convert_video3d_to_2d(Tk_hat)); colormap('gray');
    % show_3dvideo(Tk.*Idx); slice3(Tk.*Idx);
    
    %%%%%%%%%%%%%%%%%%%%%%%% Tensor Completion %%%%%%%%%%%%%%%%%%%%%%%%
    params_tc.T = Tk;
    params_tc.Idx = Idx;  % show_3dvideo(Idx); slice3(Idx);
    Tk_hat = run_tc(params_tc); % show_3dvideo(Tk_hat); slice3(Tk_hat);
    %Tk_hat = zeros(size(Tk));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Build background model for each channel
    I_bg = mean(Tk_hat,3);
    %I_bg = median(Tk_hat,2);
    
    BG_hat(:,:,k) = I_bg;
    
    % debug
    %clf,imagesc(I_bg); break;
  end
  elapsedTime = toc(timerVal);
  disp(['Elapsed time: ' num2str(elapsedTime)]);
  %% BG model
  if(params.debug)
    clf,imshow(BG_hat);
    pause(1);
  end
  %% Save results
  results_path = fullfile(params.results_path,current_alg_name);
  if(~exist(results_path,'dir'))
    mkdir(results_path);
  end
  result_file = fullfile(results_path,strcat('BC_',current_seq_name,'.',params.sequences_format));
  displog(['Saving results at: ' result_file]);
  imwrite(BG_hat,result_file);
  %% End
  rmpath(genpath(current_alg_path));
  clear variables;
end

