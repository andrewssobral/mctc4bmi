function perform_matrix_completion(params)
  %% Load pre-processed sequence
  current_seq_name = params.sequences_name(params.current_sequence).name;
  current_seq_path = fullfile(params.sequences_mat,strcat(current_seq_name,'.mat'));
  displog(['Loading sequence: ' current_seq_name]);
  load(current_seq_path,'V2','O2');
  %V2 = V2(:,:,:,2:end);
  %O2 = O2(:,:,2:end);
  % show_4dvideo(V2);
  % show_3dvideo(O2); slice3(O2);
  
  %%% Load algorithm
  current_alg_name = params.algs_name(params.current_algorithm).name;
  current_alg_path = fullfile(params.algs_path,current_alg_name);
  displog(['Loading algorithm: ' current_alg_name]);
  addpath(genpath(current_alg_path));
  
  %%% Matrix completion
  displog('Performing matrix completion');
  V2p = permute(V2,[1 2 4 3]);
  BG_hat = zeros(size(V2(:,:,:,1)));
  invO2 = double(~O2); % show_3dvideo(O2), show_3dvideo(invO2), slice3(invO2);
  Minv = convert_video3d_to_2d(invO2); % imagesc(Minv)
  Minv2 = ones(size(Minv)); % imagesc(Minv2)
  clear invO2;
  %%
  timerVal = tic;
  for i = 1:size(V2p,4)
    displog(['Processing channel: ' num2str(i)]);
    % 4D > 3D
    Vi = V2p(:,:,:,i); % show_3dvideo(Vi);
    % 3D > 2D
    Mi = convert_video3d_to_2d(Vi); % imagesc(Mi)
    % Zero values to 1e-6
    Mi(Mi == 0) = 1e-3;  % imagesc(Mi)
    % Multiply each channel by mask
    %Mi = Mi.*Minv; % imagesc(Minv)
    
    %%%%%%%%%%%%%%%%%%%%%%%% Matrix Completion %%%%%%%%%%%%%%%%%%%%%%%%
    Idx = Minv;
    %if(strcmp(current_alg_name,'GROUSE')) Idx = Minv2; end
    params_mc.M = Mi;
    params_mc.Idx = Idx;
    Mi_hat = run_mc(params_mc); % imagesc(Mi_hat)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Build background model for each channel
    Mi_bg = mean(Mi_hat,2); % imagesc(Mi_bg)
    %Mi_bg = median(Mi_hat,2);
    I_bg = reshape(Mi_bg,size(V2,1),size(V2,2)); % imagesc(I_bg)
    BG_hat(:,:,i) = I_bg;
  end
  elapsedTime = toc(timerVal); clear i timerVal;
  displog(['Elapsed time: ' num2str(elapsedTime)]);
  %%% BG model
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
