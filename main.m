%%% Setup
% restoredefaultpath;
close, clear, clc;

%%% Path
addpath(fullfile(pwd,'utils'));

%%% Params
params.sequences_path = fullfile(pwd,'sequences');
params.sequences_mat = fullfile(pwd,'sequences_mat');
params.sequences_name = get_sequences_name(params);
params.sequences_format = 'PNG';
params.sequences_ext = strcat('*.',params.sequences_format);
params.debug = 0;
%params.current_sequence = 1;
clc;

%% Frame selection
displog('Performing frame selection');
for seq = 1:size(params.sequences_name,2)
    %seq = 1;
    displog(['Processing sequence: ' params.sequences_name(seq).name]);
    params.current_sequence = seq;
    perform_frame_selection(params);
    %break;
end
clear seq;

%% Matrix Completion
clc;
params.algs_path = 'algs_mc';
params.algs_name = get_algs_name(params);
params.results_path = 'results_mc';
displog('--- Matrix Completion ---');
for alg = 1:size(params.algs_name,2)
  %alg = 1;
  displog(['Current algorithm: ' params.algs_name(alg).name]);
  params.current_algorithm = alg;
  for seq = 1:size(params.sequences_name,2)
    %seq = 1;
    displog(['Processing sequence: ' params.sequences_name(seq).name]);
    params.current_sequence = seq;
    perform_matrix_completion(params);
    %break;
  end
  %break;
end
clear alg seq;

%% Tensor Completion
clc;
params.algs_path = 'algs_tc';
params.algs_name = get_algs_name(params);
params.results_path = 'results_tc';
displog('--- Tensor Completion ---');
for alg = 1:size(params.algs_name,2)
  %alg = 1;
  displog(['Current algorithm: ' params.algs_name(alg).name]);
  params.current_algorithm = alg;
  for seq = 1:size(params.sequences_name,2)
    %seq = 1;
    displog(['Processing sequence: ' params.sequences_name(seq).name]);
    params.current_sequence = seq;
    perform_tensor_completion(params);
    %break;
  end
  %break;
end
clear alg seq;
