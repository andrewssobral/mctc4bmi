% A demo of Bayesian CP Factorization on synthetic data
% Written by  Qibin Zhao 2013 RIKEN BSI
%
% In this demo, we provide two Bayesian CP Factorization algorithms: one for
% incomplete tensor and tensor completion ("BCPF_TC.m") and the other one for
% complete tensor, which is a more efficient implementation when tensor is 
% fully observed ("BCPF.m").  
% The parameter settings are optional, please refer to the help by command 
% >> help BCPF_TC  
% Two important settings are initialization method and initilization of maximum rank.
% 
% This demo is used for evaluation of CP Rank determination and predictive
% perofrmance under varying conditions including Tensor size, True rank, 
% Signal type, Observation ratio, and Noise SNRs.   

% After the model learning, the results can be shown by
% 1. Visualization of true latent factors and the estimated factors
% 2. Visualization of observed tensor Y,  the estimated low-CP-rank tensor X


close all;
randn('state',1); rand('state',1); %#ok<RAND>
%% Generate a low-rank tensor
DIM = [30,30,30];     % Dimensions of data
R = 3;                % True CP rank
DataType = 2;         % 1: Random factors   2: The deterministic factors (sin, cos, square)

Z = cell(length(DIM),1);   
if DataType ==1
    for m=1:length(DIM)
          Z{m} =  gaussSample(zeros(R,1), eye(R), DIM(m));  
    end
end
if DataType == 2
    for m=1:length(DIM)
        temp = linspace(0, m*2*pi, DIM(m));
        part1 = [sin(temp); cos(temp); square(linspace(0, 16*pi, DIM(m)))]';
        part2 = gaussSample(zeros(DIM(m),1), eye(DIM(m)), R-size(part1,2))';
        Z{m} = [part1 part2];
        Z{m} = Z{m}(:,1:R);
        Z{m} = zscore(Z{m});
    end
end
% Generate tensor by factor matrices
X = double(ktensor(Z));

%% Random missing values
ObsRatio = 0.2;               % Observation rate: [0 ~ 1]
Omega = randperm(prod(DIM)); 
Omega = Omega(1:round(ObsRatio*prod(DIM)));
O = zeros(DIM); 
O(Omega) = 1;

%% Add noise
SNR = 10;                     % Noise levels
sigma2 = var(X(:))*(1/(10^(SNR/10)));
GN = sqrt(sigma2)*randn(DIM);

%% Generate observation tensor Y
Y = X + GN;
Y = O.*Y;

%% Run BayesCP
fprintf('------Bayesian CP Factorization---------- \n');
ts = tic;
if ObsRatio~=1 
    % Bayes CP algorithm for incomplete tensor and tensor completion    
    [model] = BCPF_TC(Y, 'obs', O, 'init', 'ml', 'maxRank', max([DIM 2*R]), 'dimRed', 1, 'tol', 1e-6, 'maxiters', 100, 'verbose', 2);
else
    % Bayes CP algorithm for fully observed tensor 
    [model] = BCPF(Y, 'init', 'ml', 'maxRank', max([DIM 2*R]), 'dimRed', 1, 'tol', 1e-6, 'maxiters', 200, 'verbose', 2);
end
t_total = toc(ts);

% Performance evaluation
X_hat = double(model.X);
err = X_hat(:) - X(:);
rmse = sqrt(mean(err.^2));
rrse = sqrt(sum(err.^2)/sum(X(:).^2));

% Report results
fprintf('\n------------Bayesian CP Factorization-----------------------------------------------------------------------------------\n')
fprintf('Observation ratio = %g, SNR = %g, True Rank=%d\n', ObsRatio, SNR, R);
fprintf('RRSE = %g, RMSE = %g, Estimated Rank = %d, \nEstimated SNR = %g, Time = %g\n', ...
    rrse, rmse, model.TrueRank, model.SNR, t_total);
fprintf('--------------------------------------------------------------------------------------------------------------------------\n')

%% Visualization of data and results
plotYXS(Y, X_hat);
factorCorr = plotFactor(Z,model.X.U);



