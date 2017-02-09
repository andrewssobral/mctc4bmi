%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tensor completion for esimtating missing values in visual data, PAMI,
% example.m file
% Date Mar. 13, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%addpath('mylib/'); % rmpath('mylib/');

T = double(imread('testImg.png')); % imshow(uint8(T))
%T = imresize(T,0.5);
Omega = (T < 254); % imagesc(Omega)

alpha = [1, 1, 1e-3];
alpha = alpha / sum(alpha);

maxIter = 500;
epsilon = 1e-5;

% "X" returns the estimation, 
% "errList" returns the list of the relative difference of outputs of two neighbor iterations 

%% High Accuracy LRTC (solve the original problem, HaLRTC algorithm in the paper)
rho = 1e-6;
[X_H, errList_H] = HaLRTC(...
    T, ...                       % a tensor whose elements in Omega are used for estimating missing value
    Omega,...               % the index set indicating the obeserved elements
    alpha,...                  % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
    rho,...                      % the initial value of the parameter; it should be small enough  
    maxIter,...               % the maximum iterations
    epsilon...                 % the tolerance of the relative difference of outputs of two neighbor iterations 
    );

%% Fast LRTC (solve the relaxed formulation, FaLRTC algorithm in the paper)
mu = 5 * alpha ./ sqrt(size(T));
C =  0.6;
L0 = 1e-5;
[X_F, errList_F] = FaLRTC(...
     T,...                    % a tensor whose elements in Omega are used for estimating missing value
     Omega,...          % the index set indicating the obeserved elements
     alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     mu,...                 % the relaxation parameters, mu >= 0. See the function for definitions
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon...          % the tolerance of the relative difference of outputs of two neighbor iterations 
     );

%% Fast LRTC no relaxation version (solve the original formuation, not appear in the paper but worthing to try it)
factor = 2;
C =  0.6;
L0 = 1e-5;
[X_Fn, errList_Fn] = FaLRTCnr(...
     T,...                    % a tensor whose elements in Omega are used for estimating missing value
     Omega,...          % the index set indicating the obeserved elements
     alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     factor,...             % control the decay speed of the relaxation parameter \mu, i.e., increase \mu by \mu = O(1 / k^factor) where k denotes the kth iteration. A reasonable "factor" could be in the range (0, 4].  
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon...          % the tolerance of the relative difference of outputs of two neighbor iterations 
    );

%% Simple LRTC (solve the relaxed formulation, SiLRTC in the paper)
beta = 0.1*ones(1, ndims(T));
[X_S, errList_S] = SiLRTC(...
    T,...                      % a tensor whose elements in Omega are used for estimating missing value
    Omega,...           % the index set indicating the obeserved elements
    alpha, ...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
    beta,...                % the relaxation parameter. The larger, the closer to the original problem. See the function for definitions.
    maxIter,...         % the maximum iterations
    epsilon...            % the tolerance of the relative difference of outputs of two neighbor iterations 
    );

%% Simple LRTC no relaxation version (solve the original problem, not appear in the paper but worthing to try it)
factor = 0.98;
 [X_Sn, errList_Sn] = SiLRTCnr(...
    T,...                      % a tensor whose elements in Omega are used for estimating missing value
    Omega,...           % the index set indicating the obeserved elements
    alpha, ...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
    factor,...             % control the increasing speed of \gamma, i.e., \gamma = \gamma * factor 
    maxIter,...         % the maximum iterations
    epsilon...            % the tolerance of the relative difference of outputs of two neighbor iterations 
     );

figure('units','normalized','outerposition',[0 0 1 1]);
% figure;
% set(0, 'DefaultFigurePosition', [0, 0, 2000, 1000]);
subplot(2,4,1);
imshow(uint8(T));
title('Original');
subplot(2,4,2);
imshow(double(1-Omega));
title('Missing Values');
subplot(2,4,3);
imshow(uint8(X_H));
title('HaLRTC');
subplot(2,4,4);
imshow(uint8(X_F));
title('FaLRTC');
subplot(2,4,5);
imshow(uint8(X_S));
title('SiLRTC');
subplot(2,4,6);
imshow(uint8(X_Fn));
title('FaLRTC without relaxation');
subplot(2,4,7);
imshow(uint8(X_Sn));
title('SiLRTC without relaxation');
h = subplot(2,4,8);
plot(1:length(errList_H), -log(errList_H), ':r', 'linewidth', 1.5); hold on;
plot(1:length(errList_F), -log(errList_F), '-b', 'linewidth', 1.5); hold on;
plot(1:length(errList_Fn), -log(errList_Fn), '--g', 'linewidth', 1.5); hold on;
plot(1:length(errList_S), -log(errList_S), '-.y', 'linewidth', 1.5); hold on;
plot(1:length(errList_Sn), -log(errList_Sn), ':k', 'linewidth', 1.5); hold on;
title('Convergence Curves');
ylabel('-log(||x_{k+1}-x_{k}||/||T||)');
xlabel('iterations');
legend('HaLRTC', 'FaLRTC', 'FaLRTC without relaxation', 'SiLRTC', 'SiLRTC without relaxation', 4);



