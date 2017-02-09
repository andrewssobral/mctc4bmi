
% A demo of Bayesian CP factorization for image completion
% Written by Qibin Zhao 2014 RIKEN BSI
%
% In this demo, we provide two algorithms including BCPF_IC and BCPF_MP.
% BCPF_IC is a Bayesian CP for image completion; BCPF_MP is a Bayesian CP
% using mixture priors, which is particularly useful for natural image
% completion. For algorithm settings, please refer to the detailed help by 
% >> help BCPF_MP

% The experimental data can be tested with
% 1) Different image files
% 2) Observation rate (1-missing rate)
% The predictive image can be online visualized during model learning. 
% The performance of RSE, PSNR, SSIM, Time Cost are evaluated and reported.


close all; clear all;
randn('state',1); rand('state',1); %#ok<RAND>
%% Load image data
filename='./TestImages/peppers.bmp';    % Image file
ObsRatio = 0.1;                      % Observation rate

X = double(imread(filename));
DIM = size(X);

Omega = randperm(prod(DIM));
Omega = Omega(1:round(ObsRatio*prod(DIM)));
O = zeros(DIM);
O(Omega) = 1;
Y = O.*X;

% plot images
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.01 0.01]);
row =1; col =2;
figure;
subplot(row,col,1);
imshow(uint8(X));
subplot(row,col,2);
imshow(uint8(Y));
drawnow;

% Initialization
TimeCost = zeros(2,1);
RSElist = zeros(2,3);
PSNRlist = zeros(2,1);
SSIMlist = zeros(2,1);
RankEst = zeros(2,1);

%% BCPF for low-rank images completion
tStart = tic;
fprintf('------Bayesian CP factorization for Image Completion---------- \n');
[model] = BCPF_IC(Y, 'obs', O, 'init', 'rand', 'maxRank', 100, 'maxiters', 20, ...
    'tol', 1e-4, 'dimRed', 1, 'verbose', 2);
X_FBCP = double(model.X);
RSElist(1,1) = perfscore(X_FBCP, X);
RSElist(1,2) = perfscore(X_FBCP(O==1), X(O==1));
RSElist(1,3) = perfscore(X_FBCP(O==0), X(O==0));

X_FBCP(O==1) = X(O==1);
PSNRlist(1) = PSNR_RGB(X_FBCP,X);
SSIMlist(1) = ssim_index(rgb2gray(uint8(X_FBCP)),rgb2gray(uint8(X)));
RankEst(1) = model.TrueRank;
TimeCost(1) = toc(tStart);
% figure; imshow(uint8(X_FBCP)); title('FBCP','FontWeight','bold'); drawnow;


%% BCPF-MP (mixture priors) for natural images
if ~isempty(strfind(filename,'facade.bmp'))
    nd=0.1;    % low-rank structural images
else
    nd=1;      % natural images
end
tStart = tic;
fprintf('------Bayesian CP with Mixture Priors for Image Completion---------- \n');
[model] = BCPF_MP(Y, 'obs', O, 'init', 'rand', 'maxRank', 100, 'maxiters', 30, ...
    'tol', 1e-4, 'dimRed', 1, 'verbose', 2, 'nd', nd);
X_FBCPS = double(model.X);

RSElist(2,1) = perfscore(X_FBCPS, X);
RSElist(2,2) = perfscore(X_FBCPS(O==1), X(O==1));
RSElist(2,3) = perfscore(X_FBCPS(O==0), X(O==0));

X_FBCPS(O==1) = X(O==1);
PSNRlist(2) = PSNR_RGB(X_FBCPS,X);
SSIMlist(2) = ssim_index(rgb2gray(uint8(X_FBCPS)),rgb2gray(uint8(X)));
RankEst(2) = model.TrueRank;
TimeCost(2) = toc(tStart);
% figure; imshow(uint8(X_FBCPS)); title('FBCP-MP','FontWeight','bold'); drawnow;

%%
RankEst
RSElist
PSNRlist
SSIMlist
TimeCost



