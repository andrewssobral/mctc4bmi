clear all; close all;

addpath Function_SPC

%% load data (you can change the file addres in this sentence.)
%addres = 'Dataset/Original/lena.bmp';
addres = 'Dataset/lena_missing_80.png';
X0 = imread(addres);
II = size(X0);

%sig = 10.0;
%noise = sig*randn(II);
%X0 = double(X0) + noise;

%% set missing indexes
Q = (X0 ~= 255);
%Q = true(II);
T = zeros(II);
T(Q) = double(X0(Q));

%T = X0;

%% hyperparameters and run SPC

%TVQV    = 'tv';        % 'tv' or 'qv' ;
%rho     = [1 1 0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.

TVQV    = 'qv';        % 'tv' or 'qv' ;
rho     = [0.5 0.5 0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
K       = 10;          % Number of components which are updated in one iteration. (typically 10)
SNR     = 25;          % error bound
nu      = 0.01;        % threshold for R <-- R + 1.
maxiter = 1000;        % maximum number of iteration
tol     = 1e-5;        % tolerance
out_im  = 1;           % you can monitor the process of 'image' completion if out == 1. 'saved' directory is necessary to save the individual rank images.

[X Z G U histo histo_R] = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);
%[X Z G U histo histo_R] = SPC_beta(T,Q,'sv',rho,K,SNR,nu,maxiter,tol,out_im);
%[X Z G U histo histo_R] = SPC_old(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);
imwrite(uint8(X),'result.png');

%% visualize distribution of G and 40 parts of rank-one tensors
[val idd] = sort(G,'descend');
UU = [];
for ii = 1:length(idd)

  Z1 = outerprod(U,idd(ii));
  Z1 = Z1/max(abs(Z1(:)))*256;
  imwrite(uint8(Z1),['saved_parts/' num2str(ii) '.png']);
  UU = [UU Z1(:)];

end
figure;
subplot(2,1,1)
plot(val,'linewidth',2);
subplot(2,1,2)
P = reshape(uint8(UU(:,1:40)),[256 256 3 10 4]);
P = permute(P,[1 5 2 4 3]);
P = reshape(P,[256*4 256*10 3]);
imagesc(P)
print -depsc parts_SPC.eps


%%%%%%%%%%%%%% test for STADARD PARAFAC
addpath('/home/yokota/toolboxes/tensorlab/');
%Tl = double(T);
%Tl(~Q) = NaN;
addres = 'Dataset/Original/lena.bmp';
X0 = imread(addres);
Tl = double(X0);
Tl = fmt(Tl);

Ufr{1} = randn(256,139);
Ufr{2} = randn(256,139);
Ufr{3} = randn(3,139);
options.Display = true;
sol = cpd_als(Tl,Ufr,options);
Xfr = cpdgen(sol);
imagesc(uint8(Xfr));
%imwrite(uint8(Xfr),'result_CPD.png');

for ii = 1:139
  Gfr(ii) = norm(sol{1}(:,ii))*norm(sol{2}(:,ii))*norm(sol{3}(:,ii));
end


