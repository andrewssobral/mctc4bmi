clear all
close all
clc

%% ====================== Load data ==============================
Video = {'basketball','windmill'};
ii = 2;
VideoName = [Video{ii} '_video.mat'];
load(VideoName);

num_chose = size(T,3);
num_chose = min(num_chose, 60);
T = T(:,:,1:num_chose);
normalize              =        max(T(:))                     ;
Xn                     =        T/normalize                   ;
[n1,n2,n3]             =        size(Xn)                      ;

p                      =        0.3                           ;
Omega                  =        zeros(size(Xn))               ;
chosen                 =        randperm(n1*n2*n3,...
                                       round(p*n1*n2*n3))     ;
Omega(chosen)          =        1                             ;

alpha                  =        1.05                             ;
maxItr                 =        200                          ; % maximum iteration

rho                    =        0.01                      ;

% rho                    =        0.01                          ;
% rho                    =        0.15;   



A                      =        diag(sparse(double(Omega(:)))); % sampling operator
b                      =        A * Xn(:)                     ; % available data
bb                     =        reshape(b,[n1,n2,n3]);

%% ================ main process of completion =======================
% X   =    tensor_cpl_admm_w( A , b , rho , alpha , ...
%                      [n1,n2,n3] , maxItr , myNorm , 0 );
[X,RSE] = LtSVD_TC(A ,b,rho,alpha ,[n1,n2,n3],maxItr, Xn(:));
X                      =        X * normalize                 ;
X                      =        reshape(X,[n1,n2,n3])         ;
%  [x,history] = WTNNM_TC(A,y,rho,alpha,sX,maxItr,QUIET,ImgShow)           
X_dif                  =        T-X                           ;
% RSE                    =        norm(X_dif(:))/norm(T(:))     ;
            
%% ======================== Result Plot =============================

figure('Name','Final_Results');
for i = 1:(size(X,3))
    subplot(221);imshow(Xn(:,:,i));title('Original');
    subplot(222);imshow(X(:,:,i)/normalize) ;title('Recovered');
    subplot(223);imshow(bb(:,:,i));title('Sampled ');
    subplot(224);imshow(X_dif(:,:,i)/normalize);title('Difference');
    pause(.1);
end
