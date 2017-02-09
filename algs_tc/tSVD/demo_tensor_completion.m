clear all
close all
clc

%% ====================== Load data ==============================
addpath('tSVD','proxFunctions','solvers')      ;
load video.mat
normalize              =        max(T(:))                     ;
Xn                     =        T/normalize                   ;
[n1,n2,n3]             =        size(Xn)                      ;

p                      =        0.5                           ;
Omega                  =        zeros(size(Xn))               ;
chosen                 =        randperm(n1*n2*n3,...
                                       round(p*n1*n2*n3))     ;
Omega(chosen)          =        1                             ;

alpha                  =        1                             ;
maxItr                 =        1000                          ; % maximum iteration
rho                    =        0.01                          ;
myNorm                 =        'tSVD_1'                      ; % dont change for now

A                      =        diag(sparse(double(Omega(:)))); % sampling operator
b                      =        A * Xn(:)                     ; % available data
bb                     =        reshape(b,[n1,n2,n3]);

%% ================ main process of completion =======================
X   =    tensor_cpl_admm( A , b , rho , alpha , ...
                     [n1,n2,n3] , maxItr , myNorm , 0 );
X                      =        X * normalize                 ;
X                      =        reshape(X,[n1,n2,n3])         ;
            
X_dif                  =        T-X                           ;
RSE                    =        norm(X_dif(:))/norm(T(:))     ;
            
%% ======================== Result Plot =============================

figure;
for i = 1:(size(X,3))
    subplot(221);imagesc(Xn(:,:,i));axis off;
    colormap(gray);title('Original Video');
    subplot(222);imagesc(X(:,:,i)) ;axis off;
    colormap(gray);title('Sampled Video');
    subplot(224);imagesc(bb(:,:,i));axis off;
    colormap(gray);title('Recovered Video');
    pause(.1);
end
