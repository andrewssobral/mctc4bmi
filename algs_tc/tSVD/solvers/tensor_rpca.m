% ===================================================================
%            Tensor Robust PCA using ADMM
% ===================================================================
%
% We are trying to solve the following problem:
%       
%          min \|L\|_{TNN} + lambda \|S\|_{1,1,2}
%          s.t.      M = L + S
%
% The Lagrangian:
% L(L,S,Y) = \|L\|_TNN + lambda * \|S\|_1 + ...
%             trace(Y^T (L+S-M)) + 0.5*rho \|L+S-M\|_F^2
%
%          = \|L\|_TNN + lambda * \|S\|_1 + ...
%             0.5*rho \|L+S-M+W\|_F^2
%          where Y = rho W is the dual
%
% So the ADMM algorithm is as follows:
%
% L^{k+1} = argmin_L \|L\|_TNN + 0.5*rho \|L+S^{k}-M+W^{k}\|_F^2
%
% S^{k+1} = argmin_S lambda \|S\|_1 + 0.5*rho\|L^{k+1}+S-M+W^{k}\|_F^2
%
% W^{k+1} = W^{k} + rho*( L^{k+1}+S^{k+1}-M )
%
%  INPUTS:
%   
%   X      :   Observed data
%   lambda :   parameter as above
%
%  OUTPUTS:
%   L      :   Low rank component
%   S      :   Sparse component
%
% Authors: G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015
%
% ====================================================================

function [ L , S ] = tensor_rpca( Xn , lambda)
%%  ============= Parameters ====================
L_new          =        ones(size(Xn))                     ;
S_new          =        Xn-L_new                           ;
Y_new          =        ones(size(Xn))                     ;

rho            =       2/3                                 ;

L_gap          =       10                                  ;
S_gap          =       10                                  ;

tol_L          =       1e-3                                ;
tol_S          =       1e-3                                ;

iter           =       0                                   ;
max_iter       =       250                                 ;
s1='%3s\t%8s\t%10s\t%10s\t%10s\t%10s\n'                    ;
s2='%3d\t%10.1f\t%10.1f\t%10.4f\t%10.4f\t%10.1e\n'         ;

fprintf(s1,'ite','||L||_*','||S||_1','L_gap','S_gap','rho');
figure;
% ============= Main Iteration ====================
while L_gap >tol_L || S_gap>tol_S

    iter = iter+1                                          ; 
    if iter == max_iter
        break                                              ;
    end
    
    L_old         =     L_new                              ;
    S_old         =     S_new                              ;
    Y_old         =     Y_new                              ;
    
    term1         =     -( S_old-Xn+Y_old)                 ;
    [L_new,L_nuc] =     proxF_tSVD_1(term1,1/rho,[])       ;
    
    term2         =     -(L_new-Xn+Y_old)                  ;
    [S_new,S_l1]  =     proxF_tube_12(term2,lambda/rho)    ;   

    Y_new         =     Y_old + L_new+S_new-Xn             ;

    L_gap         =     norm(L_new(:)-L_old(:))            ;
    S_gap         =     norm(S_new(:)-S_old(:))            ;
    
    if mod(iter,2)== 0
        fprintf(s2,...
            iter,L_nuc,S_l1,L_gap,S_gap,rho)               ;
    end   

    subplot(121);imagesc(L_new(:,:,10));colormap(gray); drawnow;
    subplot(122);imagesc(S_new(:,:,10));colormap(gray); drawnow;
end
L = L_new;
S = S_new;

end

