function X=fold(X,n,D)

% Inverse operation of matrizicing, i.e. reconstructs the multi-way array
% from it's matriziced version.
%
% Input:
% D is vector containing original dimensions
% n is the dimension along which the matrizicing was originally performed
% X is a matrix
%
% Output:
% X (multi-way array)

if n==1
    perm=[1:length(D)];
else
    perm=[2:n 1 n+1:length(D)];
end
  
X=permute(reshape(X,D([n 1:n-1 n+1:length(D)])),perm);
