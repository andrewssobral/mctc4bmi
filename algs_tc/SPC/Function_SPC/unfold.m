function Y=unfold(X,n)

% The matrizicing operation, i.e.
% matrizicing(X,n) turns the multi-way array X into a matrix of dimensions:
% (size(X,n)) x (size(X,1)*...*size(X,n-1)*size(X,n+1)*...*size(X,N))
% where N=ndims(X)
%
% Input:
% X     Multi-way array
% n     Dimension along which matrizicing is performed
%
% Output:
% Y     The corresponding matriziced version of X

N=ndims(X);
Y=reshape(permute(X, [n 1:n-1 n+1:N]),size(X,n),prod(size(X))/size(X,n));
