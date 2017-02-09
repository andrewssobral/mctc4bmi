function C=tprod(A,B)
% Computes the tensor-product of two P dimensional tensors (A & B)
% All dimensions of the inputed tensor except for the 2nd must be the same.
% The algorithm follows from paper by ...
%
% INPUTS:
%
% A - n1 x n2 x n3 .... nP tensor
% B - n2 x n4 x n3 .... nP tensor
%
% OUTPUTS: 
%
% C - n1 x n4 x n3 .... nP tensor
%
% Original author :  Misha Kilmer, Ning Hao
% Edited by       :  G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015


sa = size(A);
sb = size(B);

faces = length(sa);

nfaces = prod(sb(3:end));

for i = 3:faces
    A = fft(A,[],i);
    B = fft(B,[],i);
end
sc = sa;
sc(2) = sb(2);
C = zeros(sc);
for i = 1:nfaces
    C(:,:,i) = A(:,:,i)*B(:,:,i);
end
for i = 3:faces
    C = ifft(C,[],i);
end



