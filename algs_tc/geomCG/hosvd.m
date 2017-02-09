function T = hosvd( A, r )
%HOSVD Compute the HOSVD truncation
%   T = HOSVD(A, R) computes the rank-R HOSVD for the
%   3-dimensional tensor A. The result T is stored as a 
%   Tucker tensor (ttensor).
%
%   See also tenRank

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

A1 = double( tenmat(A,1) );
A2 = double( tenmat(A,2) );
A3 = double( tenmat(A,3) );

[U1,~,~] = svd( A1, 0 );
U1 = U1( :, 1:r(1) );
[U2,~,~] = svd( A2, 0 );
U2 = U2( :, 1:r(2) );
[U3,~,~] = svd( A3, 0 );
U3 = U3( :, 1:r(3) );

C = ttm( A, {U1,U2,U3}, 't' );
T = ttensor( C, {U1,U2,U3} );

end
