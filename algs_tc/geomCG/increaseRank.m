function X = increaseRank( X, k )
%INCREASERANK Rank increase procedure for rank adaptation heuristic
%   X = INCREASERANK(X, K) moves the ttensor X from the old rank-R
%   manifold to the new rank-(R+K) manifold. Here, K is the vector of
%   rank increments.
%
%   See also tenRank
%
    
%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    n = size(X);
    m = size(X.core);
    [U1, R1] = qr( [X.U{1}, rand( n(1), k(1))], 0 );
    r1 = diag(R1);
    r1 = [ r1( 1:m(1) ); ones( k(1), 1 ) ];
    U1 = U1 * diag(r1);

    [U2, R2] = qr( [X.U{2}, rand( n(2), k(2))], 0 );
    r2 = diag(R2);
    r2 = [ r2( 1:m(2) ); ones( k(2), 1 ) ];
    U2 = U2 * diag(r2);

    [U3, R3] = qr( [X.U{3}, rand( n(3), k(3))], 0 );
    r3 = diag(R3);
    r3 = [ r3( 1:m(3) ); ones( k(3), 1 ) ];
    U3 = U3 * diag(r3);

    A = eps*tenones( m + k );
    A( 1:m(1), 1:m(2), 1:m(3) ) = X.core; 

    X = ttensor( A, {U1, U2, U3} );

end
