function res = retraction( X, xi, t )
%RETRACTION Retraction back to the manifold
%   RES = RETRACTION( X, XI, T) retracts an element of the tangent space
%   at the point X of the manifold, (X + T*XI), back to the manifold.  
%   XI is a tangent vector given in the factorized form
%   T is the step length.
%   
%   See also calcProjection, calcNewDirection, addFactorized
%   

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    [Q1, R1] = qr( [ X.U{1}, xi.U1_tilde ], 0 );
    [Q2, R2] = qr( [ X.U{2}, xi.U2_tilde ], 0 );
    [Q3, R3] = qr( [ X.U{3}, xi.U3_tilde ], 0 );
    
    k = size(X.core);

    S = tenzeros( 2*k );

    S(1:k(1), 1:k(2), 1:k(3)) = X.core + t*xi.Y_tilde;
    S(k(1)+1:end, 1:k(2), 1:k(3)) = t*X.core;
    S(1:k(1), k(2)+1:end, 1:k(3)) = t*X.core;
    S(1:k(1), 1:k(2), k(3)+1:end) = t*X.core;

    S = ttm( S, {R1, R2, R3});
    
    HO = hosvd(S, k);

    Q1 = Q1 * HO.U{1};
    Q2 = Q2 * HO.U{2};
    Q3 = Q3 * HO.U{3};

    res = ttensor( HO.core, {Q1, Q2, Q3});

end
