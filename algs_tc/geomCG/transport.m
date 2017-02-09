function nu = transport( X, Y, xi )
%TRANSPORT Calculates vector transport
%   NU = TRANSPORT( X, Y, XI ) calculates the vector transport from
%   
%         T_X --> X_Y
%   
%   where X and Y are given in Tucker format (ttensor) and applies it 
%   to the tangent vector XI. XI is decomposed into
%     XI = Y_tilde x U1 x U2 x U3 
%          + S x U1_tilde x U2 x U3 
%          + S x U1 x U2_tilde x U3 
%          + S x U1 x U2 x U3_tilde 
%   
%   where only Y_tilde, U1_tilde, U2_tilde and U3_tilde are stored
%   as a struct in factorized form.
%
%   See also calcProjection, addFactorized
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    V1 = Y.U{1}'*X.U{1};
    V2 = Y.U{2}'*X.U{2};
    V3 = Y.U{3}'*X.U{3};

    V1_tilde = Y.U{1}'*xi.U1_tilde;
    V2_tilde = Y.U{2}'*xi.U2_tilde;
    V3_tilde = Y.U{3}'*xi.U3_tilde;

    M1 = ttm( xi.Y_tilde, {V1, V2, V3} );
    M2 = ttm( X.core, {V1_tilde, V2, V3 } );
    M3 = ttm( X.core, {V1, V2_tilde, V3 } );
    M4 = ttm( X.core, {V1, V2, V3_tilde } );

    %first part:
    nu.Y_tilde = M1 + M2 + M3 + M4;

    %second part;
        
    Y1 = ttm(xi.Y_tilde, {X.U{1}, V2, V3} );
    Y2 = ttm(xi.Y_tilde, {V1, X.U{2}, V3} );
    Y3 = ttm(xi.Y_tilde, {V1, V2, X.U{3}} );

    G1_1 = ttm( X.core, {xi.U1_tilde,   V2,         V3} );
    G1_2 = ttm( X.core, {X.U{1},        V2_tilde,   V3} );
    G1_3 = ttm( X.core, {X.U{1},        V2,         V3_tilde} );

    G2_1 = ttm( X.core, {V1_tilde,  X.U{2},         V3} );
    G2_2 = ttm( X.core, {V1,        xi.U2_tilde,    V3} );
    G2_3 = ttm( X.core, {V1,        X.U{2},         V3_tilde} );

    G3_1 = ttm( X.core, {V1_tilde,  V2,         X.U{3}} );
    G3_2 = ttm( X.core, {V1,        V2_tilde,   X.U{3}} );
    G3_3 = ttm( X.core, {V1,        V2,         xi.U3_tilde} );

    S1_inv = pinv( double( tenmat( Y.core, 1 ) ));
    S2_inv = pinv( double( tenmat( Y.core, 2 ) ));
    S3_inv = pinv( double( tenmat( Y.core, 3 ) ));

    U1_tilde = double( tenmat(Y1 + G1_1 + G1_2 + G1_3, 1) ) * S1_inv;
    U2_tilde = double( tenmat(Y2 + G2_1 + G2_2 + G2_3, 2) ) * S2_inv;
    U3_tilde = double( tenmat(Y3 + G3_1 + G3_2 + G3_3, 3) ) * S3_inv;

    nu.U1_tilde = U1_tilde - Y.U{1} * ( Y.U{1}' * U1_tilde );
    nu.U2_tilde = U2_tilde - Y.U{2} * ( Y.U{2}' * U2_tilde );
    nu.U3_tilde = U3_tilde - Y.U{3} * ( Y.U{3}' * U3_tilde );

end





