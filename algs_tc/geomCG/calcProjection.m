function xi = calcProjection( X, Y )
%CALCPROJECTION Projection into the tangent space.
%   XI = CALCPROJECTION( X, Y ) computes the orth. projection of the Tucker 
%   tensor (ttensor) Y into the tangent space at X (ttensor).
%   Of the resulting tangent tensor XI, only the variations are stored.
%   Hence, XI is a struct with the fields
%
%       XI.Y_tilde  (var. in the core)
%       XI.U1_tilde (var. in the factors...
%       XI.U2_tilde
%       XI.U3_tilde  ...)
%   
%   For efficiency, parts of the calculation are performed by the mex
%   routine calcProjection_mex
%
%   See also calcGradient.
%
    
%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    [temp,temp1,temp2,temp3] = calcProjection_mex( Y.subs', Y.vals, X.U{1}', X.U{2}', X.U{3}' );
    
    [n1,k1] = size(X.U{1});
    [n2,k2] = size(X.U{2});
    [n3,k3] = size(X.U{3});

    xi.Y_tilde = tensor(reshape( temp, [k1 k2 k3]), [k1, k2, k3]);

    Y23 = tensor(reshape(temp1, [n1 k2 k3]));
    Y13 = tensor(reshape(temp2, [k1 n2 k3]));
    Y12 = tensor(reshape(temp3, [k1 k2 n3]));

    S1_inv = pinv( double( tenmat( X.core, 1 ) ));
    S2_inv = pinv( double( tenmat( X.core, 2 ) ));
    S3_inv = pinv( double( tenmat( X.core, 3 ) ));

    U1_tilde = double( tenmat( Y23, 1 )) * S1_inv; 
    U2_tilde = double( tenmat( Y13, 2 )) * S2_inv; 
    U3_tilde = double( tenmat( Y12, 3 )) * S3_inv; 

    xi.U1_tilde = U1_tilde - X.U{1} * ( X.U{1}' * U1_tilde ); 
    xi.U2_tilde = U2_tilde - X.U{2} * ( X.U{2}' * U2_tilde ); 
    xi.U3_tilde = U3_tilde - X.U{3} * ( X.U{3}' * U3_tilde ); 

end

