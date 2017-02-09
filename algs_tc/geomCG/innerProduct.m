function res = innerProduct( X, xi, nu )
%INNERPRODUCT Calculates inner product of two tangent tensors.
%   RES = INNERPRODUCT( X, XI, NU ) calculates the inner product of two
%   tangent tensors xi and nu given in the factored form. Both lie in the 
%   same tangent space to the point X on the manifold.
%   
%   See also addFactorized
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    R1 = ttm( X.core, xi.U1_tilde'*nu.U1_tilde, 1 );
    R2 = ttm( X.core, xi.U2_tilde'*nu.U2_tilde, 2 );
    R3 = ttm( X.core, xi.U3_tilde'*nu.U3_tilde, 3 );
    
    res = xi.Y_tilde(:)'*nu.Y_tilde(:) ...
            + X.core(:)'*R1(:) ...
            + X.core(:)'*R2(:) ...
            + X.core(:)'*R3(:);

end
