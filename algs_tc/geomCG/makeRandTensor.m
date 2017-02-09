function X = makeRandTensor( n, k )
%MAKERANDTENSOR Create a random Tucker tensor
%   X = MAKERANDTENSOR( N, K ) creates a random Tucker tensor X stored as a 
%   ttensor object. The entries of the core tensor the basis factors are chosen
%   independently from the uniform distribution on [0,1]. Finally, the basis factors
%   are orthogonalized using a QR procedure.
%
%   See also makeOmegaSet
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    [U1,R1] = qr( rand( n(1), k(1) ), 0);
    [U2,R2] = qr( rand( n(2), k(2) ), 0);
    [U3,R3] = qr( rand( n(3), k(3) ), 0);

    C  = tenrand( k );
    C = ttm( C, {R1,R2,R3},[1,2,3]);

    X = ttensor( C, {U1, U2, U3} );
end
