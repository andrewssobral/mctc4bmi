function res = addTT( A, B )
%ADDTT Add two Tucker tensors in the ttensor class format
%   Z = ADDTT(X, Y) adds two ttensors X and Y and returns the
%   result.
%
%   See also addFactorized

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    resU1 = [A.U{1}, B.U{1}];
    resU2 = [A.U{2}, B.U{2}];
    resU3 = [A.U{3}, B.U{3}];

    As = A.core.size;
    Bs = B.core.size;
    resCore = tenzeros( A.core.size + B.core.size );
    resCore( 1:As(1), 1:As(2), 1:As(3) ) = A.core;
    resCore( As(1)+1:end, As(2)+1:end, As(3)+1:end ) = B.core;

    res = ttensor( resCore, {resU1, resU2, resU3} );

end
