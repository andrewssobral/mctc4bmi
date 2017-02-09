function r = tenRank( A, tol )
%TENRANK Calculates the multilinear rank.
%   R = TENRANK( A, TOL ) calculates the multilinear rank R of the
%   given tensor A. All singular values below the specified tolerance
%   TOL are considered to be zero.
%   
%   R = TENRANK( A ) calculates the multilinear rank R of the given 
%   tensor A. The default tolerance 1e-5 is chosen for the cutoff.
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

   
    if ~exist('tol','var')
        tol = 1e-5;
    end

    r = zeros( 1, ndims( A ) );

    for i = 1:ndims( A )
        r(i) = rank( double( tenmat( A, i ) ), tol);
    end
end
