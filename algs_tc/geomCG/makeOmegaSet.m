function subs = makeOmegaSet( n, sizeOmega )
%MAKEOMEGASET Create a sampling set 
%   SUBS = MAKEOMEGASET( N, SIZEOMEGA ) creates a sampling set of SIZEOMEGA entries chosen uniformly from
%   the total amount of PROD( N ) entries.
%
%   See also makeRandTensor
%
    
%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    if sizeOmega > prod(n)
        error('makeOmegaSet:sizeOmegaTooHigh', 'Requested size of Omega is bigger than the tensor itself!')
    end

    idx = randi( prod(n), sizeOmega, 1 );
    Omega = unique(idx);

    while length(Omega) < sizeOmega
        idx = [ Omega; randi( prod(n) , sizeOmega-length(Omega), 1 )]; 
        Omega = unique(idx);
    end
    
    Omega = sort( Omega(1:sizeOmega) );

    [i,j,k] = ind2sub( n, Omega );

    subs = [i,j,k];
end
