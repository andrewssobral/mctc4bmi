function t = calcInitial( Deriv, X, xi) 
%CALCINITIAL Calculate the initial guess for the line search.
%
%   Wrapper function for calcInitial_mex.c
%
%   See also calcGradient, calcFunction   
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    t = -calcInitial_mex( Deriv.subs', Deriv.vals, ...
                          X.core.data, X.U{1}', X.U{2}', X.U{3}',...
                          xi.Y_tilde.data, xi.U1_tilde', xi.U2_tilde', xi.U3_tilde');

end
