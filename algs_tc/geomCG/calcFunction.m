function res = calcFunction( A_Omega, X )
%CALCFUNCTION Calculate the value of the objective function
%   Wrapper function for calcFunction_mex.c
%
%   Computes the value of the objective Function
%
%       0.5 * || X_Omega - A_Omega ||^2
%
%   See also calcGradient

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    res = 0.5 * calcFunction_mex( A_Omega.subs', A_Omega.vals, ...
                            X.core.data, X.U{1}', X.U{2}', X.U{3}');
end
