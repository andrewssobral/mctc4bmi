function res = calcGradient( A_Omega, X )
%CALCGRADIENT Calculate the euclid. gradient of the obj. function
%   Wrapper function for calcGradient_mex.c
%   
%   Computes the euclid. gradient of the objective function
%   
%         X_Omega - A_Omega
%   
%   between a sparse tensor A_Omega and a Tucker tensor X.

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    vals = calcGradient_mex( A_Omega.subs', A_Omega.vals, ...
                            X.core.data, X.U{1}', X.U{2}', X.U{3}');

    res = sptensor( A_Omega.subs, vals, A_Omega.size, 0);

end
