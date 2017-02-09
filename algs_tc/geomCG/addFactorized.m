function res = addFactorized( X, Y, alpha)
%ADDFACTORIZED Addition of two tangent tensors
%   Z = ADDFACTORIZED(X, Y, ALPHA) adds two tangent tensors in the same
%   tangent space, where Y is scaled by the factor ALPHA. 
%       Z = X + alpha*Y
%
%   See also uminusFactorized, addTT

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    res.Y_tilde = X.Y_tilde + alpha*Y.Y_tilde;

    res.U1_tilde = X.U1_tilde + alpha*Y.U1_tilde;
    res.U2_tilde = X.U2_tilde + alpha*Y.U2_tilde;
    res.U3_tilde = X.U3_tilde + alpha*Y.U3_tilde;

end
