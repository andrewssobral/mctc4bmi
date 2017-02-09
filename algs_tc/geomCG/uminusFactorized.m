function res = uminusFactorized( xi )
%UMINUSFACTORIZED Unary minus for a tangent tensor.
%   RES = UMINUSFACTORIZED( XI ) performs unary minus operation
%   on the tangent tensor XI given in factorized form
%
%   See also addFactorized
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

        res.Y_tilde =  -xi.Y_tilde;
        res.U1_tilde = -xi.U1_tilde;
        res.U2_tilde = -xi.U2_tilde;
        res.U3_tilde = -xi.U3_tilde;
end
