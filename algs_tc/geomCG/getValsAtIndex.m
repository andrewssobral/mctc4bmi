function vals = getValsAtIndex( X, subs )
%GETVALSATINDEX Evaluates a Tucker tensor at given subscript indices
%   VALS = GETVALSATINDEX( X, SUBS ) evaluates the Tucker tensor (ttensor) X
%   at the vector of m subscript indices subs, where subs is a (m x 3) vector.
%

%   GeomCG Tensor Completion. Copyright 2013 by
%   Michael Steinlechner
%   Questions and contact: michael.steinlechner@epfl.ch
%   BSD 2-clause license, see LICENSE.txt

    vals = getValsAtIndex_mex( subs', X.core.data, X.U{1}', X.U{2}', X.U{3}');

end
