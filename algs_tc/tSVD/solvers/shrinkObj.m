function [y,objV] = shrinkObj(a,rho,norm,sX,varargin)
% Shrinkage funcion
% Authors: G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015

if ~exist('varargin','var')
   opts = [];
else
    opts = varargin;
end

a=reshape(a,sX);
norm = lower(norm);

switch norm
    
    case 'tsvd_1'
         [y,objV] =proxF_tSVD_1(a,rho,opts);

    otherwise
        error('Invalid Norm')

end

y = y(:);

end
 