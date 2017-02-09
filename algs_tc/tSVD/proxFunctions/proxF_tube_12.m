function [y, objV] = proxF_tube_12(x,t) 
% Proximal function of tubal norm
% Authors: G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015

norm_matrix = sum( x.*x , 3) ;

s = pos( 1 - t./norm_matrix );
v = repmat(s,[1,1,size(x,3)]);
y = v .*  x;

objV = sum(sum( sum(y.*y,3)  ));


end

