function [x, objV] = proxF_l1(x,t) 
%% Usage proxF_l1(x,t)
% ProxF_l1(x,t) applies L1 shrinkage to n-dimensional matrix x for
% a shrinkage factor of t.  Returns the shrunk object.
%
%
% INPUTS:
%
% x: n-dimensional matrix
% t: shrinkage threshold
%
% OUTPUTS:
%
% x: n-dimensional matrix after shrinkage.
%
% Authors: G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015

tq = t * 1;
s  = 1 - min( tq./abs(x), 1 );
x  = x .* s;
objV = sum(x(:));
end
