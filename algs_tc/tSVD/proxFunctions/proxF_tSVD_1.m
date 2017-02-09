function [X, objV] = proxF_tSVD_1(Y,rho,opts)
% proximal function of tensor nuclear norm
% Authors: G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015

%% Input checking
parOP = false;
if ~exist('opts','var') && ~isempty(opts)
    nARG = length(opts);
    if nARG > 0
        parOP=opts{1};
    end    
end

%%
sa = size(Y);
la = length(sa);
[n1,n2,n3] = size(Y);

% this returned this STILL FFT'd version
[U,S,V]=ntsvd(Y,1,parOP);

sTrueV =zeros(min(n1,n2),n3);
for i = 1:n3
    s = S(:,:,i);  % tube
    s = diag(s);
    sTrueV(:,i) = s;
end
%%
[sTrueV] = proxF_l1(sTrueV,rho);
%objV = sum(sTrueV(:));
objV = sum(sTrueV(:));
%%
for i = 1:min(n1,n2) 
    for j = 1:n3
        S(i,i,j) = sTrueV(i,j);
    end
end
    
for i = la:-1:3
    U = ifft(U,[],i);
    S = ifft(S,[],i);
    V = ifft(V,[],i);
end

% X = tprod( tprod(U,S), ttrans_HO(V));

X = tprod( tprod(U,S), tran(V));
