%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Low Rank Tensor Completion no relaxation version (SiLRTCnr)
% Time: 03/11/2012
% Reference: "Tensor Completion for Estimating Missing Values 
% in Visual Data", PAMI, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, errList] = SiLRTCnr(T, Omega, alpha, factor, maxIter, epsilon, X)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_X:  \alpha1||X_(1)||_* + \alpha2||X_(2)||_* + \alpha3||X_(3)||_* + ....
%         s.t.  X_\Omega = T_\Omega
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    X = T;
    X(logical(1-Omega)) = mean(T(Omega));
end

dim = size(T);
N = ndims(T);
for i = 1:N
    % tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.4;
    tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) / max(alpha) * 0.3;
end

errList = zeros(maxIter, 1);
%L = errList;
dim = size(T);
M = cell(ndims(T), 1);
tau = tmp;
normT = norm(T(:));

%normT = norm(T(:));
for k = 1:maxIter
    if mod(k, 20) == 0
        fprintf('SiLRTCnr: iterations = %d   difference=%f\n', k, errList(k-1));
    end
    
    tau = tau*factor;
    % disp(k);
    Xsum = 0;
    for i = 1:ndims(T)
        M{i} = Fold(Pro2TraceNorm(Unfold(X, dim, i), tau(i)), dim, i);
        Xsum = Xsum +  M{i};
    end
    Xlast = X;
    X = Xsum / N;
    X(Omega) = T(Omega);
    errList(k) = norm(X(:)-Xlast(:)) / normT;
    if (errList(k) < epsilon)
        errList = errList(1:k);
        break;
    end
    %L(k) = norm(X(:)-T(:)) / normT;
end
fprintf('SiLRTCnr ends: total iterations = %d   difference=%f\n\n', k, errList(k));