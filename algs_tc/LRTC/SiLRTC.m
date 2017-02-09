%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Simple Low Rank Tensor Completion (SiLRTC) 
% Time: 03/11/2012
% Reference: "Tensor Completion for Estimating Missing Values 
% in Visual Data", PAMI, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, errList] = SiLRTC(T, Omega, alpha, gamma, maxIter, epsilon, X)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% min(X, M1, M2, M3,... Mn): (\gamma1||X_(1)-M1||^2 + \gamma2||X_(2)-T2||^2 + \gamma3||X_(3)-T3||^2 + ...)/2 + 
%               \alpha1||M1||_* + \alpha2||M2||_* + \alpha3||M3||_* + ....
%         s.t.  X_\Omega = T_\Omega
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    X = T;
    X(logical(1-Omega)) = mean(T(Omega));
end

errList = zeros(maxIter, 1);
normT = norm(T(:));
%L = errList;
dim = size(T);
M = cell(ndims(T), 1);
gammasum = sum(gamma);
tau = alpha./ gamma;

%normT = norm(T(:));
for k = 1:maxIter
    if mod(k, 20) == 0
        fprintf('SiLRTC: iterations = %d   difference=%f\n', k, errList(k-1));
    end
    Xsum = 0;
    for i = 1:ndims(T)
        M{i} = Fold(Pro2TraceNorm(Unfold(X, dim, i), tau(i)), dim, i);
        Xsum = Xsum + gamma(i) * M{i};
    end
    Xlast = X;
    X = Xsum / gammasum;
    X(Omega) = T(Omega);
    errList(k) = norm(X(:)-Xlast(:)) / normT;
    if (errList(k) < epsilon)
        errList = errList(1:k);
        break;
    end
    %L(k) = norm(X(:)-T(:)) / normT;
end
fprintf('SiLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));