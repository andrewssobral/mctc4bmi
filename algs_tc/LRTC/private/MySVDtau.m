function [S, V, D, Sigma2] = MySVDtau(A, tau)
[m, n] = size(A);
if 2*m < n
    AAT = A*A';
    [S, Sigma2, D] = svd(AAT);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    tol = max(size(A)) * eps(max(V));
    R = sum(V > max(tol, tau));

    %tol = min(size(A)) * eps(max(V));
    % R = sum(V > 1/tau)
        
%     tol = min(size(A)) * eps(max(V));
%     R = sum(V > tau)
    
    V = V(1:R);
    S = S(:,1:R);
    D = A'*S*diag(1./V);
    V = diag(V);
    return;
end
if m > 2*n
    [S, V, D, Sigma2] = MySVDtau(A', tau);
    mid = D;
    D = S;
    S = mid;
    return;
end
[S,V,D] = svd(A);
Sigma2 = diag(V).^2;
R = sum(diag(V) > tau);
S = S(:, 1:R);
V = V(1:R, 1:R);
D = D(:, 1:R);