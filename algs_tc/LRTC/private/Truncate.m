function [X, Sigma2] = Truncate(Z, tau)
% [Y, n, Sigma2] = Pro2TraceNorm(Z, tau);
% X = Z-Y;

%% new
[m, n] = size(Z);
if 2*m < n
    AAT = Z*Z';
    [S, Sigma2, D] = svd(AAT);
    Sigma2 = diag(Sigma2);
    V = sqrt(Sigma2);
    tol = max(size(Z)) * eps(max(V));
    n = sum(V > max(tol, tau));
    mid = max(V(1:n)-tau, 0) ./ V(1:n) ;
    X = (eye(m)-S(:, 1:n) * diag(mid) * S(:, 1:n)') * Z;
    return;
end
if m > 2*n
    [X, Sigma2] = Truncate(Z', tau);
    X = X';
    return;
end
[S,V,D] = svd(Z, 0);
Sigma2 = diag(V).^2;
n = sum(diag(V) > tau);
X = Z - S(:, 1:n) * diag(diag(V(1:n,1:n))-tau) * D(:, 1:n)';