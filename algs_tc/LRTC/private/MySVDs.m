function [S, V, D] = MySVDs(A, R)
[m, n] = size(A);
if 10*m < n
    AAT = A*A';
    [S, V, D] = svd(AAT);
    V = diag(V) .^ 0.5;
    tol = max(size(A)) * eps(max(V));
    R = min([R, sum(V > tol)]);
    V = V(1:R);
    S = S(:,1:R);
    D = A'*S*diag(1./V);
    V = diag(V);
    return;
end
if m > 10*n
    [S, V, D] = MySVD(A');
    mid = D;
    D = S;
    S = mid;
    return;
end
[S,V,D] = svd(A);
S = S(:, 1:R);
V = V(1:R, 1:R);
D = D(:, 1:R);