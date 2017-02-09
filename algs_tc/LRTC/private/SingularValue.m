function V = SingularValue(A)
[m, n] = size(A);
if 2*m < n
    AAT = A*A';
    [S, V, D] = svd(AAT);
    V = sqrt(diag(V));
    return;
end
if m > 2*n
    AAT = A'*A;
    [S, V, D] = svd(AAT);
    V = sqrt(diag(V));
    return;
end
[S,V,D] = svd(A);
V = diag(V);