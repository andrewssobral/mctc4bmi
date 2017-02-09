function TX = Trans_Faces(X)
[l, m, n] = size(X);
TX = zeros(m, l, n);
for i = 1 : n
    TX(:,:,i)  =  X(:,:,i)';
end
end
