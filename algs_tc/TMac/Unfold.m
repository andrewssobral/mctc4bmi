function W = Unfold(W, dim, i)
    W = reshape(shiftdim(W,i-1), dim(i), []);
end