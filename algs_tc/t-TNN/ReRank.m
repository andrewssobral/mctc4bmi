function Y = ReRank(X, r)

Y3 = X2Yi(X,3);

Y3hat = fft(Y3,[],3);

[n1,n2,n3] = size(Y3hat);

nn = min(n1,n2);
ss = zeros(nn*n3,1);

for i = 1 : n3
    [uhat,shat,vhat] = svd(full(Y3hat(:,:,i)), 'econ');
    ss((i-1)*nn+1 : i*nn) = diag(shat);
end

ss = sort(abs(ss), 'descend');

tau = ss(round(r*nn*n3));


for i = 1 : n3
    [uhat,shat,vhat] = svd(full(Y3hat(:,:,i)), 'econ');
    shat(abs(shat)<tau) = 0;
    Y3hat(:,:,i) = uhat*shat*vhat';  
end
Y3 = ifft(Y3hat,[],3);

Y = Yi2X(Y3,3);
end
