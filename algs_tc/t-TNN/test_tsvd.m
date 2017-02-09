Y = rand(10,5,7);

% [U,S,V] = ntsvd(Y,1,0);

[n1,n2,n3] = size(Y);
Yhat = fft(Y,[],3);
if isinteger(n3/2)
    endValue = int16(n3/2+1);
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
        end
    end
    [uhat,shat,vhat] = svd(full(Yhat(:,:,endValue+1)),'econ');
    Yhat(:,:,endValue+1) = uhat*shat*vhat';
else
   endValue = int16(n3/2+1);
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
        end
    end 
end

Yslide = ifft(Yhat,[],3);

diff = Y - Yslide;

