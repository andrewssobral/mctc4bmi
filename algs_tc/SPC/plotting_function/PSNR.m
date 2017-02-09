function psnr = PSNR(X,G);

X = double(uint8(X));
G = double(uint8(G));

psnr = 10*log10( (255^2)/( mean((X(:)-G(:)).^2) ) );

end
