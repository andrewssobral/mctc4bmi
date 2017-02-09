function Z = outerprod(U,r);

clear kron

N = length(U);
for n = 1:N
  II(n) = size(U{n},1);
end
NN = prod(II);
Z  = zeros(II);

if N == 3

  for i3 = 1:II(3)

    Z(:,:,i3) = (U{3}(i3,r)*U{1}(:,r))*U{2}(:,r)';

  end

else

%%%%%%

uu = U{1}(:,r);
for n = 2:N

  b = U{n}(:,r);
  uu = kron(b,uu);

end
Z = reshape(uu,II);

%%%%%%

end
