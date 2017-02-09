% Sample Code to test the Tensor Completion program :

clear;

n1 = 50; 
n2 = 50; 
n3 = 50; 
r = 5;
tol   = [];
nitr  = []; 
ninit = [];

% generate random low-rank orthogonal tensor
U01 = randn(n1,r);U02 = randn(n2,r);U03 = randn(n3,r);
[U1 temp] = qr(U01);[U2 temp] = qr(U02);[U3 temp] = qr(U03);
U1=U1(:,1:r);U2=U2(:,1:r);U3=U3(:,1:r);
T = zeros(n1,n2,n3);
for i=1:n3
    T(:,:,i) = U1*diag(U3(i,:))*U2';
end

% sample entries 
p = 2*(r^0.5*log(n1*n2*n3))/sqrt(n1*n2*n3); 
E = ceil(rand(n1,n2,n3)-1+p);
TE = T.*E;

% exact tensor completion
fprintf(1,'exact tensor completion using alternating minimization\n');
[V1 V2 V3 S dist] = TenALS(TE, E, r, ninit, nitr, tol);
rmse=0;
for i3=1:n3 
   A1 = U1;
   A2 = U2*diag(U3(i3,:));
   B1 = V1;
   B2 = V2*diag(V3(i3,:).*S(:)');
   rmse = rmse + trace(A1'*A1*A2'*A2) + trace(B1'*B1*B2'*B2) -2*trace(B1'*A1*A2'*B2) ;
end
fprintf(1,'root mean squared error = %e\n',sqrt(rmse/r));

% noisy tensor completion
fprintf(1,'\nnoisy tensor completion using alternating minimization\n');
[V1 V2 V3 S dist] = TenALS(TE+(0.0001/sqrt(n1*n2*n3)*randn(n1,n2,n3).*E), E, r, ninit, nitr, tol);
rmse=0;
for i3=1:n3 
   A1 = U1;
   A2 = U2*diag(U3(i3,:));
   B1 = V1;
   B2 = V2*diag(V3(i3,:).*S(:)');
   rmse = rmse + trace(A1'*A1*A2'*A2) + trace(B1'*B1*B2'*B2) -2*trace(B1'*A1*A2'*B2) ;
end
fprintf(1,'root mean squared error = %e\n',sqrt(rmse/r));

