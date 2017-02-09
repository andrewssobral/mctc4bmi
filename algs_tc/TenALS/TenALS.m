function [V1 V2 V3 S dist] = TenAls(TE, E, r, ninit, nitr, tol)
% An algorithm for Low-rank Tensor Reconstruction from a partially revealed set. 
% See "Provable Tensor Factorization with Missing Data"(http://arxiv.org/abs/1406.2784) for details
% Usage :
% [V1 V2 V3 S dist] = TenAls(TE, E, r, ninit, nitr, tol)
% [V1 V2 V3 S dist] = TenAls(TE, E, r)
% [V1 V2 V3 S dist] = TenAls(TE, E)
% 
% INPUT :
% TE    :  The partially revealed 3rd-order Tensor.
%          Sparse Tensor with zeroes at the unrevealed indices.
% E     :  E is the sparse 0-1 matrix indicating which eantries are revealed.
%          E(i,j,k) = 1 if T(i,j,k) is revealed.
% r     :  The target rank to be used for reconstruction. Use [] to guess the rank.
% ninit :  The number of initial vectors to use in the initialization step.
%          Larger value gives better estimate at the expense of longer
%          processing time. Use [] to use defalut (50).
% niter :  The max. no. of iterations. Use [] to use default (50).
% tol   :  Stop iterations if || (estimated T).*E - TE ||_F/||TE||_F < tol, where
%          Use [] to use the default (1e-8)
%
% OUTPUT :
% V1     : A size(TE,1)*r matrix
% V2     : A size(TE,2)*r matrix
% V3     : A size(TE,3)*r matrix
% S      : An r-dim vector
% such that T_hat(i,j,k) = \sum_{a=1}^r S(a)*V1(i,a)*V2(i,a)*V3(i,a) 
% dist   : A vector containing  || T_hat.*E - TE ||_F/||TE||_F at each
%          successive iteration
%
% Date : 16st June, 2014
% COPYRIGHT 2014 Prateek Jain, Sewoong Oh

    if (isempty(ninit)==1) ninit=5; end; % ninit=10;
    if (isempty(nitr)==1) nitr=5;  end; % nitr=50;
    if (isempty(tol)==1)  tol = 1e-4; end; % tol = 1e-8;

    [n1 n2 n3] = size(TE);
    p = nnz(E)/(n1*n2*n3);
    normTE=0;
    for i3=1:n3
        normTE = normTE + norm(TE(:,:,i3),'fro')^2;
    end

% initialization by Robust Tensor Power Method (modified for non-symmetric tensors)
    U01 = zeros(n1,r);
    U02 = zeros(n2,r);
    U03 = zeros(n3,r);
    S0 = zeros(r,1);
    for i=1:r 
        tU1 = zeros(n1,ninit);
        tU2 = zeros(n2,ninit);
        tU3 = zeros(n3,ninit);
        tS  = zeros(ninit,1);
        for init=1:ninit
            [tU1(:,init) tU2(:,init) tU3(:,init)] = RTPM(TE-CPcomp(S0,U01,U02,U03), nitr);  
            tU1(:,init) = tU1(:,init)./norm(tU1(:,init));
            tU2(:,init) = tU2(:,init)./norm(tU2(:,init));
            tU3(:,init) = tU3(:,init)./norm(tU3(:,init));
            tS(init) = TenProj(TE-CPcomp(S0,U01,U02,U03),tU1(:,init),tU2(:,init),tU3(:,init) );
        end
        [C I] = max(tS);   
        U01(:,i) = tU1(:,I)/norm(tU1(:,I));
        U02(:,i) = tU2(:,I)/norm(tU2(:,I));
        U03(:,i) = tU3(:,I)/norm(tU3(:,I));
        S0(i) = TenProj(TE-CPcomp(S0,U01,U02,U03),U01(:,i),U02(:,i),U03(:,i));
    end
% apply alternating least squares (or equivalently alternating minimization)
    V1 = U01;V2 = U02;V3 = U03; 
    S = S0; 
    for itr1=1:nitr
        V1_=V1;V2_=V2;V3_=V3;
        for q=1:r
            S_ = S;
            S_(q)=0;
            A = CPcomp(S_,V1,V2,V3).*E;
            v1=V1(:,q);v2=V2(:,q);v3=V3(:,q);
            V1(:,q)=0;V2(:,q)=0;V3(:,q)=0;
            den1=zeros(n1,1);den2=zeros(n2,1);
            s = S(q); 
            for i3=1:n3
                V1(:,q) = V1(:,q) + v3(i3)*(TE(:,:,i3)-A(:,:,i3))*v2;
                den1    = den1    + v3(i3)^2*E(:,:,i3)*(v2.*v2);
            end
            v1 = V1(:,q)./den1;
            v1=v1/norm(v1);
            for i3=1:n3
                V2(:,q) = V2(:,q) + v3(i3)*(TE(:,:,i3)-A(:,:,i3))'*v1;
                den2    = den2    + v3(i3)^2*E(:,:,i3)'*(v1.*v1);
            end
            v2 = V2(:,q)./den2;
            v2=v2/norm(v2);
            for i3=1:n3
                V3(i3,q) = ( v1'*(TE(:,:,i3)-A(:,:,i3))*v2 ) /( (v1.*v1)'*E(:,:,i3)*(v2.*v2) );
                if (isnan(V3(i3,q))) 
                    fprintf(1,'NaN in als_tensor, denomenator=0\n'); 
                    break; 
                end 
            end
            if (nnz(den1)~=n1 || nnz(den2)~=n2)
                fprintf(1,'NaN in als_tensor, denomenator=0.\n'); 
                break;                     
            end
            if (norm(V1(:,q))==0 || norm(V2(:,q))==0 || norm(V3(:,q))==0) fprintf(1,'ERROR: estimate zero!'); end
            V1(:,q) = v1;
            V2(:,q) = v2;
            S(q)    = norm(V3(:,q)); 
            V3(:,q) = V3(:,q)/norm(V3(:,q)); 
        end
        ERR = (TE - E.*CPcomp(S,V1,V2,V3));
        normERR=0;
        for i3=1:n3
            normERR = normERR + norm(ERR(:,:,i3),'fro')^2;
        end
        fprintf(1,'iter %i fiterror = %e \n',itr1,sqrt(normERR/normTE));
        if (sqrt(normERR/normTE)<tol) break; end
    end
    dist = sqrt(normERR/normTE);
end

function T = CPcomp(S,U1,U2,U3)
   n1 = size(U1,1);n2 = size(U2,1);n3 = size(U3,1);
   r = min([length(S) size(U1,2) size(U2,2) size(U3,2)]);
   T = zeros(n1,n2,n3);
   for i=1:n3
        T(:,:,i) = U1*diag(U3(i,:).*S(:)')*U2';
   end
end

function M = TenProj(T,U1,U2,U3)
   n1=size(U1,1);n2=size(U2,1);n3=size(U3,1);
   r1=size(U1,2);r2=size(U2,2);r3=size(U3,2);
   M =zeros(r1,r2,r3);
   for i=1:r3
       A = zeros(n1,n2);
       for j=1:n3
           A = A+T(:,:,j)*U3(j,i);
       end
       M(:,:,i) = U1'*A*U2;
   end
end

function [u1 u2 u3] = RTPM(T, mxitr)
    n1=size(T, 1);n2=size(T, 2);n3=size(T, 3);
    uinit=randn(n1,1);u1=uinit/norm(uinit);
    uinit=randn(n2,1);u2=uinit/norm(uinit);
    uinit=randn(n3,1);u3=uinit/norm(uinit);
    for itr=1:mxitr
        v1=zeros(n1,1);v2=zeros(n2,1);v3=zeros(n3,1);
        for i3=1:n3
            v3(i3)=u1'*T(:, :, i3)*u2; 
            v1 = v1 + u3(i3)*T(:, :, i3)*u2;
            v2 = v2 + u3(i3)*T(:, :, i3)'*u1;
        end
        u10 = u1;
        u1 = v1/norm(v1);
        u20 = u2;
        u2 = v2/norm(v2);
        u30 = u3;
        u3 = v3/norm(v3);
        if(norm(u10-u1)+norm(u20-u2)+norm(u30-u3)<1e-7) break; end
    end
end


