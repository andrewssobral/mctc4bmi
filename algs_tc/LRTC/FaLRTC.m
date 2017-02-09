%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Low Rank Tensor Completion (FaLRTC)
% Time: 03/11/2012
% Reference: "Tensor Completion for Estimating Missing Values 
% in Visual Data", PAMI, 2012.
% min_{X} : \Psi(X) 
% s.t.         : X_\Omega = M_\Omega
% \Psi(X) = max_{Z_{i(i)} <= 1}: <X, \sum_i Y_i> - 0.5 \mu_i \|Y_i\|_F^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% code 
function [Y, errList] = FaLRTC(M, Omega, alpha, mu, L, C, maxIter, epsilon, X)

% initialization
if(nargin < 9 || isempty(X)) disp('FaLRTC.initialization');
%if(nargin < 9)
    X = M;
    X(logical(1-Omega)) = mean(M(Omega));
end

Y = X;
Z = X;
B = 0;

N = ndims(M);
dim = size(M);
Gx = zeros(dim);

%funValueList = zeros(K+1,1);
%funValueList(1) = inf;

%%%%%%
normM = norm(M(:));
%errList(1) = norm(Y(:)-M(:)) / normM;
errList = zeros(maxIter, 1);
%%%%%%

% a2m = alpha.^2 ./ mu;
% ma = mu ./ alpha;
% 
% fylast = inf;

Lmax = 10*sum(1./mu);

tmp = zeros(1, N);
for i = 1:N
    % tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.4;
    tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.3;
end
P = 1.15;
flatNum = 15;
slope = (tmp - mu) / (1-(maxIter-flatNum)^(-P));
offset = (mu*(maxIter-flatNum)^P - tmp) / ((maxIter-flatNum)^P-1); 
% offset = (mu*K^P - tmp) / ((K-flatNum)^P-1); 

mu0 = mu;
for k = 1:maxIter
    if mod(k, 20) == 0
         fprintf('FaLRTC: iterations = %d   difference=%f\n', k, errList(k-1));
    end
    
    % update mu
    mu = max(slope / (k^P) +offset, mu0);
    % mu = mu0;
    % mu = mu0 * (K-k+1)^2
    
    a2m = alpha.^2 ./ mu;
    ma = mu ./ alpha;
    
    Ylast = Y;
    %%  test L
    %L = L*C;
    while true
        b = (1+sqrt(1+4*L*B)) / (2*L);
        X = b/(B+b) * Z + B/(B+b) * Ylast;
        % compute f'(x) namely "Gx" and f(x) namely "fx"
        Gx = Gx * 0;
        fx = 0;
        for i = 1 : N
            [temp, sigma2] = Truncate(Unfold(X, dim, i), ma(i));
            temp = Fold(temp, dim, i);
            Gx = Gx + a2m(i) * temp;
            fx = fx + a2m(i)*(sum(sigma2) - sum(max(sqrt(sigma2)-ma(i), 0).^2));
        end
        Gx(Omega) = 0;

        % compute f(Ytest) namely fy
        Y = X - Gx / L;
        fy = 0;
        for i = 1 : N
            [sigma] = SingularValue(Unfold(Y, dim, i));
            fy = fy + a2m(i)*(sum(sigma.^2) - sum(max(sigma-ma(i), 0).^2));
        end
        % test if L(fx-fy) > \|Gx\|^2
        if (fx - fy)*L < sum(Gx(:).^2)
            if L > Lmax
                k
                % funValueList = funValueList(2:k);
                % disp('Exceed the Maximum Lipschitiz Constant');
                fprintf('FaLRTC: iterations = %d   difference=%f\n Exceed the Maximum Lipschitiz Constan\n\n', k, errList(k-1));
                errList = errList(1:k);
                return;
                % break;
            end
            L = L/C;
        else
             break;
        end
    end
    
    errList(k) =  norm(Y(:)-Ylast(:)) / normM;
    if errList(k) < epsilon
        break;
    end
    
    %% test Y = ? and if return
    % test whether the fylast satisfy the following condition
    % funValueList(k+1) = fy / 2;
    % Y = Ytest;
    % fylast = fy;
    % if abs(funValueList(k+1) - funValueList(k)) < epsilon
%         disp('Exceed the Minimum Tolerance');
%         break;
%     end
    
%     if fylast < fy
%         funValueList(k+1) = fylast / 2;
%     else
%         funValueList(k+1) = fy / 2;
%         Y = Ytest;
%         fylast = fy;
%         if abs(funValueList(k+1) - funValueList(k)) < epsilon
%             disp('Exceed the Minimum Tolerance');
%             break;
%         end
%     end
    
    %% update Z, Y, and B
    Z = Z - b*Gx;
    B = B+b;
    
    %%%%%%%%%%%%
    % errList(end+1) = norm(Y(:)-M(:)) / normM;
    %%%%%%%%%%%%
end
% k
% funValueList = funValueList(2:k+1);

%%%%%%%%
%L = errList;
%%%%%%%%
errList = errList(1:k);
fprintf('FaLRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));
