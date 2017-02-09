%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: FRSPC
% Time: Sep. 7th, 2015.
% 
% Smooth modeling for incomplete tensor by the PARAFAC decomposition.
% This algorithm is proposed in
%  "Yokota, Tatsuya, et al. "Smooth PARAFAC Decomposition for Tensor Completion." arXiv:1505.06611 (2015)."
%
% minimize || T.*Q - Z.*Q ||_F^2 + penalty_func(U),
%
%    s.t.  Z = [G; U{1}, U{2}, ..., U{N}]
%
% Inputs
% - T       : N-way incomplete tensor
% - Q       : binary tensor which represents elements are available or not (available:1, missing:0)
% - R       : Number of components
% - TV_QV   : it take 'tv' or 'qv' for selecting types of smoothing
% - rho     : N-dimensional vector which represents smoothness of individual modes
% - K       : accelerate parameter for fast optimization (small for fast <--> large for slow, typically 10)
% - SNR     : error threshold via signal-to-noise ratio
% - maxiter : maximum number of iteration
% - tol     : tolerance parameter for convergence evaluation
% - out     : 1 for image completion, 0 for the others
%
% Outputs
% - X       : Results of tensor completion X, where X(Q) = T and X(~Q) = Z
% - Z       : Results of smooth PARAFAC decomposition Z
% - G       : Results of core values G
% - U       : Results of factor matrices U
% - histo   : optimization behavior of error || T.*Q - Z.*Q ||_F^2
% - histo_R : optimization behavior of number of components R
%
% This code was implemented by Tatsuya Yokota
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X Z G U histo histo_R] = FRSPC(T,Q,R,TV_QV,rho,K,SNR,maxiter,tol,out)

  warning off; 

  histo=[];
  histo_R=[];
  N  = ndims(T);
  II = size(T);
  NN = sum(Q(:));
  epsiron = 10^(-SNR/10) * sum_all(T(Q).^2);

  %% initialization
  E  = zeros(II);
  X  = zeros(II);
  Z  = zeros(II);
  ONE= zeros(II);

  X(Q) = T(Q);
  X(~Q) = sum(X(:))/NN;
  for n = 1:N
    Pp  = eye(II(n)-1,II(n));
    Pm  = eye(II(n)); Pm(1,:) = [];
    P{n} = Pp - Pm;
    PtP{n} = P{n}'*P{n};
    PtPI{n} = rho(n)*PtP{n} + eye(II(n));
    Pr{n} = inv(rho(n)*PtP{n} + eye(II(n)));
    U{n} = Pr{n}*randn(II(n),R);
    for r = 1:R
      U{n}(:,r) = U{n}(:,r)/norm(U{n}(:,r));
    end
  end

  G = 1e-16*ones(1,R);
  for r = 1:R
    Z = Z + G(r)*outerprod(U,r);
  end
  obj = sum_all((T(Q) - Z(Q)).^2);
  X(~Q) = Z(~Q);
  E = X - Z;

  %% output some figures (only for color image)
  h1 = figure();clf;
  subplot(1,2,1);cla;hold on;
  subplot(1,2,2);
  drawnow;
  if out == 1
    h2 = figure();clf;
    imagesc(uint8(X));drawnow;
    imwrite(uint8(X),['saved/' TV_QV '_iter_0.png']);
  end

  %% start main algorithm
  for iter = 1:maxiter

    [val ID] = sort(abs(G));
    for k = ID(1:min(end,K))
  
      ONE = G(k)*outerprod(U,k);
      Z = Z - ONE;
      E = E + ONE;

      div = 1;
      for n = 1:N

        if strcmp(TV_QV,'tv') 

          u = innerprod_one_exc(E,U,k,n);
          u = u(:);
          %initilization
          a = u/norm(u);
          % main iteration for constrained version
          object = 0.5*G(k)^2*rho(n)*sum(abs(P{n}*a)) - G(k)*a'*u;
          for nn = 1:1000
            df = P{n}'*sign(P{n}*a);
            dL = (0.5*G(k)^2*rho(n)*df - G(k)*u + G(k)^2*a)/G(k)^2;
            al = [0 0.1 0.01 0.001 0.0001 0.00001];
            for ai = 1:length(al)
              a2 = a - al(ai)*dL;
              a2 = a2/norm(a2);
              score(ai) = 0.5*G(k)^2*rho(n)*sum(abs(P{n}*a2)) - G(k)*a2'*u;
            end
            [object2 ai] = min(score);
            a2 = a - al(ai)*dL;
            a2 = a2/norm(a2);
            if abs(object2 - object)/II(n) < 1e-3
              break;
            else
              a = a2;
              object = object2;
            end
          end
          lam = norm(a);
          u = a/lam;
          U{n}(:,k) = u;
          v{n} = u;
          div  = div + rho(n)*sum(abs(P{n}*u));

        elseif strcmp(TV_QV,'qv')

          u = innerprod_one_exc(E,U,k,n);
          u = u(:);
          % initialization
          a = Pr{n}*u;
          a = a/norm(a);
          % main iteration for constrained version
          object = 0.5*G(k)^2*rho(n)*a'*PtP{n}*a - G(k)*a'*u;
          mu = 0.1;
          for nn = 1:1000
            dL = (G(k)^2*rho(n)*PtPI{n}*a - G(k)*u'*a)/G(k)^2;
            al = [0 0.1 0.01 0.001 0.0001 0.00001];
            for ai = 1:length(al)
              a2 = a - al(ai)*dL;
              a2 = a2/norm(a2);
              score(ai) = 0.5*G(k)^2*rho(n)*a2'*PtP{n}*a2 - G(k)*a2'*u;
            end
            [object2 ai] = min(score);
            a2 = a - al(ai)*dL;
            if abs(object2 - object)/II(n) < 1e-3
              break;
            else
              a = a2;
              object = object2;
            end
          end

          u = a2;
          lam = norm(u);
          u = u/lam;
          U{n}(:,k) = u;
          v{n} = u;
          div  = div + rho(n)*u'*PtP{n}*u;

        else
          error('2rd input is ''tv'' or ''qv'' ');
        end

      end

      G(k) = tensor_allprod(E,v,1)/div;
      if G(k) < 0
        G(k) = -G(k);
        U{1}(:,k) = - U{1}(:,k);
      end

      ONE = G(k)*outerprod(U,k);
      E = E - ONE;
      E(~Q) = 0;
      Z = Z + ONE;
      X(~Q) = Z(~Q);

    end
   
    %% calculate MSE
    obj2 = sum_all(E.^2);

    %% convergence speed
    speed = abs(obj2 - obj)/abs(epsiron - obj2);

    %% checking convergence
    if obj2 < epsiron || abs(obj2-obj)/NN < tol
      break;
    else
      obj = obj2;
      if mod(iter,5)==0
        fprintf('%d:  %f :: %f :: %f :: Pid %d \n',iter,obj2/NN,epsiron/NN,speed,R);
      end
      histo(iter) = obj;
      histo_R(iter,:) = size(G);
    end

    %% output figures
    set(0,'CurrentFigure',h1);
    subplot(1,2,1);cla;hold on;
    plot((histo));
    plot((epsiron)*ones(1,length(histo)))
    grid on;
    set(gca,'YScale','log')
    title('MSE')
    subplot(1,2,2);
    plot(histo_R(:,2));
    title('number of components R')
    
    drawnow;
    if out == 1
      set(0,'CurrentFigure',h2);
      imagesc(uint8(X));
      title(['Number of R = ' num2str(R)]);
      drawnow;
    end

    if mod(iter,10) == 0 & out == 1
      imwrite(uint8(Z),['saved/' TV_QV '_iter_' num2str(iter) '.png']);
    end

    if mod(iter,100) == 0
      pack;
    end

  end


