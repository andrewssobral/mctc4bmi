function [X Z G U histo histo_R] = SPC(T,Q,TV_SV,rho,K,SNR,rate,maxiter,tol,out)

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

  R  = 1;

  X(Q) = T(Q);
  X(~Q) = sum(X(:))/NN;
  for n = 1:N
    U{n} = randn(II(n),R);
    U{n} = U{n}/norm(U{n});
    Pp  = eye(II(n)-1,II(n));
    Pm  = eye(II(n)); Pm(1,:) = [];
    P{n} = Pp - Pm;
    PtP{n} = P{n}'*P{n};
    Pr{n} = inv(rho(n)*PtP{n} + eye(II(n)));
  end
  G = tensor_allprod(X,U,1);

  for r = 1:R
    Z = Z + G(r)*outerprod(U,r);
  end
  obj = sum_all((T(Q) - Z(Q)).^2);
  X(~Q) = Z(~Q);
  E = X - Z;

  %% output some figures

  h1 = figure();clf;
  subplot(1,2,1);cla;hold on;
  subplot(1,2,2);
  drawnow;
  if out == 1
    h2 = figure();clf;
    imagesc(uint8(X));drawnow;
    imwrite(uint8(X),['saved/' TV_SV '_iter_0.png']);
  end

  %% start main algorithm
  for iter = 1:maxiter

    [val ID] = sort(G);
    for k = ID(1:min(end,K))
  
      ONE = G(k)*outerprod(U,k);
      Z = Z - ONE;
      E = E + ONE;

      div = 1;
      for n = 1:N

        if strcmp(TV_SV,'tv') 

          u = innerprod_one_exc(E,U,k,n);
          u = u(:);
          a = u(:);
          object = 0.5*G(k)*sum(abs(P{n}*a)) + a'*u + 0.5*G(k);
          for nn = 1:1000
            df = P{n}'*sign(P{n}*a);
            dL = (0.5*G(k)*rho(n)*df - u + G(k)*a)/G(k);
            al = [0 0.1 0.01 0.001 0.0001 0.00001];
            for ai = 1:length(al)
              a2 = a - al(ai)*dL;
              score(ai) = 0.5*G(k)*sum(abs(P{n}*a2)) + a2'*u + 0.5*G(k);
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
          lam = norm(a);
          u = a/lam;
          U{n}(:,k) = u;
          v{n} = u;
          div  = div + rho(n)*sum(abs(P{n}*u));

        elseif strcmp(TV_SV,'sv')

          u = innerprod_one_exc(E,U,k,n);
          u = Pr{n} * u(:);
          lam = norm(u);
          u = u/lam;
          U{n}(:,k) = u;
          v{n} = u;
          div  = div + rho(n)*u'*PtP{n}*u;

        else
          error('2rd input is ''tv'' or ''qc'' ');
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

    %% step for R <= R + 1
    if speed < rate || abs(obj2-obj)/NN < tol

      R = R + 1;
      for n = 1:N
        u = randn(II(n),1);
        u = u/norm(u);
        U{n}(:,R) = u;
        v{n} = u;
      end
      G(R) = tensor_allprod(E,v,1);
      E = E - G(R)*outerprod(U,R);
      E(~Q) = 0;
      
    end

    %% checking convergence
    if obj2 < epsiron
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
      imwrite(uint8(Z),['saved/' TV_SV '_iter_' num2str(iter) '.png']);
    end

    if mod(iter,100) == 0
      pack;
    end

  end


