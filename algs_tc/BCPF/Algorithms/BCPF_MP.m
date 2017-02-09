function [model] = BCPF_MP(Y, varargin)
% Bayesian CP Factorization using Gaussian Mixture Priors for Image Completion
% Author : Qibin Zhao  2014
%
% -----------------------------------------------------------------------
%  [model] = BCPF_MP(Y, 'PARAM1', val1, 'PARAM2', val2, ...)
%
%  INPUTS
%     Y              - Input tensor
%     'obs'          - Binary (0-1) tensor indicating missing entries
%                      (0: missing; 1: observed)
%     'init'         - Initialization method
%                     - 'ml'  : SVD initilization (default)
%                     - 'rand': Random matrices
%     'maxRank'      - The initialization of rank (larger than true rank)
%     'dimRed'       - 1: Remove unnecessary components automaticly (default)
%                    - 0: Not remove
%     'maxiters'     - max number of iterations (default: 100)
%     'tol'          - lower band change tolerance for convergence dection
%                      (default: 1e-5)
%     'noise'        - whether noise is updated
%                        - 'on': update noise parameter (default)
%                        - 'off': fixed noise parameter (1e-5)
%     'predVar'      - Predictive distribution
%                         - 1:  compute and output
%                         - 0:  doesnot compute  (default)
%     'verbose'      - visualization of results
%                       - 0: no
%                       - 1: text (default)
%                       - 2: online display image
%                       - 3: show factors by image
%                       - 4: show factors by hinton plot (very slow)
%   OUTPUTS
%      model         - Model parameters and hyperparameters
% -----------------------------------------------------------------------
%
%   Example:
%
%     [model] = BCPF_MP(Y, 'obs', O, 'init', 'rand', 'maxRank', 10, 'dimRed', 1, 'maxiters', 100, ...
%                                'tol', 1e-6, 'verbose', 3);
%
% < Bayesian CP Factorization of Incomplete Image using Gaussian Mixture Priors >
% Copyright (C) 2014  Qibin Zhao
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
warning off; %#ok<WNOFF>
randn('state',1); rand('state',1); %#ok<RAND>
dimY = size(Y);
N = ndims(Y);

%% Set parameters from input or by using defaults
ip = inputParser;
ip.addParamValue('obs', ones(dimY), @(x) (isnumeric(x) || islogical(x)) );
ip.addParamValue('init', 'rand', @(x) (ismember(x,{'ml','rand'})));
ip.addParamValue('maxRank', max(dimY), @isscalar);
ip.addParamValue('maxiters', 100, @isscalar);
ip.addParamValue('tol', 1e-5, @isscalar);
ip.addParamValue('verbose', 1, @isscalar);
ip.addParamValue('noise', 'on', @(x)ismember(x,{'on','off'}));
ip.addParamValue('dimRed', 1, @isscalar);
ip.addParamValue('predVar', 0, @isscalar);
ip.addParamValue('nd', 1, @isscalar);
ip.parse(varargin{:});

O     = ip.Results.obs;
init  = ip.Results.init;
R   = ip.Results.maxRank;
maxiters  = ip.Results.maxiters;
tol   = ip.Results.tol;
verbose  = ip.Results.verbose;
DIMRED   = ip.Results.dimRed;
noise = ip.Results.noise;
predVar = ip.Results.predVar;
nd = ip.Results.nd;

%% Initialization
Y = tensor(Y.*O);
O = tensor(O);
nObs = sum(O(:));

a_gamma0     = 1e-6;
b_gamma0     = 1e-6;
if  strcmp(noise,'on')
    a_beta0      = 1e-6;
    b_beta0      = 1e-6;
else
    a_beta0      = 1e-1;
    b_beta0      = 1e-6;
end
gammas = ones(R,1);
beta = 1e4;
dscale =1;

W = cell(N,1);
W1 = cell(N,1);
oR = nObs/prod(dimY);
for n=1:N
    W{n} = zeros(dimY(n),dimY(n));
    for i=1:dimY(n)
        for j=1:dimY(n)
            W{n}(i,j) = exp(-2*oR*abs(i-j)^2);
        end
    end
    W1{n} = W{n} - diag(diag(W{n}));
    W1{n} = nd*5*bsxfun(@times, W1{n}, 1./sum(W1{n},2));
end

switch init,
    case 'ml'    % Maximum likelihood
        Z = cell(N,1);
        ZSigma = cell(N,1);
        if ~isempty(find(O==0))
            %   Y(find(O==0)) = sum(Y(:))/nObs;
            Y1 = Y;
            for n=1:N
                Y1 = ttm(Y1, W{n}, n);
            end
            Y(find(O==0)) = Y1(find(O==0));
        end
        for n = 1:N
            ZSigma{n} = (repmat(eye(R), [1 1 dimY(n)]));
            [U, S, V] = svd(double(tenmat(Y,n)), 'econ');
            if R <= size(U,2)
                Z{n} = U(:,1:R)*(S(1:R,1:R)).^(0.5);
            else
                Z{n} = [U*(S.^(0.5)) randn(dimY(n), R-size(U,2))];
            end
        end
        Y = Y.*O;
    case 'rand'   % Random initialization
        Z = cell(N,1);
        ZSigma = cell(N,1);
        for n = 1:N
            Z{n} = rand(dimY(n),R);
            ZSigma{n} = (repmat(eye(R), [1 1 dimY(n)]));
            if n<3
                %   Z{n} = cholcov(W{n})'*Z{n};
                Z{n} = Z{n} +  W1{n}*Z{n};
            end
        end
end

% --------- E(aa') = cov(a,a) + E(a)E(a')----------------
EZZT = cell(N,1);
for n=1:N
    EZZT{n} = (reshape(ZSigma{n}, [R*R, dimY(n)]))';
end

Fit =0;
LB = 0;
X = double(ktensor(Z));

%% Create figures
if verbose >2,
    scrsz = get(0,'ScreenSize');
    h1 = figure('Position',[scrsz(3)*0.2 scrsz(4)*0.3 scrsz(3)*0.6 scrsz(4)*0.4]);
    figure(h1);
    switch verbose,
        case 4,
            subplot(2,3,1); hintonDiagram(Z{1}); title('Mode-1'); ylabel('Length of #-mode');
            subplot(2,3,2); hintonDiagram(Z{2}); title('Mode-2'); xlabel('Latent dimensions');
            if N>=3, subplot(2,3,3); hintonDiagram(Z{3}); title('Mode-3'); end
        case 3,
            subplot(2,3,1); imagesc(Z{1}); title('Mode-1'); ylabel('Length of #-mode');
            subplot(2,3,2); imagesc(Z{2}); title('Mode-2'); xlabel('Latent dimensions');
            if N>=3, subplot(2,3,3); imagesc(Z{3}); title('Mode-3');end
    end
    subplot(2,3,4); bar(gammas); title('Posterior mean of \lambda'); xlabel('Latent components'); ylabel(''); axis tight;
    subplot(2,3,5); plot(LB, '-r.','LineWidth',1.5,'MarkerSize',10 ); title('Lower bound'); xlabel('Iteration');  grid on;
    subplot(2,3,6); plotGamma(a_beta0, a_beta0); title('Posterior pdf'); xlabel('Noise precision \tau');grid on;
    set(findall(h1,'type','text'),'fontSize',12);
    drawnow;
end
if verbose ==2;
    h3 = figure;
    temp = 255.*(X-min(X(:)))/(max(X(:))-min(X(:)));
    imshow(uint8(temp));
    title(['Iter.= '  num2str(0),',  Rank = ' num2str(R)],'FontSize', 13, 'color','b');
    tic;
    xlabel(['(BCPF-MP)  Time: ' num2str(round(toc)), ' seconds'],'FontSize', 13, 'color','b');
    drawnow;
end

%% Model learning
for it=1:maxiters,
    %% Update factor matrices
    Aw = diag(gammas);
    for n=1:N
        % compute E(Z_{\n}^{T} Z_{\n})
        ENZZT = reshape(khatrirao_fast(EZZT{[1:n-1, n+1:N]},'r')' * double(tenmat(O,n)'), [R,R,dimY(n)]);
        % compute E(Z_{\n})
        FslashY = khatrirao_fast(Z{[1:n-1, n+1:N]},'r')' * tenmat(Y.*O, n)';
        for i=1:dimY(n)
            ZSigma{n}(:,:,i) = (beta * ENZZT(:,:,i) + Aw )^(-1);
            Z{n}(i,:) = (beta * ZSigma{n}(:,:,i) * FslashY(:,i))';
        end
        if n<3
            Z{n} = Z{n} + W1{n}*Z{n};
        end
        EZZT{n} = (reshape(ZSigma{n}, [R*R, dimY(n)]) + khatrirao_fast(Z{n}',Z{n}'))';
    end
    
    %% Update latent tensor X
    X = double(ktensor(Z));
    
    %% Update hyperparameters gamma
    a_gammaN = (0.5*sum(dimY) + a_gamma0)*ones(R,1);
    b_gammaN = 0;
    for n=1:N
        b_gammaN = b_gammaN + diag(Z{n}'*Z{n}) + diag(sum(ZSigma{n},3));
    end
    b_gammaN = b_gamma0 + 0.5.* b_gammaN;
    gammas = a_gammaN./b_gammaN;
    
    %% update noise beta
    %  The most time and space consuming part
    if 0 % save time but large space needed
        EX2 =  O(:)' * khatrirao_fast(EZZT,'r') * ones(R*R,1);
    else  % save space but slow
        temp1 = cell(N,1);
        EX2 =0;
        for i =1:R
            for n=1:N
                temp1{n} = EZZT{n}(:,(i-1)*R+1: i*R);
            end
            EX2 = EX2 + O(:)' * khatrirao_fast(temp1,'r')* ones(R,1);
        end
    end
    err = Y(:)'*Y(:) - 2*Y(:)'*X(:) + EX2;
    if  strcmp(noise,'on')
        a_betaN = a_beta0 + 0.5*nObs;
        b_betaN = b_beta0 + 0.5*err;
    else
        a_betaN = a_beta0;
        b_betaN = b_beta0;
    end
    beta = a_betaN/b_betaN;
    Fit = 1 - sqrt(sum(err(:)))/norm(Y(:));
    
    %% Lower bound
    temp1 = -0.5*nObs*safelog(2*pi) + 0.5*nObs*(psi(a_betaN)-safelog(b_betaN)) - 0.5*(a_betaN/b_betaN)*err;
    temp22 =0;
    for n=1:N
        temp22= temp22 + Z{n}'*Z{n} + sum(ZSigma{n},3);
    end
    temp2 = -0.5*R*sum(dimY)*safelog(2*pi) + 0.5*sum(dimY)*sum(psi(a_gammaN)-safelog(b_gammaN)) -0.5*trace(diag(gammas)* temp22);
    temp3 = sum(-safelog(gamma(a_gamma0)) + a_gamma0*safelog(b_gamma0) -  b_gamma0.*(a_gammaN./b_gammaN) + (a_gamma0-1).*(psi(a_gammaN)-safelog(b_gammaN)));
    temp4 = -safelog(gamma(a_beta0)) + a_beta0*safelog(b_beta0) + (a_beta0-1)*(psi(a_betaN)-safelog(b_betaN)) - b_beta0*(a_betaN/b_betaN);
    temp5=0;
    for n=1:N
        for i=1:size(ZSigma{n},3)
            temp5 = temp5 + 0.5*safelog(det(ZSigma{n}(:,:,i))) + 0.5*R*(1+safelog(2*pi));
        end
    end
    temp6 = sum(safelog(gamma(a_gammaN)) - (a_gammaN-1).*psi(a_gammaN) -safelog(b_gammaN) + a_gammaN);
    temp7 = safelog(gamma(a_betaN)) - (a_betaN-1)*psi(a_betaN) -safelog(b_betaN) + a_betaN;
    LB(it) = temp1 + temp2 + temp3 + temp4 + temp5 + temp6 + temp7;
    
    %% Prune irrelevant dimensions?
    Zall = cell2mat(Z);
    comPower = diag(Zall' * Zall);
    comTol = sum(dimY)*eps(norm(Zall,'fro'));
    rankest = sum(comPower> comTol );
    if max(rankest)==0
        disp('Rank becomes 0 !!!');
        break;
    end
    if DIMRED==1  && it >=2,
        if R~= max(rankest)
            indices = comPower > comTol;
            gammas = gammas(indices);
            temp = ones(R,R);
            temp(indices,indices) = 0;
            temp = temp(:);
            for n=1:N
                Z{n} = Z{n}(:,indices);
                ZSigma{n} = ZSigma{n}(indices,indices,:);
                EZZT{n} = EZZT{n}(:, temp == 0);
            end
            R = max(rankest);
        end
    end
    
    %% Display progress
    if it>2
        LBRelChan = abs(LB(it) - 2*LB(it-1) + LB(it-2))/-LB(2);
    else
        LBRelChan = NaN;
    end
    if verbose,
        fprintf('Iter. %d: RelChan = %g, Fit = %g, R = %d \n', it, LBRelChan, Fit, rankest);
    end
    
    %% visualize online results
    if verbose >2 ,
        switch verbose,
            case 4,
                set(0,'CurrentFigure',h1);
                subplot(2,3,1); hintonDiagram(Z{1}); title('Mode-1'); ylabel('Length of #-mode');
                subplot(2,3,2); hintonDiagram(Z{2}); title('Mode-2'); xlabel('Latent dimensions');
                if N>=3, subplot(2,3,3); hintonDiagram(Z{3}); title('Mode-3'); end
            case 3,
                set(0,'CurrentFigure',h1);
                subplot(2,3,1); imagesc(Z{1}); title('Mode-1'); ylabel('Length of #-mode');
                subplot(2,3,2); imagesc(Z{2}); title('Mode-2'); xlabel('Latent dimensions');
                if N>=3, subplot(2,3,3); imagesc(Z{3}); title('Mode-3'); end
        end
        subplot(2,3,4); bar(gammas); title('Posterior mean of \lambda'); xlabel('Latent components'); ylabel(''); axis tight;
        subplot(2,3,5); plot(LB, '-r.','LineWidth',1.5,'MarkerSize',10 ); title('Lower bound'); xlabel('Iteration');  grid on;
        subplot(2,3,6); plotGamma(a_betaN, b_betaN); title('Posterior pdf'); xlabel('Noise precision \tau');grid on;
        set(findall(h1,'type','text'),'fontSize',12);
        drawnow;
    end
    if verbose==2
        set(0,'CurrentFigure',h3);
        figure(h3);
        %        temp = (X-min(X(:)))/(max(X(:))-min(X(:)));
        %        image(temp);
        %        axis off;
        imshow(uint8(X));
        title(['Iter.= '  num2str(it),',  Rank = ' num2str(max(rankest))],'FontSize', 13, 'color','b');
        xlabel(['(BCPF-MP)  Time: ' num2str(round(toc)), ' seconds'],'FontSize', 13, 'color','b');
        drawnow;
    end
    
    %% Convergence check
    if it>5 && abs(LBRelChan) < tol
        disp('\\\======= Converged===========\\\');
        break;
    end
end

%% Predictive distribution
if predVar==1
    Xvar =  tenzeros(size(Y));
    for n=1:N
        Xvar = tenmat(Xvar,n);
        Fslash = khatrirao_fast(Z{[1:n-1, n+1:N]},'r');
        if 1
            temp1 = double(tenmat(tensor(ZSigma{n}),3));
            temp2 = khatrirao_fast(Fslash', Fslash');
            Xvar(:,:) = Xvar(:,:) + temp1*temp2;
        else
            % ---  slow computation ------
            for i=1:size(Xvar,1)     %#ok
                Xvar(i,:) = Xvar(i,:) + diag(Fslash * ZSigma{n}(:,:,i) *Fslash')';
            end
            % ---  slow computation ------
        end
        Xvar = tensor(Xvar);
    end
    Xvar = Xvar + beta^(-1);
    Xvar = Xvar.*(2*a_betaN)/(2*a_betaN-2);
    Xvar = Xvar.*(dscale^2);
else
    Xvar =[];
end

%% Prepare the results
SNR = 10*log10(var(X(:))*beta);
X = ktensor(Z)*dscale;
X = arrange(X);

%% Output
model.X = X;
model.ZSigma = ZSigma;
model.gammas = gammas;
model.Fit = Fit;
model.SNR = SNR;
model.Xvar = double(Xvar);
model.TrueRank = rankest;
model.LowBound = max(LB);


function y = safelog(x)
x(x<1e-300)=1e-200;
x(x>1e300)=1e300;
y=log(x);

