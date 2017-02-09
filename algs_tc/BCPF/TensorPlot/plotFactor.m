function  perf = plotFactor(Ztruth, varargin)

% For example:
%      Ztruth = cell(1,3);
%      Ztruth{1} = randn(10,3);
%      plotFactor(Ztruth, Zest1, Zest2,... )

nMethod = length(varargin);
Zest = varargin;

for i=1:nMethod
    [Ztruth,Zest{i}] = rescaleorder(Ztruth, Zest{i});
end

N = length(Ztruth);
nCom = size(Ztruth{1},2);

perf = zeros(nCom, N);
for n=1:N
    perf(:,n) = diag(corr(Zest{1}{n},Ztruth{n}));
end



%% plot the factors
figure;
str = {'-r*','-bs','-go','-y+'};
k=1;
for r=1:nCom
    for n =1:N
        subplot(nCom, N,k);
        plot(Ztruth{n}(:,r),str{1});
        hold on;
        for m=1:nMethod
            plot(Zest{m}{n}(:,r),str{m+1});
        end
        
        axis tight; hold off;
        
        if k==1, title('Mode #1'); ylabel('Com. #1'); end
        if k==2, title('Mode #2'); end
        if k==3, title('Mode #3'); end
        if k==4, ylabel('Com. #2'); end
        if k==7, ylabel('Com. #3'); end
        
        k = k+1;
    end
end
legend({'True factors', 'Inferred factors'});

%% plot the factors
% figure;
% k=1;
% for n=1:N
%     for m =0:nMethod
%         subplot(N, nMethod+1,k);
%         if m==0
%            plot(Ztruth{n}); ylabel(['Mode-' num2str(n)]); title('Ground truth');
%            axis tight; box on;
%         else
%            plot(Zest{m}{n});
%             axis tight; box on;
%         end
%         k = k+1;
%     end
% end









