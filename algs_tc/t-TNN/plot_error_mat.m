clear all
close all
clc

h1 = figure(1);

% load numerical_results_50_50_40
load numerical_results
% ErrorMat = ErrorMat./max(abs(ErrorMat(:)));
Temp = ErrorMat(:,1:end-1);
ErrorMat(:,end) = max(Temp(:));
% ErrorMat  = ErrorMat/max(ErrorMat(:));
ErrorMat(ErrorMat<0) = 0;
LogErrorMat = log(ErrorMat);
LogErrorMat(LogErrorMat<0) = 0;
LogErrorMat = LogErrorMat/max(LogErrorMat(:));

ps = [0.01:0.01:1];
pr = [0.01:0.01:1];

% subplot(231);
% imagesc(ps,pr,ErrorMat, [-.1,1.2]); 
imagesc(ps,pr,LogErrorMat, [0,1]); 
figure(gcf);colormap(gray(256)); colorbar;
set(gca, 'XTick', [0.1:0.2:1]);
set(gca, 'YTick', [0.1:0.2:1]);
axis xy; 
xlabel('p','fontsize',20);
ylabel('r','fontsize',20);

print(h1,'-djpeg','-r200','recovery1.jpg')

h2 = figure(2);
ErrorMat(ErrorMat<30) = 0;
ErrorMat(ErrorMat>=30) = 1;
ErrorMat = ErrorMat/max(ErrorMat(:));
imagesc(ps,pr,ErrorMat, [0,1]); 
figure(gcf);colormap(gray(256)); colorbar;
set(gca, 'XTick', [0.1:0.2:1]);
set(gca, 'YTick', [0.1:0.2:1]);
axis xy; 
xlabel('p','fontsize',20);
ylabel('r','fontsize',20);

print(h2,'-djpeg','-r200','recovery2.jpg')


