clear all
close all
clc

addpath('Data','MyLib','ShowImg','Solvers');
addpath(genpath('Results/.'));
Video = {'basketball','people','mountainbike','escalator','building','windmill','snow','sailingship','led'};
for ii = 1:9
VideoName = [Video{ii} '_Rec_SR_0_3.mat'];
load(VideoName); 
Tn = T./max(T(:));
% X_AtSVD = reshape(X_AtSVD,size(T));
% X_LtSVD = reshape(X_LtSVD,size(T));
% X_FtSVD = reshape(X_FtSVD,size(T));
[n1,n2,n3] = size(T);
frse_ltsvd = zeros(n3,1);
frse_ftsvd = zeros(n3,1);
frse_snn = zeros(n3,1);
frse_mnn = zeros(n3,1);
frse_tmac = zeros(n3,1);
frse_nn = zeros(n3,1);

for i = 1 : n3
    frse_ltsvd(i) = -rse(X_LtSVD(:,:,i),Tn(:,:,i));
    frse_ftsvd(i) = -rse(X_FtSVD(:,:,i),Tn(:,:,i)); 
    frse_snn(i) = -rse(X_SNN(:,:,i),Tn(:,:,i)*255);
    frse_mnn(i) = -rse(X_MNN(:,:,i),Tn(:,:,i)*255);
    frse_tmac(i) = -rse(X_TMac(:,:,i),Tn(:,:,i));
    frse_nn(i) = -rse(X_NN(:,:,i),Tn(:,:,i)*255); 
end

h = figure('Name',[Video{ii} ' RSE OF EACH FRAME']);
plot(1:n3,frse_ltsvd,'-r','linewidth',3);hold on;
plot(1:n3,frse_ftsvd,'--b','linewidth',3);hold on;
plot(1:n3,frse_snn,'-.g','linewidth',3);hold on;
plot(1:n3,frse_mnn,'-.','color',[1 0.5 0], 'linewidth',3);hold on;
plot(1:n3,frse_tmac,'--m','linewidth',3);hold on;
plot(1:n3,frse_nn,':k','linewidth',3);hold on;
%     title('RSE Agaist Sample Rate');
grid on;
xlabel('Frame');
ylabel('RSE(dB)');
% legend('NtSVD','GtSVD','TMac','SNN','NN','Location','NorthEast');
saveas(h,[Video{ii} '_rse_of_each_frame'], 'fig');
end
