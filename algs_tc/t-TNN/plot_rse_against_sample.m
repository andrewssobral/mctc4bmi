clear all
close all
clc
addpath('Data','MyLib','ShowImg','Solvers');
addpath(genpath('Results'));
Video = {'basketball','people','mountainbike','airport','escalator','fountain'};
sr = 0.1:0.1:0.9;
for ii = 1 : 6
    VideoName = [Video{ii} '_RSE.mat'];
    load(VideoName); 
%     Temp = VRSE_AtSVD - VRSE_LtSVD;
%     Mask = (Temp>0);
%     VRSE_LtSVD(Mask) = VRSE_AtSVD(Mask);
    h = figure('Name',[Video{ii} ' RSE against sample rate']);
    plot(sr,-VRSE_LtSVD,'-r','linewidth',2.5,'Marker','o','MarkerSize',8.5);hold on;
    plot(sr,-VRSE_FtSVD,'--b','linewidth',2.5,'Marker','d','MarkerSize',8.5);hold on;    
    plot(sr,-VRSE_SNN,'-.g','linewidth',2.5,'Marker','^','MarkerSize',8.5);hold on;
    plot(sr,-VRSE_MNN,'-.c','linewidth',2.5,'Marker','>','MarkerSize',8.5);hold on;
    plot(sr,-VRSE_TMac,'--m','linewidth',2.5,'Marker','s','MarkerSize',8.5);hold on;
    plot(sr,-VRSE_NN,':k','linewidth',2.5,'Marker','v','MarkerSize',8.5);hold on;
%     title('RSE Agaist Sample Rate');
    grid on;
    xlabel('Sample Rate');
    ylabel('RSE(dB)');
    legend('t-TNN','GTNN','SNN','LNN','TMac','MNN','Location','NorthEast');
    saveas(h,[Video{ii} '_rse_against_sr'], 'fig');
end    