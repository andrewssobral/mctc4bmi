function  h = plotYXS(varargin)

nTensor = length(varargin);
% scrsz = get(0,'ScreenSize');
h = figure('Position',[100 250  200*nTensor 200]);

strName = cell(1,4);
strName{1} = 'Incomplete Tensor';
strName{2} = 'Low Rank Approximation';
strName{3} = 'Sparse';
strName{4} = 'Noise';

parThresh = zeros(1,nTensor);
parThresh(3) =0.2;

for i = 1:nTensor
    subplot(1,nTensor,i);
    voxel3(abs(varargin{i}),'thresh',parThresh(i), 'degree',5);
    title(strName{i},'Fontsize',11,'FontWeight','bold');
    xlabel(''); ylabel(''); zlabel('');
    set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
    box on;
end

if nTensor==4
    annotation(h,'textbox',[0.29 0.45 0.045 0.165],...
        'String',{'='}, 'FontSize', 20, 'LineStyle', 'none', ...
        'FitBoxToText','on');
    annotation(h,'textbox',[0.49 0.45 0.045 0.165],...
        'String',{'+'}, 'FontSize', 20, 'LineStyle', 'none', ...
        'FitBoxToText','on');
    annotation(h,'textbox',[0.7 0.45 0.045 0.165],...
        'String',{'+'}, 'FontSize', 20, 'LineStyle', 'none', ...
        'FitBoxToText','on');
end