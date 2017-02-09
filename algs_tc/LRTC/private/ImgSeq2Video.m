function ImgSeq2Video(preFileName, formatFile, nS, nE, isColor)
% transform the image sequence into the .mat format
if isColor == 0
    for i = nS:nE
        fileName = [preFileName, int2str(i), formatFile];
        data(:,:,i-nS+1) = double(imread(fileName));
    end
    save(['preFileName', 'gray.mat'], 'data');
end

if isColor == 1
    for i = nS:nE
        fileName = [preFileName, int2str(i), formatFile];
        I = double(imread(fileName));
        dataR(:,:,i-nS+1) = I(:,:,1);
        dataG(:,:,i-nS+1) = I(:,:,2);
        dataB(:,:,i-nS+1) = I(:,:,3);
    end
    save([preFileName, 'R.mat'], 'dataR');
    save([preFileName, 'G.mat'], 'dataG');
    save([preFileName, 'B.mat'], 'dataB');
end