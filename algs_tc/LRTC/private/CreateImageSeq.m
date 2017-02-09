function CreateImageSeq(inData, filePath, format)
[height, width, layer] = size(inData{1});
for n = 1:layer
    Img(:,:,1) = inData{1}(:,:,n);
    Img(:,:,2) = inData{2}(:,:,n);
    Img(:,:,3) = inData{3}(:,:,n);
    fileName = [filePath, int2str(n), format];
    imwrite(Img, fileName);
end