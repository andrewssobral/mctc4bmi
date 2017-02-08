%%% [3dmatrix] = convert_video4d_to_3d(4dvideo,string)
% format = '' | 'uint8' | 'double'
function [V] = convert_video4d_to_3d(video,format)
  if(nargin < 2) format = ''; end
  
  [vidHeight,vidWidth,~,nFrames] = size(video);
  
  if(strcmp(format,'uint8'))
    V = uint8(zeros(vidHeight,vidWidth,nFrames));
  end
  if(strcmp(format,'double'))
    V = zeros(vidHeight,vidWidth,nFrames);
  end
  
  for i = 1 : nFrames
    im = video(:,:,:,i); % video(m,n,p,i)
    imgray = rgb2gray(im);
    if(strcmp(format,'double'))
      imgray = im2double(imgray);
    end
    V(:,:,i) = imgray;
  end
end
