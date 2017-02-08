%%% void show_4dvideo(4dmatrix)
%
function show_4dvideo(video)
  for i = 1 : size(video,4)
    if(mod(i,10) == 0) fprintf('%d ',i); end
    frame = video(:,:,:,i);
    imshow(frame,[],'InitialMagnification','fit');
    pause(0.01);
  end
  disp(i);
end
