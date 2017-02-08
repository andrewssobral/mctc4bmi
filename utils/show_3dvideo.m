%%% void show_3dvideo(3dmatrix)
%
function show_3dvideo(V)
  for i = 1 : size(V,3)
    if(mod(i,10) == 0) fprintf('%d ',i); end
    frame = V(:,:,i);
    imshow(frame,[],'InitialMagnification','fit');
    %imshow(frame,'InitialMagnification','fit');
    pause(0.01);
  end
  disp(i);
end
