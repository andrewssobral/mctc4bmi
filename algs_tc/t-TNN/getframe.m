function frame = getframe(y, sY, i)

temp = reshape(y, sY);
frame = temp(:,:,i);
