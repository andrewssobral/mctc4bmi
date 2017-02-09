function CreateMovie(T, fileName)
[height, width, layer] = size(T{1});
aviobj = avifile(fileName, 'fps', 30, 'compression', 'None');
for i = 1:layer
    img(:,:,1) = T{1}(:,:,i);
    img(:,:,2) = T{2}(:,:,i);
    img(:,:,3) = T{3}(:,:,i);
    aviobj = addframe(aviobj, img);
end
aviobj = close(aviobj);
