function V = load_images(params)
  sequence_name = params.sequences_name(params.current_sequence).name;
  path = fullfile(params.sequences_path,sequence_name,params.sequences_format);
  list = dir(fullfile(path,params.sequences_ext));
  names = char({list.name}');
  frames = size(names,1);
  V = [];
  for i = 1:frames
    if(mod(i,10) == 0) fprintf('%d ',i); end
    I = imread(fullfile(path,names(i,:)));
    I = im2double(I);
    V(:,:,:,i) = I;
  end
  disp(i);
end
