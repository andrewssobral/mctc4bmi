function [algs] = get_algs_name(params)
  algs = struct('name', {});
  list = dir(params.algs_path);
  j = 1;
  for i = 1:size(list,1)
    if(strcmp(list(i).name,'.') || strcmp(list(i).name,'..'))
      continue;
    else
      algs(j).name = list(i).name;
      j = j + 1;
    end
  end
end

