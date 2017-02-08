function [sequences] = get_sequences_name(params)
  sequences = struct('name', {});
  list = dir(params.sequences_path);
  j = 1;
  for i = 1:size(list,1)
    if(strcmp(list(i).name,'.') || strcmp(list(i).name,'..'))
      continue;
    else
      sequences(j).name = list(i).name;
      j = j + 1;
    end
  end
end

