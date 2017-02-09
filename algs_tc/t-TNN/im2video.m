function vi = im2video(file_path, img_format, issave, vi_name, isshow)

if ~exist('issave', 'var')
    issave = 1;
    vi_name = 'temp_name';
end
vi_name = [vi_name '.mat'];
if ~exist('ishow', 'var')
    isshow = 1;
end

img_format = ['*.' img_format];
file_path = ['.\' file_path '\'];
img_path_list = dir(strcat(file_path, img_format));
img_num = length(img_path_list);
img_name = img_path_list(1).name;
img = imread(strcat(file_path, img_name));
[n, m, c] = size(img);
vi = zeros(n,m,img_num);
if img_num > 0
    for i = 1 : img_num
        img_name = img_path_list(i).name;
        img = imread(strcat(file_path, img_name));
        if c == 3
            img = rgb2gray(img);
        end
        img = double(img);
        if max(img(:))>5
            img = img/255;
        end
        vi(:,:,i) = img;
    end
end

if issave == 1
    save(vi_name, 'vi');
end

if isshow == 1
    showvideo(vi, 0.1);
end
       