function h = showvideo(X, tp)
nf = size(X,3);
fn = sprintf('Show Video : Pause %.2fs', tp);
h = figure('Name', fn);
for i = 1 : nf
    imshow(X(:,:,i),[]);
    tn = sprintf('The %sth Frame', num2str(i));
    title(tn);
    pause(tp);
end
