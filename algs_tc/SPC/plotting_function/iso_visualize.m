function [p] = iso_visualize(Z,iso);

[d1 d2 d3] = size(Z);
[x, y, z] = meshgrid (1:d1, 1:d2, 1:d3);

view (120, 30);
p = patch(isosurface (x, y, z, Z, iso));
isonormals(x,y,z, Z, p)
set(p, 'FaceColor', 'green', 'EdgeColor', 'none');
daspect([1 1 1]); axis tight; 
colormap(prism(28))
camlight; lighting phong

