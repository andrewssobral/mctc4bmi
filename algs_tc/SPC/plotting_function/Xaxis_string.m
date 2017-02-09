function t = Xaxis_string(Xt,bnames,R,fs);

ax = axis;
set(gca,'XLimMode','manual')
Yl = ax(3:4);
t = text(Xt,Yl(1)*ones(1,length(Xt)),bnames);
if R == 0
  set(t,'HorizontalAlignment','center','VerticalAlignment','top', 'Rotation',R,'fontsize',fs);
elseif R == 90
  set(t,'HorizontalAlignment','right','VerticalAlignment','middle', 'Rotation',R,'fontsize',fs);
else
  set(t,'HorizontalAlignment','right','VerticalAlignment','top', 'Rotation',R,'fontsize',fs);
end
set(gca,'XTickLabel','')

