function t = Yaxis_string(Xt,bnames,fs);

%ax = axis;
%set(gca,'YLimMode','manual')
%Yl = ax(1:2);
%t = text(Yl(1)*ones(1,length(Xt)),Xt,bnames);
%set(t,'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',fs);
set(gca,'YDir','normal');
set(gca,'YTick',Xt,'YTickLabel',bnames,'fontsize',fs)
t = gca;
%t=1;
