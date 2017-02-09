function h = setpapersize(s)

h = 1;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 s(1) s(2)]);


