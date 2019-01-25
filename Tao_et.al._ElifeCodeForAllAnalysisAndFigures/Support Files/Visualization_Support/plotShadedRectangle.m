function [] = plotShadedRectangle(x,y,c1,c2)

ROIx = [x(1) x(2) x(2) x(1)];
ROIy = [y(1) y(1) y(2) y(2)];

% Add lines
h1 = line([x(1) x(1)],[y(1) y(2)]);
h2 = line([x(2) x(2)],[y(1) y(2)]);
% Set properties of lines
set([h1 h2],'Color',c1,'LineStyle','none')
% Add a patch
patch(ROIx,ROIy,c2,'EdgeColor','none')
% The order of the "children" of the plot determines which one appears on top.
% I need to flip it here.
set(gca,'children',flipud(get(gca,'children')))
end