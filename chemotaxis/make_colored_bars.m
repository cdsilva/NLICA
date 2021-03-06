function make_colored_bars(bar_height, bar_color)
% bar_height is the vector of bar values
% bar_color is the vector of bar colors (assumed to be between 0 and 1)

for i=1:length(bar_height)
    bar(i, bar_height(i), 'facecolor', (1-min(bar_color(i),1))*[1 1 0]+[0 0 1])
    hold on
end

colormap([linspace(1, 0, 64)' linspace(1, 0, 64)' ones(64, 1)])
cbar = colorbar('peer',gca);
%set(get(cbar,'ylabel'),'String','CV error');
set(get(cbar,'xlabel'),'String','r_k');
caxis([0 1])
set(gca, 'xtick', 1:length(bar_height))