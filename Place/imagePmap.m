function imagePmap(rate_map,oc_map)
axis equal off;
%ratemap
clims = [0 max(rate_map(:))];
no_occupancy=(oc_map==0);
[cbar, clims] = SmartColorbar(clims, 'jet(255)');
rate_map(no_occupancy) = clims(1);
imagesc(rate_map, clims); 
colormap(cbar);
axis equal off
peak=max(rate_map(:));
gTitle=sprintf('Peak:%1.2f',peak);
title(gTitle);

return;