function draw_single_raincloud(ax, y)
y = y(:); y = y(~isnan(y));
if isempty(y), return; end
ymin = min(y); ymax = max(y); yr = ymax-ymin; if yr<=0, yr=1; end
ygrid = linspace(ymin-0.05*yr, ymax+0.05*yr, 200);
try
    [fhat, yvals] = ksdensity(y, ygrid);
catch
    [counts, edges] = histcounts(y, max(10, round(sqrt(numel(y)))));
    yvals = movmean(edges(1:end-1) + diff(edges)/2, 2, 'Endpoints','shrink');
    fhat  = movmean(counts, 3, 'Endpoints','shrink');
end
if all(~isfinite(fhat)) || max(fhat)<=0
    fhat = ones(size(ygrid)); yvals = ygrid;
end
fhat = fhat ./ max(fhat);
w = fhat * 0.35;
x0 = 0;
patch(ax, [x0-w, fliplr(x0+w)], [yvals, fliplr(yvals)], [0.7 0.7 0.7], ...
      'FaceAlpha', 0.25, 'EdgeColor','none');
xj = (rand(numel(y),1)-0.5)*0.10;
plot(ax, xj, y, 'o', 'MarkerSize',4, 'LineWidth',0.5);
ymed = median(y,'omitnan'); plot(ax, [-0.20 0.20], [ymed ymed], '-', 'LineWidth',2.5);
xlim(ax, [-0.5 0.5]); ax.YLim = [min([-180,ymin]) max([180,ymax])];
ax.XTick = []; box(ax,'on');
end
