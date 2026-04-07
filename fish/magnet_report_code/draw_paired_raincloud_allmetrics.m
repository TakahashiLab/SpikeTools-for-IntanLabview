function fileMap = draw_paired_raincloud_allmetrics(PF, conditions_in, commonPairs, METRICS, outdir, ws, wl, tiles, writePNG, pngRes)
fileMap = struct();
if isempty(tiles), cols=2; rows=ceil(numel(METRICS)/cols); else, rows=tiles(1); cols=tiles(2); end

pdfPath = fullfile(outdir, 'raincloud_stats_allmetrics.pdf');
pngPath = fullfile(outdir, 'raincloud_stats_allmetrics.png');
if exist(pdfPath,'file'), try, delete(pdfPath); end, end
if exist(pngPath,'file'), try, delete(pngPath); end, end

C = string(conditions_in(:)); nC = numel(C);
P = string(commonPairs(:)); nF = numel(P);

f = figure('Color','w','Position',[100 100 480*cols 380*rows]);
tlo = tiledlayout(rows, cols, 'Padding','compact','TileSpacing','compact');
safe_sgtitle(tlo, sprintf('Paired raincloud (Window %gs-%gs)', ws, finite_end(ws, wl)));

for mi = 1:numel(METRICS)
    metric = METRICS{mi};
    ax = nexttile(tlo); hold(ax,'on');

    X = nan(nF, nC);
    for cj = 1:nC
        Tj = PF(strcmp(PF.Condition, C(cj)), {'Pair', metric});
        [tf,loc] = ismember(P, string(Tj.Pair));
        vals = nan(nF,1); if ~isempty(Tj); vals(tf) = Tj.(metric)(loc(tf)); end
        X(:, cj) = vals;
    end

    jitter = 0.06;
    for i = 1:nF
        yi = X(i, :);
        good = ~isnan(yi);
        xi = 1:nC; xi = xi(good); yi = yi(good);
        if numel(xi) >= 2
            plot(ax, xi, yi, '-', 'LineWidth', 0.8, 'Color', [0 0 0 0.25]);
        end
    end

    for cj = 1:nC
        yi = X(:, cj); yi = yi(~isnan(yi));
        if isempty(yi), continue; end
        ygrid = linspace(min(yi), max(yi), 200);
        try
            [fhat, yvals] = ksdensity(yi, ygrid);
        catch
            [counts, edges] = histcounts(yi, max(10, round(sqrt(numel(yi)))));
            yvals = movmean(edges(1:end-1)+diff(edges)/2,2,'Endpoints','shrink');
            fhat  = movmean(counts,3,'Endpoints','shrink');
        end
        if all(~isfinite(fhat)) || max(fhat)<=0, fhat=ones(size(ygrid)); yvals=ygrid; end
        fhat = fhat/max(fhat); w = fhat*0.35;
        patch(ax, [cj-w, fliplr(cj+w)], [yvals, fliplr(yvals)], [0.7 0.7 0.7], 'FaceAlpha',0.25,'EdgeColor','none');
        xj = cj + (rand(numel(yi),1)-0.5)*2*jitter;
        plot(ax, xj, yi, 'o', 'MarkerSize',4, 'LineWidth',0.5);
        ymed = median(yi, 'omitnan'); plot(ax, [cj-0.25 cj+0.25], [ymed ymed], '-', 'LineWidth',2.5);
    end

    xlim(ax, [0.5 nC+0.5]); xticks(ax,1:nC); xticklabels(ax,cellstr(C));
    ax.XTickLabelRotation = 15; grid(ax,'on'); ylabel(ax, metric, 'Interpreter','none');
    title(ax, metric, 'Interpreter','none');

% --- NEW: significance annotations between conditions using paired tests over X ---
try
    yl = ylim(ax); y_top = yl(2); y_step = 0.06*(yl(2)-yl(1));
    ypos = y_top + y_step*0.2; hold(ax,'on');
    pair_idx = nchoosek(1:nC,2); row_offset = 0;
    for pp = 1:size(pair_idx,1)
        i = pair_idx(pp,1); j = pair_idx(pp,2);
        vi = X(:, i); vj = X(:, j); ok = isfinite(vi) & isfinite(vj);
        vi = vi(ok); vj = vj(ok);
        pval = NaN;
        if numel(vi) >= 3 && numel(vj) >= 3
            try, pval = signrank(vi, vj); catch, try, [~,pval] = ttest(vi, vj); end, end
        end
        stars = ''; if ~isnan(pval)
            if pval < 0.001, stars='***'; elseif pval < 0.01, stars='**'; elseif pval < 0.05, stars='*'; end
        end
        if ~isempty(stars)
            x1=i; x2=j;
            plot(ax, [x1 x1 x2 x2], [ypos+row_offset, ypos+row_offset+y_step*0.5, ypos+row_offset+y_step*0.5, ypos+row_offset], '-','LineWidth',1);
            text(ax, (x1+x2)/2, ypos+row_offset+y_step*0.55, stars, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12);
            row_offset = row_offset + y_step*0.9;
        end
    end
    ylim(ax, [yl(1), yl(2) + max(row_offset,0) + y_step*0.5]);
catch
end
% --- END NEW ---

end

try, exportgraphics(f, pdfPath, 'ContentType','image'); catch, print(f, pdfPath, '-dpdf', '-painters'); end
if writePNG
    try, exportgraphics(f, pngPath, 'Resolution', pngRes); catch, print(f, pngPath, sprintf('-r%d', pngRes), '-dpng'); end
end
close(f);

fileMap.pdf = pdfPath; fileMap.png = ternary(writePNG, pngPath, '');
end
