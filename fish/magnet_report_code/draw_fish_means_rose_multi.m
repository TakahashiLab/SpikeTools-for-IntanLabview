function fileMap = draw_fish_means_rose_multi(S, conditions_in, outdir, ws, wl, pngRes, drawMeanVector, uniformTbl, alpha, tiles, writePNG)
fileMap = struct(); rows = tiles(1); cols = tiles(2);

pdfPath = fullfile(outdir, 'rose_fishMeans_all.pdf');
pngPath = fullfile(outdir, 'rose_fishMeans_all.png');
if exist(pdfPath,'file'), try, delete(pdfPath); end, end
if exist(pngPath,'file'), try, delete(pngPath); end, end

C = cellstr(string(conditions_in)); N = numel(C);
pageCount = 0; nPerPage = rows*cols;

for j = 1:N
    if mod(j-1, nPerPage)==0
        if exist('f','var'), finalize_page(f, pdfPath, pngPath, pngRes, pageCount>0, writePNG); close(f); end
        f = figure('Color','w','Position',[100 100 360*cols 360*rows]); %#ok<NASGU>
        tlo = tiledlayout(rows, cols, 'Padding','compact','TileSpacing','compact');
        pageCount = pageCount + 1;
        safe_sgtitle(tlo, sprintf('Fish-means Rose (Window %gs-%gs) page %d', ws, finite_end(ws, wl), pageCount));
    end

    cj = C{j}; Sj = S(S.Condition==string(cj), :);
    tileIdx = mod(j-1, nPerPage) + 1;
    ax = polaraxes(tlo); ax.Layout.Tile = tileIdx; hold(ax,'on');

    if isempty(Sj)
        title(ax, sprintf('%s (No data)', safe_text(cj)), 'Interpreter','none');
        continue;
    end

    G    = findgroups(Sj.Pair);
    cosm = splitapply(@(x) mean(x,'omitnan'), Sj.cosMag, G);
    sinm = splitapply(@(x) mean(x,'omitnan'), Sj.sinMag, G);
    ok   = ~isnan(cosm) & ~isnan(sinm);
    alpha_i = atan2(sinm(ok), cosm(ok));

    if ~isempty(alpha_i)
        r = ones(size(alpha_i));
        polarplot(ax, alpha_i, r, 'o', 'MarkerSize',5, 'LineStyle','none');
        rticks(ax, [0 0.5 1]); ax.RLim=[0 1];
        if drawMeanVector
            Cbar = mean(cos(alpha_i)); Sbar = mean(sin(alpha_i));
            mu = atan2(Sbar, Cbar);
            polarplot(ax, [mu mu], [0 1], 'LineWidth', 2);
        end
    else
        text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center','FontSize',10);
    end

    pv = NaN;
    if ~isempty(uniformTbl)
        ii = find(uniformTbl.Condition==string(cj),1);
        if ~isempty(ii), pv = uniformTbl.p(ii); end
    end
    ttl = sprintf('Fish-means: %s (p=%.3g)', safe_text(cj), pv);
    if ~isnan(pv) && pv < alpha, ttl = [ttl ' *']; end
    title(ax, ttl, 'Interpreter','none');
    ax.ThetaZeroLocation='top'; ax.ThetaDir='clockwise';
end

if exist('f','var')
    finalize_page(f, pdfPath, pngPath, pngRes, pageCount>0, writePNG); close(f);
end

fileMap.all = struct('pdf', pdfPath, 'png', ternary(writePNG, pngPath, ''));
end
