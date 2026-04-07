function fileMap = plot_circ_pairwise_overview(S, conditions_in, circPair, outdir, ws, wl, alpha, tiles, writePNG, pngRes)
if nargin<7 || isempty(tiles)
    nPairs = size(nchoosek(1:numel(conditions_in),2),1);
    cols = min(2, nPairs); rows = ceil(nPairs/cols);
else
    rows = tiles(1); cols = tiles(2);
end
if nargin<8, writePNG = false; end
if nargin<9, pngRes = 300; end

pdfPath = fullfile(outdir, 'circ_pairwise_overview.pdf');
pngPath = fullfile(outdir, 'circ_pairwise_overview.png');
if exist(pdfPath,'file'), try, delete(pdfPath); end, end
if exist(pngPath,'file'), try, delete(pngPath); end, end

conds = string(conditions_in(:));
pairsIdx = nchoosek(1:numel(conds),2);
nTiles = size(pairsIdx,1);

if ~ismember('Pair', S.Properties.VariableNames)
    S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
end
condCol = string(S.Condition);

f = figure('Color','w','Position',[100 100 520*cols 380*rows]);
tlo = tiledlayout(rows, cols, 'Padding','compact','TileSpacing','compact');
safe_sgtitle(tlo, sprintf('Circular pairwise Δθ (deg) | Window %gs-%gs', ws, finite_end(ws, wl)));

for ti = 1:nTiles
    ax = nexttile(tlo); hold(ax,'on');
    ia = pairsIdx(ti,1); ib = pairsIdx(ti,2);
    cA = conds(ia); cB = conds(ib);

    d = compute_deltas_rad(S, condCol, cA, cB);
    d = d(~isnan(d));
    if isempty(d)
        title(ax, sprintf('%s vs %s (No data)', safe_text(cA), safe_text(cB)), 'Interpreter','none');
        grid(ax,'on'); xlim(ax,[-180 180]); xlabel(ax,'Δθ (deg)'); continue;
    end

    ddeg = rad2deg(d);
    draw_single_raincloud(ax, ddeg);
    star = '';
    if ~isempty(circPair)
        row = find( (circPair.CondA==cA & circPair.CondB==cB) | (circPair.CondA==cB & circPair.CondB==cA), 1,'first');
        if ~isempty(row)
            pv = circPair.p_adj(row);
            if ~isnan(pv) && pv < alpha, star = '*'; end
        end
    end
    title(ax, sprintf('%s vs %s %s', safe_text(cA), safe_text(cB), star), 'Interpreter','none');
    xlabel(ax, 'Δθ (deg)'); xlim(ax,[-180 180]); grid(ax,'on');

    pax = polaraxes('Position', inset_rect(ax, 0.70, 0.60, 0.25, 0.35)); hold(pax,'on');
    C = mean(cos(d)); Ssin = mean(sin(d));
    mu = atan2(Ssin, C); Rbar = hypot(C,Ssin);
    polarplot(pax, [mu mu], [0 Rbar], 'LineWidth', 2);
    pax.RLim = [0 1]; pax.ThetaZeroLocation='top'; pax.ThetaDir='clockwise';
    title(pax, 'mean Δ', 'FontSize',9);
end

try, exportgraphics(f, pdfPath, 'ContentType','image'); catch, print(f, pdfPath, '-dpdf', '-painters'); end
if writePNG
    try, exportgraphics(f, pngPath, 'Resolution', pngRes); catch, print(f, pngPath, sprintf('-r%d', pngRes), '-dpng'); end
end
close(f);

fileMap = struct('pdf', pdfPath, 'png', ternary(writePNG, pngPath, ''));
end
