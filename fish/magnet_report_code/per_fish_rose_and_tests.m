function [rayTbl, fileMap] = per_fish_rose_and_tests(S, conditions_in, outdir, ws, wl, opt)
fileMap = struct();
rayC = strings(0,1); rayPair = strings(0,1); rayN = []; rayMean = []; rayR = []; rayZ = []; rayP = [];

tiles = opt.perFishTiles; rows = tiles(1); cols = tiles(2);
binEdges = linspace(-pi, pi, opt.roseBinN+1);

for j = 1:numel(conditions_in)
    cj = conditions_in(j);
    Sj = S(S.Condition==cj, :);
    if isempty(Sj), continue; end
    pairs = unique(Sj.Pair);
    if isempty(pairs), continue; end

    pdfPath = fullfile(outdir, sprintf('rose_perFish_%s.pdf', safe_text(cj)));
    pngPath = fullfile(outdir, sprintf('rose_perFish_%s.png', safe_text(cj)));
    if exist(pdfPath,'file'), try, delete(pdfPath); end, end
    if exist(pngPath,'file'), try, delete(pngPath); end, end
    pageCount = 0;

    for k = 1:numel(pairs)
        if mod(k-1, rows*cols)==0
            if exist('f','var'), finalize_page(f, pdfPath, pngPath, opt.pdfRasterRes, pageCount>0, opt.perFishWritePNG); close(f); end
            f = figure('Color','w','Position',[100 100 360*cols 360*rows]); %#ok<NASGU>
            tlo = tiledlayout(rows, cols, 'Padding','compact','TileSpacing','compact');
            pageCount = pageCount + 1;
            safe_sgtitle(tlo, sprintf('Per-fish Rose: %s (Window %gs-%gs) page %d', ...
                safe_text(cj), ws, finite_end(ws, wl), pageCount));
        end
        tileIdx = mod(k-1, rows*cols) + 1;
        ax = polaraxes(tlo); ax.Layout.Tile = tileIdx; hold(ax,'on');

        p  = pairs(k);
        Sk = Sj(Sj.Pair==p, :);
        ang = deg2rad(Sk.HeadingMag); ang = ang(~isnan(ang));

        rmax = 1;
        if ~isempty(ang)
            h = polarhistogram(ax, ang, binEdges);
            if ~isempty(h.Values), rmax = max(1, max(h.Values)); end
        else
            text(ax,0.5,0.5,'No data','Units','normalized','HorizontalAlignment','center','FontSize',10);
        end

        if opt.drawMeanVector && ~isempty(ang)
            C = mean(cos(ang)); Ssin = mean(sin(ang));
            mu = atan2(Ssin, C);
            polarplot(ax, [mu mu], [0 rmax], 'LineWidth', 2);
        end

        if opt.doPerFishRayleigh
            [pval, z, Rbar, ~, n] = rayleigh_on_angles(ang);
            if opt.annotatePerFishP
                txt = sprintf('N=%d, R=%.3f, Z=%.2f, p=%.3g', n, Rbar, z, pval);
                annotation_on_axes(ax, txt);
            end
            if ~isnan(pval) && pval < opt.alpha
                add_sig_star(ax);
            end
            rayC(end+1,1)    = cj;                 
            rayPair(end+1,1) = string(p);          
            rayN(end+1,1)    = n;                  
            rayMean(end+1,1) = NaN;                
            rayR(end+1,1)    = Rbar;               
            rayZ(end+1,1)    = z;                  
            rayP(end+1,1)    = pval;               
        end

        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir          = 'clockwise';
        title(ax, sprintf('%s', safe_text(p)), 'Interpreter','none');
    end

    if exist('f','var')
        finalize_page(f, pdfPath, pngPath, opt.pdfRasterRes, pageCount>0, opt.perFishWritePNG); close(f);
        clear f;
    end

    fileMap.(safe_field(cj)) = struct('pdf', pdfPath, 'png', ternary(opt.perFishWritePNG, pngPath, ''));
end

rayTbl = table(rayC, rayPair, rayN, rayMean, rayR, rayZ, rayP, ...
    'VariableNames', {'Condition','Pair','N','mean_angle_deg','Rbar','Z','p'});
end
