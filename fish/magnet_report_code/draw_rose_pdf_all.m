function draw_rose_pdf_all(S, conditions_in, out_pdf, out_png, ws, wl, pngRes, drawMeanVector, drawTargetOnRose, targetDirDeg, binN, writePNG)
if isempty(conditions_in), return; end
conditions_in = cellstr(string(conditions_in));
N = numel(conditions_in);
f = figure('Color','w','Position',[100 100 max(360*N,500) 360]);
tlo = tiledlayout(1, N, 'Padding','compact','TileSpacing','compact');
binEdges = linspace(-pi, pi, binN+1);
theta0 = normalize_target_angles(targetDirDeg, N);

for j = 1:N
    pax = polaraxes(tlo); pax.Layout.Tile = j; hold(pax,'on');
    condMask = strcmpi(string(S.Condition), conditions_in{j});
    ang = deg2rad(S.HeadingMag(condMask)); ang = ang(~isnan(ang));
    h = []; rmax = 1;
    if ~isempty(ang)
        h = polarhistogram(pax, ang, binEdges);
        if ~isempty(h.Values), rmax = max(1, max(h.Values)); end
    else
        text(pax, 0.5, 0.5, 'No data', 'Units','normalized', 'HorizontalAlignment','center','FontSize',10);
    end
    if drawMeanVector && ~isempty(ang)
        C = mean(cos(ang)); Ssin = mean(sin(ang));
        mu = atan2(Ssin, C);
        polarplot(pax, [mu mu], [0 rmax], 'LineWidth', 2);
    end
    if drawTargetOnRose && ~isnan(theta0(j))
        mu0 = deg2rad(theta0(j));
        polarplot(pax, [mu0 mu0], [0 rmax], '--', 'LineWidth', 1);
    end
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir          = 'clockwise';
    title(pax, safe_text(conditions_in{j}), 'Interpreter','none');
end
safe_sgtitle(tlo, safe_window_title('HeadingMag Rose', ws, wl));
okPDF = true;
try, exportgraphics(f, out_pdf, 'ContentType','image'); catch, okPDF = false; end
if writePNG
    try, exportgraphics(f, out_png, 'Resolution', pngRes); catch, print(f, out_png, sprintf('-r%d', pngRes), '-dpng'); end
end
if ~okPDF, print(f, out_pdf, '-dpdf', '-painters'); end
close(f);
end
