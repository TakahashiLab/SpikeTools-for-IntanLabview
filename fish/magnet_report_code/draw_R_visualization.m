
function fileMap = draw_R_visualization(S, conditions_in, outdir, ws, wl, writePNG, pngRes)
% R図の紙面拡大＋annotationで重なりゼロ化
if nargin<6, writePNG=false; end
if nargin<7, pngRes=300; end
pdfPath = fullfile(outdir, 'R_visualization.pdf');
pngPath = fullfile(outdir, 'R_visualization.png');
if exist(pdfPath,'file'), try, delete(pdfPath); end, end
if exist(pngPath,'file'), try, delete(pngPath); end, end

Rgrp = across_fish_R(S, conditions_in);

% 紙面拡大
f = figure('Color','w','Units','pixels','Position',[100 100 1200 520]);
set(f,'Units','normalized');
tlo = tiledlayout(1,2, 'Padding','compact','TileSpacing','loose');

% 左：バー＋CI
ax1 = nexttile(tlo,1); hold(ax1,'on');
x = 1:height(Rgrp);
bar(ax1, x, Rgrp.R_group);
for i=1:numel(x)
    if isfinite(Rgrp.R_CI_lo(i)) && isfinite(Rgrp.R_CI_hi(i))
        plot(ax1, [x(i) x(i)], [Rgrp.R_CI_lo(i) Rgrp.R_CI_hi(i)], '-','LineWidth',2);
    end
end
xlim(ax1, [0.5, numel(x)+0.5]); xticks(ax1,x); xticklabels(ax1,cellstr(Rgrp.Condition));
ylabel(ax1,'R (group)'); title(ax1,'Trajectory concentration (R) with 95% CI'); grid(ax1,'on');

% 右：極座標（上に持ち上げ、下にannotationでタイトル）
ax2 = nexttile(tlo,2); pos = ax2.Position; delete(ax2);
pax = polaraxes('Position', [pos(1), pos(2)+0.12*pos(4), pos(3), 0.80*pos(4)]);
hold(pax,'on');
for i=1:height(Rgrp)
    th = deg2rad(Rgrp.mu_group_deg(i)); Rv = Rgrp.R_group(i);
    polarplot(pax, [th th], [0 Rv], 'LineWidth', 3);
    text(pax, th, min(Rv+0.08, 1.02), char(Rgrp.Condition(i)), ...
        'HorizontalAlignment','center','Interpreter','none','FontSize',10);
end
pax.RLim = [0 1.05]; pax.ThetaZeroLocation='top'; pax.ThetaDir='clockwise';
% annotationのタイトル（右パネル下側）
annPos = [pos(1), pos(2)+0.02*pos(4), pos(3), 0.08*pos(4)];
annotation(f, 'textbox', annPos, 'String','Mean direction vectors (length = R)', ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',11);

try, exportgraphics(f, pdfPath, 'ContentType','image'); catch, print(f, pdfPath, '-dpdf', '-painters'); end
if writePNG
    try, exportgraphics(f, pngPath, 'Resolution', pngRes); catch, print(f, pngPath, sprintf('-r%d', pngRes), '-dpng'); end
end
close(f);

fileMap = struct('pdf', pdfPath, 'png', ternary(writePNG, pngPath, ''));
end
