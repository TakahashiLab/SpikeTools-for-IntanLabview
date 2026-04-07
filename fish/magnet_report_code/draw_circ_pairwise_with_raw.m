function fileMap = draw_circ_pairwise_with_raw(S, conditions_in, outdir, ws, wl, alpha, writePNG, pngRes)
% Pairwise display: [Cond A fish-mean] | [Delta(A-B)] | [Cond B fish-mean]
% - タイルは明示 index 指定 (nexttile(tlo, r)) し、削除しない（重なり防止）
% - ラベル/統計は annotation('textbox', ...) で軸外に配置
% - 左右パネルにも赤の平均ベクトルを追加

if nargin<6 || isempty(alpha), alpha=0.05; end
if nargin<7, writePNG=false; end
if nargin<8, pngRes=300; end

pdfPath = fullfile(outdir, 'circ_pairwise_with_raw.pdf');
pngPath = fullfile(outdir, 'circ_pairwise_with_raw.png');
if exist(pdfPath,'file'), try, delete(pdfPath); end, end
if exist(pngPath,'file'), try, delete(pngPath); end, end

conds = string(conditions_in(:));
pairs = nchoosek(1:numel(conds),2);
nPairs = size(pairs,1);

% 図サイズを行数に応じて拡大
pxPerRow  = 440;                   % 行あたり高さ（必要に応じて 480–520 に）
figWidth  = 1280;
figHeight = max(400, pxPerRow*nPairs + 80);
f = figure('Color','w','Units','pixels','Position',[100 80 figWidth figHeight]);
set(f,'Units','normalized');
tlo = tiledlayout(nPairs, 1, 'Padding','compact','TileSpacing','loose');

for r = 1:nPairs
    cA = conds(pairs(r,1)); 
    cB = conds(pairs(r,2));

    % Δ分布と各条件の魚平均角
    d = compute_deltas_rad(S, string(S.Condition), cA, cB);
    d = d(~isnan(d));
    SA = S(S.Condition==cA, :); 
    SB = S(S.Condition==cB, :);
    [~, muA] = fish_mean_angles(SA);   % 各魚の平均角（A）
    [~, muB] = fish_mean_angles(SB);   % 各魚の平均角（B）

    % 行タイル（削除しない -> 重なり防止）
    axpanel = nexttile(tlo, r); 
    set(axpanel, 'Visible','off');  % 位置だけ使う
    pos = axpanel.Position;

    % タイル内に 3 極座標を手配置（下にラベル帯を確保）
    margin = 0.06;
    w = pos(3); h = pos(4);
    w3 = (w - 4*margin)/3;
    y0 = pos(2) + 0.12*h;        % 下側にラベル/統計帯
    hh = h * 0.76;

    leftPos   = [pos(1)+margin,          y0, w3, hh];
    centerPos = [pos(1)+2*margin + w3,   y0, w3, hh];
    rightPos  = [pos(1)+3*margin + 2*w3, y0, w3, hh];

    red = [0.85 0.20 0.20];

    % ===== 左: 条件A =====
    paxA = polaraxes('Position', leftPos); hold(paxA,'on');
    polarhistogram(paxA, muA, 18, 'Normalization','probability');
    paxA.ThetaZeroLocation='top'; 
    paxA.ThetaDir='clockwise';

    % A の平均ベクトル（赤）
    if ~isempty(muA)
        CA = mean(cos(muA)); SA_ = mean(sin(muA));
        muAbar = atan2(SA_, CA); RA = hypot(CA, SA_);
        polarplot(paxA, [muAbar muAbar], [0 RA], 'Color', red, 'LineWidth', 2);
    end

    % A ラベル（軸外下段）
    aPos = [leftPos(1), pos(2)+0.02*h, leftPos(3), 0.07*h];
    annotation(f, 'textbox', aPos, 'String', char(cA), 'EdgeColor','none', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle', 'Interpreter','none');

    % ===== 中央: Δ分布 =====
    paxC = polaraxes('Position', centerPos); hold(paxC,'on');
    if isempty(d), d = 0; end
    polarhistogram(paxC, d, 18, 'Normalization','probability');

    Cbar = mean(cos(d)); Sbar = mean(sin(d));
    mu  = atan2(Sbar, Cbar); 
    Rbar = hypot(Cbar,Sbar);
    polarplot(paxC, [mu mu], [0 Rbar], 'Color', red, 'LineWidth', 2);  % Δ の平均ベクトル

    paxC.RLim = [0 1.05]; 
    paxC.ThetaZeroLocation='top'; 
    paxC.ThetaDir='clockwise';

    % 中央見出し＋統計（軸外下段 2 行）
    cHeadPos = [centerPos(1), pos(2)+0.045*h, centerPos(3), 0.07*h];
    annotation(f, 'textbox', cHeadPos, 'String', sprintf('%s vs %s', cA, cB), ...
        'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none');

    pv = NaN; 
    try, [~,pv] = rayleigh_on_angles(d); catch, end
    statStr = sprintf('\\Delta\\mu=%.1f^\\circ   R=%.2f   p=%.3g', rad2deg(mu), Rbar, pv);
    cStatPos = [centerPos(1), pos(2)+0.005*h, centerPos(3), 0.07*h];
    annotation(f, 'textbox', cStatPos, 'String', statStr, ...
        'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle');

    % ===== 右: 条件B =====
    paxB = polaraxes('Position', rightPos); hold(paxB,'on');
    polarhistogram(paxB, muB, 18, 'Normalization','probability');
    paxB.ThetaZeroLocation='top'; 
    paxB.ThetaDir='clockwise';

    % B の平均ベクトル（赤）
    if ~isempty(muB)
        CB = mean(cos(muB)); SB_ = mean(sin(muB));
        muBbar = atan2(SB_, CB); RB = hypot(CB, SB_);
        polarplot(paxB, [muBbar muBbar], [0 RB], 'Color', red, 'LineWidth', 2);
    end

    % B ラベル（軸外下段）
    bPos = [rightPos(1), pos(2)+0.02*h, rightPos(3), 0.07*h];
    annotation(f, 'textbox', bPos, 'String', char(cB), 'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none');
end

% 出力
try, exportgraphics(f, pdfPath, 'ContentType','image'); catch, print(f, pdfPath, '-dpdf', '-painters'); end
if writePNG
    try, exportgraphics(f, pngPath, 'Resolution', pngRes); catch, print(f, pngPath, sprintf('-r%d', pngRes), '-dpng'); end
end
close(f);

fileMap = struct('pdf', pdfPath, 'png', ternary(writePNG, pngPath, ''));
end
