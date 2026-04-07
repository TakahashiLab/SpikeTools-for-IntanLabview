function OUT = report_magnet_stats(csvPath, conditions, winStartSec, winLenSec, varargin)
% report_magnet_stats  v5.4 (circular pairwise 可視化 修正版)
%  - 条件別ローズ（平均ベクトル/目標角）
%  - 個体別ローズ（Rayleigh, *付与, タイル出力）
%  - fishMeans ローズ（サブプロット、across-fish Rayleigh、*付与）
%  - fish mean pairwise（ΔθのRayleigh）: Holm/FDR補正 + 概観図 + summary出力
%  - 連続・R・cos/sinのt検定、多重補正、GLME、Across-fish R、V-test
%  - PNGは既定ですべて出力しない（PDFのみ）

%% -------- オプション --------
p = inputParser;
addParameter(p, 'outdir', 'magnet_report', @(x)ischar(x)||isstring(x));
addParameter(p, 'alpha', 0.05, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'adjust', 'holm', @(x) any(strcmpi(x,{'holm','fdr'})));
addParameter(p, 'winlist', [], @(x) (isnumeric(x)&&isempty(x)) || (isnumeric(x)&&size(x,2)==2));

% 高速化
addParameter(p, 'doBootstrapDeltaR', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'bootstrapB', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p, 'useParallel', false, @(x)islogical(x)&&isscalar(x));

% 図/出力（PNGは既定false）
addParameter(p, 'pdfRasterRes', 300, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'drawMeanVector', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'drawTargetOnRose', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'roseBinN', 30, @(x)isnumeric(x)&&isscalar(x)&&x>=4);
addParameter(p, 'roseWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'perFishWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'fishMeansWritePNG', false, @(x)islogical(x)&&isscalar(x));

% 方向性（非一様性）
addParameter(p, 'doUniformityTest', true, @(x)islogical(x)&&isscalar(x));

% 目標方向 V-test
addParameter(p, 'doTargetTest', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'targetDirDeg', NaN, @(x)isnumeric(x));

% 個体別ローズ & fishMeans ローズ
addParameter(p, 'doPerFishPlots', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'perFishTiles', [3 4], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p, 'annotatePerFishP', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doPerFishRayleigh', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doFishMeansPlots', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'fishMeansTiles', [2 3], @(x)isnumeric(x)&&numel(x)==2);

% Paired raincloud（連続指標の可視化）
addParameter(p, 'doStatsRaincloud', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'rainTiles', [], @(x)isnumeric(x)&&(isempty(x)||(numel(x)==2)));
addParameter(p, 'rainWritePNG', false, @(x)islogical(x)&&isscalar(x));

% Δθ（fishmeanのpairwise可視化）
addParameter(p, 'doCircPairPlot', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'circPairTiles', [], @(x)isnumeric(x)&&(isempty(x)||(numel(x)==2)));
addParameter(p, 'circPairWritePNG', false, @(x)islogical(x)&&isscalar(x));

parse(p, varargin{:});
opt = p.Results;

if ischar(conditions) || isstring(conditions), conditions = cellstr(conditions); end
assert(numel(conditions)>=2, 'conditions は2つ以上を指定してください。');

% --- 複数窓モード ---
if ~isempty(opt.winlist)
    if ~isfolder(opt.outdir), mkdir(opt.outdir); end
    OUT = struct(); OUT.conditions = string(conditions);
    OUT.winlist = opt.winlist; OUT.results = cell(size(opt.winlist,1),1);
    for i = 1:size(opt.winlist,1)
        ws = opt.winlist(i,1); wl = opt.winlist(i,2);
        subOutDir = fullfile(opt.outdir, sprintf('win_%gs_%gs', ws, finite_end(ws, wl)));
        res = report_magnet_stats_core(csvPath, conditions, ws, wl, opt, subOutDir);
        OUT.results{i} = res;
    end
    return
end

% --- 単一窓モード ---
OUT = report_magnet_stats_core(csvPath, conditions, winStartSec, winLenSec, opt, string(opt.outdir));
end

%% ================= コア処理 =================
function OUT = report_magnet_stats_core(csvPath, conditions, ws, wl, opt, outdir)
if ~isfolder(outdir), mkdir(outdir); end

% -------- 読み込み --------
T = readtable(csvPath);

% 必須列チェック
needCols = {'Level','Session','Condition','Cohort','Fish','Track','Time', ...
            'Speed','HeadingMag','cosMag','sinMag','TurnRate','Curvature'};
for k = 1:numel(needCols)
    assert(ismember(needCols{k}, T.Properties.VariableNames), ...
        'CSVに必要列がありません: %s', needCols{k});
end

% sample 行のみ & 時間窓
S = T(strcmp(T.Level,'sample'), :);
if isfinite(wl)
    S = S(S.Time >= ws & S.Time <= ws + wl, :);
else
    S = S(S.Time >= ws, :);
end

% 条件名の正規化
S.Condition = lower(string(S.Condition));
conditions_in = lower(string(conditions));
conditions_in = strrep(conditions_in, "ohotukusea", "okhotsksea");
S = S(ismember(S.Condition, conditions_in), :);

% 追加指標
S = add_distance10s(S);
S = add_tortuosity(S);

% 共通Fish（対応あり）
S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
commonPairs = [];
for i = 1:numel(conditions_in)
    Pi = unique(S.Pair(S.Condition==conditions_in(i)));
    if i==1, commonPairs = Pi; else, commonPairs = intersect(commonPairs, Pi); end
end
S = S(ismember(S.Pair, commonPairs), :);

%% ---- ペアワイズ・集計 ----
METRICS = {'Speed','Distance10s','TurnRate','Curvature','Tortuosity'};
PF  = per_fish_cond_means(S, METRICS);
PFc = per_fish_cond_circular(S);
pairsCond = nchoosek(conditions_in, 2);

% 対応t（連続）
tt_cont = table();
for p = 1:size(pairsCond,1)
    cA = pairsCond(p,1); cB = pairsCond(p,2);
    A = PF(strcmp(PF.Condition,cA),:);
    B = PF(strcmp(PF.Condition,cB),:);
    [A,B] = align_pairs(A,B);
    for m = 1:numel(METRICS)
        a = A.(METRICS{m}); b = B.(METRICS{m});
        valid = ~isnan(a) & ~isnan(b);
        if any(valid)
            [~,pval,ci,stats] = ttest(a(valid), b(valid));
            dz = stats.tstat / sqrt(sum(valid));
            tt_cont = [tt_cont; table(string(cA),string(cB),string(METRICS{m}), ...
                mean(a(valid),'omitnan'), mean(b(valid),'omitnan'), ...
                mean(a(valid)-b(valid),'omitnan'), pval, stats.tstat, dz, ci(1), ci(2), ...
                'VariableNames', {'CondA','CondB','Metric','MeanA','MeanB','MeanDiff_AminusB','p','t','dz','CI_lo','CI_hi'})]; %#ok<AGROW>
        end
    end
end

% 対応t（R）
tt_R = table();
for p = 1:size(pairsCond,1)
    cA = pairsCond(p,1); cB = pairsCond(p,2);
    A = PFc(strcmp(PFc.Condition,cA),:);
    B = PFc(strcmp(PFc.Condition,cB),:);
    [A,B] = align_pairs(A,B);
    a = A.R; b = B.R;
    valid = ~isnan(a) & ~isnan(b);
    if any(valid)
        [~,pval,ci,stats] = ttest(a(valid), b(valid));
        dz = stats.tstat / sqrt(sum(valid));
        tt_R = [tt_R; table(string(cA),string(cB), ...
            mean(a(valid),'omitnan'), mean(b(valid),'omitnan'), ...
            mean(a(valid)-b(valid),'omitnan'), pval, stats.tstat, dz, ci(1), ci(2), ...
            'VariableNames', {'CondA','CondB','RA','RB','RDiff_AminusB','p','t','dz','CI_lo','CI_hi'})]; %#ok<AGROW>
    end
end

% 対応t（Fish平均 cos/sin）
tt_proj = table();
for p = 1:size(pairsCond,1)
    cA = pairsCond(p,1); cB = pairsCond(p,2);
    Z = per_fish_means_cos_sin(S, cA, cB);
    if ~isempty(Z)
        [row_cos, row_sin] = paired_t_for_means(Z, cA, cB);
        tt_proj = [tt_proj; row_cos; row_sin]; %#ok<AGROW>
    end
end

% 多重比較補正（連続・R・proj）
switch lower(string(opt.adjust))
    case "holm"
        if ~isempty(tt_cont), tt_cont.p_adj = holm_adjust(tt_cont.p); end
        if ~isempty(tt_R),    tt_R.p_adj    = holm_adjust(tt_R.p);    end
        if ~isempty(tt_proj), tt_proj.p_adj = holm_adjust(tt_proj.p); end
    case "fdr"
        if ~isempty(tt_cont), tt_cont.p_adj = mafdr(tt_cont.p,'BHFDR',true); end
        if ~isempty(tt_R),    tt_R.p_adj    = mafdr(tt_R.p,'BHFDR',true);    end
        if ~isempty(tt_proj), tt_proj.p_adj = mafdr(tt_proj.p,'BHFDR',true); end
end

% LME
S_LME = S; 
S_LME.Condition = categorical(S_LME.Condition);
S_LME.Condition = reordercats(S_LME.Condition, cellstr(conditions_in));
lme_cont = struct();
for m = 1:numel(METRICS)
    y = METRICS{m};
    mdl = fitglme_safe(S_LME, sprintf('%s ~ 1 + Condition + (1|Cohort) + (1|Fish)', y), ...
                  'Distribution','normal','Link','identity');
    if isempty(mdl), lme_cont.(y) = struct('anova', [], 'model', []);
    else,            lme_cont.(y) = struct('anova', anova(mdl), 'model', mdl);
    end
end
mdl_cos = fitglme_safe(S_LME, 'cosMag ~ 1 + Condition + (1|Cohort) + (1|Fish)');
mdl_sin = fitglme_safe(S_LME, 'sinMag ~ 1 + Condition + (1|Cohort) + (1|Fish)');
if isempty(mdl_cos), anova_cos = []; else, anova_cos = anova(mdl_cos); end
if isempty(mdl_sin), anova_sin = []; else, anova_sin = anova(mdl_sin); end

% Across-fish R
Rgrp = across_fish_R(S, conditions_in);

% ΔR ブート（任意）
if opt.doBootstrapDeltaR && opt.bootstrapB>0
    dR = delta_R_bootstrap(S, conditions_in, opt.bootstrapB, opt.useParallel);
else
    dR = table();
end

% 全体方向性（個体平均方向の Rayleigh）
uniformTbl = table();
if opt.doUniformityTest
    uniformTbl = rayleigh_uniformity_across_fish(S, conditions_in);
end

% 目標方向 V-test（任意）
targetTbl = table();
if opt.doTargetTest
    tgt = normalize_target_angles(opt.targetDirDeg, numel(conditions_in));
    targetTbl = target_direction_vtest(S, conditions_in, tgt);
end

% ★ fish mean の circular pairwise 検定（Δθ Rayleigh, 補正付）
circPair = fishmeans_circ_pairwise_tests(S, conditions_in, opt.adjust);

%% --------- 図出力 ---------
files = struct();

% 条件別・全体ローズ図（PNGは既定false）
files.rose_pdf = fullfile(outdir,'rose_headingMag.pdf');
files.rose_png = ternary(opt.roseWritePNG, fullfile(outdir,'rose_headingMag.png'), "");
draw_rose_pdf_all(S, conditions_in, files.rose_pdf, files.rose_png, ...
    ws, wl, opt.pdfRasterRes, opt.drawMeanVector, opt.drawTargetOnRose, opt.targetDirDeg, opt.roseBinN, opt.roseWritePNG);

% 個体別ローズ
perFishRayleighTbl = table();
perFishFiles = struct();
if opt.doPerFishPlots || opt.doPerFishRayleigh
    [perFishRayleighTbl, perFishFiles] = per_fish_rose_and_tests( ...
        S, conditions_in, outdir, ws, wl, opt);
end
files.perFish = perFishFiles;

% fishMeans ローズ（複数サブプロット）
fishMeansFiles = struct();
if opt.doFishMeansPlots
    fishMeansFiles = draw_fish_means_rose_multi(S, conditions_in, outdir, ws, wl, ...
        opt.pdfRasterRes, opt.drawMeanVector, uniformTbl, opt.alpha, opt.fishMeansTiles, opt.fishMeansWritePNG);
end
files.fishMeans = fishMeansFiles;

% 連続メトリクス：paired raincloud
statRainFiles = struct();
if opt.doStatsRaincloud
    statRainFiles = draw_paired_raincloud_allmetrics(PF, conditions_in, commonPairs, METRICS, outdir, ws, wl, opt.rainTiles, opt.rainWritePNG, opt.pdfRasterRes);
end
files.statPairs = statRainFiles;

% Δθ（fish mean pairwise）：可視化
circPairFiles = struct();
if opt.doCircPairPlot
    circPairFiles = plot_circ_pairwise_overview(S, conditions_in, circPair, outdir, ws, wl, opt.alpha, opt.circPairTiles, opt.circPairWritePNG, opt.pdfRasterRes);
end
files.circPairwise = circPairFiles;

%% --------- サマリ ---------
SUM = summarize_results(opt.alpha, conditions_in, ws, wl, ...
    tt_cont, tt_R, tt_proj, anova_cos, anova_sin, Rgrp, dR, targetTbl, uniformTbl, ...
    perFishRayleighTbl, circPair);
files.summary_txt = fullfile(outdir,'summary.txt');
save_lines(files.summary_txt, SUM);

%% --------- 出力 ---------
OUT = struct();
OUT.conditions           = conditions_in;
OUT.N_pairs_common       = numel(commonPairs);
OUT.window               = [ws, finite_end(ws, wl)];
OUT.pairwise.ttest_cont  = tt_cont;
OUT.pairwise.ttest_R     = tt_R;
OUT.pairwise.ttest_proj  = tt_proj;
OUT.lme.continuous       = lme_cont;
OUT.lme.cos              = struct('anova', anova_cos, 'model', mdl_cos);
OUT.lme.sin              = struct('anova', anova_sin, 'model', mdl_sin);
OUT.groupR               = Rgrp;
OUT.deltaR               = dR;
OUT.uniformityRayleigh   = uniformTbl;
OUT.targetTest           = targetTbl;
OUT.perFish.rayleigh     = perFishRayleighTbl;
OUT.circPair             = circPair;
OUT.files                = files;
OUT.summary              = SUM;

% 画面サマリ
fprintf('Common paired Fish (Cohort×Fish): %d\n', OUT.N_pairs_common);
disp('--- Paired t-tests (continuous) ---'); if ~isempty(tt_cont), disp(tt_cont); end
disp('--- Paired t-tests (R) ---');        if ~isEmptyTable(tt_R),    disp(tt_R);    end
disp('--- Fish-mean paired t-tests (cos/sin) ---'); if ~isEmptyTable(tt_proj), disp(tt_proj); end
disp('--- LME anova (cos) ---'); if ~isempty(anova_cos), disp(anova_cos); else, disp('N/A'); end
disp('--- LME anova (sin) ---'); if ~isempty(anova_sin), disp(anova_sin); else, disp('N/A'); end
disp('--- Across-fish R ---');   disp(Rgrp);
if ~isEmptyTable(uniformTbl), disp('--- Rayleigh across fish means ---'); disp(uniformTbl); end
if ~isEmptyTable(perFishRayleighTbl), disp('--- Per-fish Rayleigh ---'); disp(perFishRayleighTbl); end
if ~isEmptyTable(circPair), disp('--- Circular pairwise (fish means Δθ, adjusted) ---'); disp(circPair); end
end

%% ================= ヘルパ関数群 =================
function PF = per_fish_cond_means(S, metricNames)
S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
G = findgroups(S.Pair, S.Condition);
PF = table();
PF.Pair      = splitapply(@(x) x(1), S.Pair, G);
PF.Condition = splitapply(@(x) x(1), S.Condition, G);
for i=1:numel(metricNames)
    m = metricNames{i};
    PF.(m) = splitapply(@(x) mean(x,'omitnan'), S.(m), G);
end
end

function PFc = per_fish_cond_circular(S)
S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
G = findgroups(S.Pair, S.Condition);
Pair = splitapply(@(x) x(1), S.Pair, G);
Cond = splitapply(@(x) x(1), S.Condition, G);
mu = nan(size(Pair)); R = nan(size(Pair)); Z = nan(size(Pair)); p = nan(size(Pair));
uG = unique(G);
for gi = 1:numel(uG)
    g = uG(gi);
    idx = (G==g);
    ang = deg2rad(S.HeadingMag(idx)); ang = ang(~isnan(ang));
    if ~isempty(ang)
        C = mean(cos(ang)); Ssin = mean(sin(ang));
        mu(gi) = rad2deg(atan2(Ssin, C));
        R(gi)  = hypot(C,Ssin);
        n      = numel(ang);
        Z(gi)  = n * R(gi)^2;
        p(gi)  = exp(sqrt(1+4*n+4*(n^2 - (n*R(gi))^2)) - (1+2*n)); % 近似
        p(gi)  = min(max(p(gi),0),1);
    end
end
PFc = table(Pair, Cond, mu, R, Z, p, 'VariableNames', ...
    {'Pair','Condition','meanAngle_deg','R','RayleighZ','RayleighP'});
end

function [A,B] = align_pairs(A,B)
[~, ia, ib] = intersect(A.Pair, B.Pair, 'stable');
A = A(ia,:); B = B(ib,:);
end

function p_adj = holm_adjust(p)
[p_sorted, idx] = sort(p(:));
m = numel(p_sorted);
p_holm = zeros(m,1);
for i = 1:m
    p_holm(i) = min(1, max(p_sorted(1:i) .* (m - (0:i-1))'));
end
for i = m-1:-1:1
    p_holm(i) = max(p_holm(i), p_holm(i+1));
end
p_adj = zeros(size(p));
p_adj(idx) = p_holm;
end

function S = add_distance10s(S)
if ~ismember('Distance10s', S.Properties.VariableNames)
    S.Distance10s = nan(height(S),1);
elseif height(S.Distance10s) ~= height(S)
    S.Distance10s = nan(height(S),1);
end
keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g);
    [~, ord] = sort(S.Time(idx)); ii = idx(ord);
    if numel(ii) < 2, continue; end
    t = S.Time(ii); v = S.Speed(ii);
    dt = median(diff(t), 'omitnan');
    win = max(1, round(10 / max(dt, eps)));
    step = v .* dt;
    D = nan(numel(ii),1);
    if numel(step) >= win
        dist10 = movsum(step, [win-1, 0], 'Endpoints','discard');
        D(win:end) = dist10;
    end
    S.Distance10s(ii) = D;
end
end

function S = add_tortuosity(S)
if ~ismember('Tortuosity', S.Properties.VariableNames)
    S.Tortuosity = nan(height(S),1);
elseif height(S.Tortuosity) ~= height(S)
    S.Tortuosity = nan(height(S),1);
end
keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g);
    [~, ord] = sort(S.Time(idx)); ii = idx(ord);
    if numel(ii) < 2, continue; end
    t  = S.Time(ii);
    v  = S.Speed(ii);
    th = deg2rad(S.HeadingMag(ii));
    dt  = median(diff(t), 'omitnan');
    win = max(1, round(10 / max(dt, eps)));
    dx = v .* dt .* cos(th);
    dy = v .* dt .* sin(th);
    X  = cumsum([0; dx]); Y = cumsum([0; dy]);
    step = hypot(diff(X), diff(Y));
    TT = nan(numel(ii),1);
    if numel(step) >= win
        idxK = win:numel(step);
        path10 = movsum(step, [win-1, 0], 'Endpoints','discard');
        disp10 = nan(size(step));
        for k = idxK
            x0 = X(k-win+1); y0 = Y(k-win+1);
            x1 = X(k+1);     y1 = Y(k+1);
            disp10(k) = hypot(x1-x0, y1-y0);
        end
        tor_k = path10 ./ disp10(idxK);
        tor_k(~isfinite(tor_k)) = NaN;
        putIdx = idxK + 1; putIdx(putIdx > numel(TT)) = [];
        if ~isempty(putIdx)
            TT(putIdx) = tor_k(1:numel(putIdx));
        end
    end
    S.Tortuosity(ii) = TT;
end
end

function mdl = fitglme_safe(T, formula, varargin)
y = strtrim(formula(1:find(formula=='~',1)-1));
need = {'Condition','Cohort','Fish'};
good = true(height(T),1);
if ismember(y, T.Properties.VariableNames), good = good & ~isnan(T.(y)); end
for k=1:numel(need)
    if ismember(need{k}, T.Properties.VariableNames)
        if isnumeric(T.(need{k})); good = good & ~isnan(T.(need{k}));
        else;                       good = good & ~ismissing(T.(need{k}));
        end
    end
end
TT = T(good,:);
if isempty(TT), mdl = []; return; end
mdl = fitglme(TT, formula, varargin{:});
end

function Rgrp = across_fish_R(S, conditions_in)
outC = strings(numel(conditions_in),1);
Nf   = zeros(numel(conditions_in),1);
Rg   = nan(numel(conditions_in),1);
mu   = nan(numel(conditions_in),1);
lo   = nan(numel(conditions_in),1);
hi   = nan(numel(conditions_in),1);
for j=1:numel(conditions_in)
    Cj = conditions_in(j);
    G = findgroups(S.Pair(S.Condition==Cj));
    cosm = splitapply(@(x) mean(x,'omitnan'), S.cosMag(S.Condition==Cj), G);
    sinm = splitapply(@(x) mean(x,'omitnan'), S.sinMag(S.Condition==Cj), G);
    cosm = cosm(~isnan(cosm)); sinm = sinm(~isnan(sinm));
    n = min(numel(cosm), numel(sinm));
    cosm = cosm(1:n); sinm=sinm(1:n);
    Cbar = mean(cosm); Sbar = mean(sinm);
    Rg(j) = hypot(Cbar,Sbar);
    mu(j) = rad2deg(atan2(Sbar,Cbar));
    Nf(j) = n;
    if n>1
        B=1000; Rs=zeros(B,1);
        for b=1:B
            idx = randsample(n,n,true);
            Rs(b) = hypot(mean(cosm(idx)), mean(sinm(idx)));
        end
        lo(j) = quantile(Rs,0.025); hi(j)=quantile(Rs,0.975);
    end
    outC(j) = Cj;
end
Rgrp = table(outC, Nf, Rg, mu, lo, hi, ...
    'VariableNames',{'Condition','N_fish','R_group','mu_group_deg','R_CI_lo','R_CI_hi'});
end

function dR = delta_R_bootstrap(S, conditions_in, B, useParallel)
pairs = nchoosek(conditions_in,2);
dR = table();
for p=1:size(pairs,1)
    A  = pairs(p,1); Bc = pairs(p,2);
    Ra = one_R_global(S, A);
    Rb = one_R_global(S, Bc);
    if isnan(Ra) || isnan(Rb)
        dR = [dR; table(string(A),string(Bc),NaN,NaN,NaN,NaN, ...
            'VariableNames',{'CondA','CondB','Delta','CI_lo','CI_hi','p_boot'})];
        continue
    end
    fish = unique(S.Pair(ismember(S.Condition,[A;Bc])));
    n = numel(fish);
    D = zeros(B,1);
    if useParallel && license('test','Distrib_Computing_Toolbox')
        parfor b=1:B
            pick = fish(randi(n,n,1));
            Ra_b = one_R_from_pairs(S, A, pick);
            Rb_b = one_R_from_pairs(S, Bc, pick);
            D(b)  = Ra_b - Rb_b;
        end
    else
        for b=1:B
            pick = fish(randi(n,n,1));
            Ra_b = one_R_from_pairs(S, A, pick);
            Rb_b = one_R_from_pairs(S, Bc, pick);
            D(b)  = Ra_b - Rb_b;
        end
    end
    dR = [dR; table(string(A),string(Bc), Ra-Rb, quantile(D,0.025), quantile(D,0.975), ...
        mean((D - (Ra-Rb))>=0), 'VariableNames',{'CondA','CondB','Delta','CI_lo','CI_hi','p_boot'})]; %#ok<AGROW>
end
end

function Rval = one_R_global(S, Cname)
Ssub = S(S.Condition==Cname,:);
Rval = one_R_core(Ssub);
end

function Rval = one_R_from_pairs(S, Cname, allowedPairs)
mask = (S.Condition==Cname) & ismember(S.Pair, allowedPairs);
Ssub = S(mask,:);
Rval = one_R_core(Ssub);
end

function Rval = one_R_core(Ssub)
if isempty(Ssub), Rval = NaN; return; end
G = findgroups(Ssub.Pair);
cosm = splitapply(@(x) mean(x,'omitnan'), Ssub.cosMag, G);
sinm = splitapply(@(x) mean(x,'omitnan'), Ssub.sinMag, G);
cosm = cosm(~isnan(cosm)); sinm = sinm(~isnan(sinm));
if isempty(cosm) || isempty(sinm), Rval = NaN; return; end
Rval = hypot(mean(cosm), mean(sinm));
end

function Z = per_fish_means_cos_sin(S, cA, cB)
Pa = unique(S.Pair(S.Condition==cA));
Pb = unique(S.Pair(S.Condition==cB));
P  = intersect(Pa,Pb,'stable');
if isempty(P), Z = table(); return; end
Z = table('Size',[0 5], ...
          'VariableTypes', {'string','double','double','double','double'}, ...
          'VariableNames', {'Pair','cosA','cosB','sinA','sinB'});
for k=1:numel(P)
    p = P(k);
    Sa = S(S.Pair==p & S.Condition==cA,:);
    Sb = S(S.Pair==p & S.Condition==cB,:);
    row = table( string(p), ...
                 mean(Sa.cosMag,'omitnan'), mean(Sb.cosMag,'omitnan'), ...
                 mean(Sa.sinMag,'omitnan'), mean(Sb.sinMag,'omitnan'), ...
                 'VariableNames', {'Pair','cosA','cosB','sinA','sinB'});
    Z = [Z; row]; %#ok<AGROW>
end
end

function [row_cos, row_sin] = paired_t_for_means(Z, cA, cB)
[~,p,ci,st] = ttest(Z.cosA, Z.cosB);
dz = st.tstat / sqrt(numel(Z.cosA));
row_cos = table(string(cA), string(cB), "cos", mean(Z.cosA), mean(Z.cosB), ...
          mean(Z.cosA-Z.cosB), p, st.tstat, dz, ci(1), ci(2), height(Z), ...
          'VariableNames', {'CondA','CondB','Var','MeanA','MeanB','MeanDiff','p','t','dz','CI_lo','CI_hi','N'});
[~,p,ci,st] = ttest(Z.sinA, Z.sinB);
dz = st.tstat / sqrt(numel(Z.sinA));
row_sin = table(string(cA), string(cB), "sin", mean(Z.sinA), mean(Z.sinB), ...
          mean(Z.sinA-Z.sinB), p, st.tstat, dz, ci(1), ci(2), height(Z), ...
          'VariableNames', {'CondA','CondB','Var','MeanA','MeanB','MeanDiff','p','t','dz','CI_lo','CI_hi','N'});
end

%% ===== 図：条件ごとの全体ローズ =====
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
try, exportgraphics(f, out_pdf, 'ContentType','vector'); catch, okPDF = false; end
if writePNG
    try, exportgraphics(f, out_png, 'Resolution', pngRes); catch, print(f, out_png, sprintf('-r%d', pngRes), '-dpng'); end
end
if ~okPDF, print(f, out_pdf, '-dpdf', '-painters'); end
close(f);
end

%% ===== 個体別ローズ + Rayleigh =====
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
            rayC(end+1,1)    = cj;                 %#ok<AGROW>
            rayPair(end+1,1) = string(p);          %#ok<AGROW>
            rayN(end+1,1)    = n;                  %#ok<AGROW>
            rayMean(end+1,1) = NaN;                %#ok<AGROW>
            rayR(end+1,1)    = Rbar;               %#ok<AGROW>
            rayZ(end+1,1)    = z;                  %#ok<AGROW>
            rayP(end+1,1)    = pval;               %#ok<AGROW>
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

function finalize_page(figH, pdfPath, pngPath, pngRes, appendPDF, writePNG)
try
    exportgraphics(figH, pdfPath, 'ContentType','vector', 'Append', appendPDF);
catch
    print(figH, pdfPath, '-dpdf', '-painters');
end
if writePNG
    try
        exportgraphics(figH, pngPath, 'Resolution', pngRes, 'Append', appendPDF);
    catch
        print(figH, pngPath, sprintf('-r%d', pngRes), '-dpng');
    end
end
end

function [pval, z, Rbar, mu_deg, n] = rayleigh_on_angles(alpha_rad)
alpha_rad = alpha_rad(:); alpha_rad = alpha_rad(~isnan(alpha_rad));
n = numel(alpha_rad);
if n==0, pval=NaN; z=NaN; Rbar=NaN; mu_deg=NaN; return; end
Cbar = mean(cos(alpha_rad)); Sbar = mean(sin(alpha_rad));
Rbar = hypot(Cbar, Sbar);
mu_deg = rad2deg(atan2(Sbar, Cbar));
if exist('circ_rtest','file') == 2
    [pval, z] = circ_rtest(alpha_rad);
else
    z = n * Rbar^2;
    pval = exp(sqrt(1+4*n+4*(n^2 - (n*Rbar)^2)) - (1+2*n));
    pval = min(max(pval,0),1);
end
end

function annotation_on_axes(ax, txt)
text(ax, 0.02, 0.02, txt, 'Units','normalized', ...
    'VerticalAlignment','bottom','HorizontalAlignment','left', ...
    'FontSize',8, 'Interpreter','none');
end

function add_sig_star(ax, nStars)
if nargin<2, nStars = 1; end
stars = repmat('*',1,nStars);
text(ax, 0.95, 0.95, stars, 'Units','normalized', 'HorizontalAlignment','right', ...
    'VerticalAlignment','top', 'FontWeight','bold', 'FontSize', 14);
end

function s = safe_field(x)
s = char(safe_text(x)); s(~isstrprop(s,'alphanum')) = '_';
end

%% ===== 個体平均方向のローズ図（複数サブプロットを1つのPDF/PNGへ） =====
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

    % 個体平均角（fish means）
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

    % Rayleigh（across fish means）のp
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

%% ===== 統計用：paired raincloud =====
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
end

try, exportgraphics(f, pdfPath, 'ContentType','vector'); catch, print(f, pdfPath, '-dpdf', '-painters'); end
if writePNG
    try, exportgraphics(f, pngPath, 'Resolution', pngRes); catch, print(f, pngPath, sprintf('-r%d', pngRes), '-dpng'); end
end
close(f);

fileMap.pdf = pdfPath; fileMap.png = ternary(writePNG, pngPath, '');
end

%% ===== Δθ 可視化：ペアごとのレインクラウド + 生点 + メディアン + コンパス =====
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

    % コンパス（インセット）: polaraxes を直接作る（<< 修正点）
    pax = polaraxes('Position', inset_rect(ax, 0.70, 0.60, 0.25, 0.35)); hold(pax,'on');
    C = mean(cos(d)); Ssin = mean(sin(d));
    mu = atan2(Ssin, C); Rbar = hypot(C,Ssin);
    polarplot(pax, [mu mu], [0 Rbar], 'LineWidth', 2);
    pax.RLim = [0 1]; pax.ThetaZeroLocation='top'; pax.ThetaDir='clockwise';
    title(pax, 'mean Δ', 'FontSize',9);
end

try, exportgraphics(f, pdfPath, 'ContentType','vector'); catch, print(f, pdfPath, '-dpdf', '-painters'); end
if writePNG
    try, exportgraphics(f, pngPath, 'Resolution', pngRes); catch, print(f, pngPath, sprintf('-r%d', pngRes), '-dpng'); end
end
close(f);

fileMap = struct('pdf', pdfPath, 'png', ternary(writePNG, pngPath, ''));
end

function d = compute_deltas_rad(S, condCol, cA, cB)
SjA = S(condCol==cA,:); SjB = S(condCol==cB,:);
[Pa, muA] = fish_mean_angles(SjA);
[Pb, muB] = fish_mean_angles(SjB);
[~, ia, ib] = intersect(Pa, Pb, 'stable');
if isempty(ia), d = []; return; end
d = local_wrapToPi(muA(ia) - muB(ib));
end

function [P, mu] = fish_mean_angles(Ssub)
if ~ismember('Pair', Ssub.Properties.VariableNames)
    Ssub.Pair = strcat("C", string(Ssub.Cohort), "_F", string(Ssub.Fish));
end
G = findgroups(Ssub.Pair);
cosm = splitapply(@(x) mean(x,'omitnan'), Ssub.cosMag, G);
sinm = splitapply(@(x) mean(x,'omitnan'), Ssub.sinMag, G);
ok = ~isnan(cosm) & ~isnan(sinm);
mu = atan2(sinm(ok), cosm(ok));
P  = splitapply(@(x) x(1), Ssub.Pair, G); P = P(ok);
end

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

%% ===== 全体の Rayleigh（個体平均方向） =====
function uniformTbl = rayleigh_uniformity_across_fish(S, conditions_in)
Cname = strings(0,1); Nfish = []; meanAng = []; Rbar = []; Zs = []; Ps = [];
for j = 1:numel(conditions_in)
    cj = conditions_in(j);
    Sj = S(S.Condition==cj, :);
    if isempty(Sj), continue; end
    G    = findgroups(Sj.Pair);
    cosm = splitapply(@(x) mean(x,'omitnan'), Sj.cosMag, G);
    sinm = splitapply(@(x) mean(x,'omitnan'), Sj.sinMag, G);
    ok   = ~isnan(cosm) & ~isnan(sinm);
    cosm = cosm(ok); sinm = sinm(ok);
    if isempty(cosm), continue; end
    alpha = atan2(sinm, cosm);
    n = numel(alpha);
    Cbar = mean(cos(alpha)); Sbar = mean(sin(alpha));
    mu   = atan2(Sbar,Cbar);
    Rb   = hypot(Cbar,Sbar);
    if exist('circ_rtest','file') == 2
        [p, z] = circ_rtest(alpha);
    else
        z = n * Rb^2;
        p = exp(sqrt(1+4*n+4*(n^2 - (n*Rb)^2)) - (1+2*n));
        p = min(max(p,0),1);
    end
    Cname(end+1,1)   = cj;               %#ok<AGROW>
    Nfish(end+1,1)   = n;                %#ok<AGROW>
    meanAng(end+1,1) = rad2deg(mu);      %#ok<AGROW>
    Rbar(end+1,1)    = Rb;               %#ok<AGROW>
    Zs(end+1,1)      = z;                %#ok<AGROW>
    Ps(end+1,1)      = p;                %#ok<AGROW>
end
uniformTbl = table(Cname, Nfish, meanAng, Rbar, Zs, Ps, ...
    'VariableNames', {'Condition','N_fish','mean_angle_deg','Rbar_fish','Z','p'});
end

%% ===== V-test（目標角） =====
function T = target_direction_vtest(S, conditions_in, targetAngles)
Cname = strings(0,1); Nfish = []; theta0= []; meanAngDeg = []; Rbar = []; meanProj = []; Vstat= []; p_one = [];
for j = 1:numel(conditions_in)
    cj  = conditions_in(j);
    th0 = targetAngles(j);
    if isnan(th0), continue; end

    Sj = S(S.Condition==cj, :); if isempty(Sj), continue; end
    G    = findgroups(Sj.Pair);
    cosm = splitapply(@(x) mean(x,'omitnan'), Sj.cosMag, G);
    sinm = splitapply(@(x) mean(x,'omitnan'), Sj.sinMag, G);
    ok   = ~isnan(cosm) & ~isnan(sinm);
    cosm = cosm(ok); sinm = sinm(ok);
    if isempty(cosm), continue; end

    alpha = atan2(sinm, cosm); n = numel(alpha); mu0 = deg2rad(th0);
    Cbar = mean(cos(alpha)); Sbar = mean(sin(alpha));
    meanAngDeg(end+1,1) = rad2deg(atan2(Sbar,Cbar)); %#ok<AGROW>
    Rbar(end+1,1)       = hypot(Cbar,Sbar);          %#ok<AGROW>

    if exist('circ_vtest','file') == 2
        [p, V] = circ_vtest(alpha, mu0); Vstat(end+1,1)=V; p_one(end+1,1)=p; meanProj(end+1,1)=mean(cos(alpha-mu0));
    else
        u = mean(cos(alpha - mu0)); V = sqrt(2*n) * u; p = 1 - normcdf(V);
        Vstat(end+1,1)=V; p_one(end+1,1)=p; meanProj(end+1,1)=u;
    end

    Cname(end+1,1)=cj; Nfish(end+1,1)=n; theta0(end+1,1)=th0;
end
T = table(Cname, Nfish, theta0, meanAngDeg, Rbar, meanProj, Vstat, p_one, ...
    'VariableNames', {'Condition','N_fish','target_deg','mean_angle_deg','Rbar','meanProjection','Vstat','p_one_sided'});
end

%% ===== fish mean circular pairwise（Δθ Rayleigh, 補正付） =====
function circPair = fishmeans_circ_pairwise_tests(S, conditions_in, adjustMethod)
conds = string(conditions_in(:));
if ~ismember('Pair', S.Properties.VariableNames)
    S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
end
condCol = string(S.Condition);

A_all = strings(0,1); B_all = strings(0,1);
N_all = []; Z_all = []; P_all = []; MeanDeltaDeg_all = []; Rbar_all = [];

pairs = nchoosek(1:numel(conds),2);
for k = 1:size(pairs,1)
    cA = conds(pairs(k,1)); cB = conds(pairs(k,2));
    d = compute_deltas_rad(S, condCol, cA, cB);
    d = d(~isnan(d));
    if isempty(d)
        A_all(end+1,1)=cA; B_all(end+1,1)=cB;
        N_all(end+1,1)=0; Z_all(end+1,1)=NaN; P_all(end+1,1)=NaN;
        MeanDeltaDeg_all(end+1,1)=NaN; Rbar_all(end+1,1)=NaN; %#ok<AGROW>
        continue;
    end
    n = numel(d);
    C = mean(cos(d)); Ssin = mean(sin(d));
    Rbar = hypot(C,Ssin); Z = n * Rbar^2;
    if exist('circ_rtest','file')==2
        [pval, ~] = circ_rtest(d);
    else
        pval = exp(sqrt(1+4*n+4*(n^2 - (n*Rbar)^2)) - (1+2*n));
        pval = min(max(pval,0),1);
    end
    mu = atan2(Ssin, C);
    A_all(end+1,1)=cA; B_all(end+1,1)=cB;
    N_all(end+1,1)=n; Z_all(end+1,1)=Z; P_all(end+1,1)=pval;
    MeanDeltaDeg_all(end+1,1)=rad2deg(mu); Rbar_all(end+1,1)=Rbar; %#ok<AGROW>
end

circPair = table(A_all, B_all, N_all, Z_all, P_all, MeanDeltaDeg_all, Rbar_all, ...
    'VariableNames', {'CondA','CondB','N','Z','p','MeanDelta_deg','Rbar'});

% 調整
switch lower(string(adjustMethod))
    case "holm"
        circPair.p_adj = holm_adjust(circPair.p);
    case "fdr"
        circPair.p_adj = mafdr(circPair.p,'BHFDR',true);
    otherwise
        circPair.p_adj = circPair.p;
end
end

%% ===== サマリ =====
function txt = summarize_results(alpha, conds, ws, wl, tt_cont, tt_R, tt_proj, anova_cos, anova_sin, Rgrp, dR, targetTbl, uniformTbl, perFishRayleighTbl, circPair)
lines = {};
lines{end+1} = sprintf('Window %gs-%gs | alpha=%.3f', ws, finite_end(ws, wl), alpha);

% 連続
if ~isEmptyTable(tt_cont)
    sig = pick_sig(tt_cont, alpha);
    if ~isempty(sig)
        lines{end+1} = 'Significant (continuous):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s vs %s: %s  Δ=%.3g (p=%.3g%s)', ...
                sig.CondA(i), sig.CondB(i), sig.Metric(i), sig.MeanDiff_AminusB(i), ...
                pick_p(sig,i), pick_padj(sig,i));
        end
    end
end

% R
if ~isEmptyTable(tt_R)
    sig = pick_sig(tt_R, alpha);
    if ~isempty(sig)
        lines{end+1} = 'Significant (R):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s vs %s: ΔR=%.3g (p=%.3g%s)', ...
                sig.CondA(i), sig.CondB(i), sig.RDiff_AminusB(i), ...
                pick_p(sig,i), pick_padj(sig,i));
        end
    end
end

% Fish-mean cos/sin
if ~isEmptyTable(tt_proj)
    sig = pick_sig(tt_proj, alpha);
    if ~isempty(sig)
        lines{end+1} = 'Significant (Fish-mean cos/sin):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s vs %s: %s Δ=%.3g (p=%.3g%s)', ...
                sig.CondA(i), sig.CondB(i), sig.Var(i), sig.MeanDiff(i), ...
                pick_p(sig,i), pick_padj(sig,i));
        end
    end
end

% LME 概要
pc = pick_p_from_anova(anova_cos,'Condition');
ps = pick_p_from_anova(anova_sin,'Condition');
lines{end+1} = sprintf('LME ProjectedNormal: cos p=%.3g | sin p=%.3g', pc, ps);

% Across-fish R
if ~isempty(Rgrp)
    lines{end+1} = 'Across-fish R:';
    for i=1:height(Rgrp)
        lines{end+1} = sprintf('  %s: R=%.3g (95%%CI %.3g-%.3g) mu=%.1f° (N=%d)', ...
            Rgrp.Condition(i), Rgrp.R_group(i), Rgrp.R_CI_lo(i), Rgrp.R_CI_hi(i), ...
            Rgrp.mu_group_deg(i), Rgrp.N_fish(i));
    end
end

% ΔR
if ~isEmptyTable(dR)
    sig = dR(dR.CI_lo>0 | dR.CI_hi<0, :);
    if ~isempty(sig)
        lines{end+1} = 'Δ Across-fish R (bootstrap, 95%CI not crossing 0):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s - %s: ΔR=%.3g (CI %.3g-%.3g), p_boot=%.3g', ...
                sig.CondA(i), sig.CondB(i), sig.Delta(i), sig.CI_lo(i), sig.CI_hi(i), sig.p_boot(i));
        end
    end
end

% 全体方向性（Rayleigh across fish means）
if ~isEmptyTable(uniformTbl)
    sig = uniformTbl(uniformTbl.p < alpha, :);
    if ~isempty(sig)
        lines{end+1} = 'Directional bias (Rayleigh across fish means):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s: mean=%.1f°, R̄_fish=%.3f, Z=%.2f, p=%.3g (N=%d)', ...
                sig.Condition(i), sig.mean_angle_deg(i), sig.Rbar_fish(i), sig.Z(i), sig.p(i), sig.N_fish(i));
        end
    end
end

% 目標方向 V-test
if ~isEmptyTable(targetTbl)
    sig = targetTbl(targetTbl.p_one_sided < alpha, :);
    if ~isempty(sig)
        lines{end+1} = 'Target direction (circular V-test, one-sided):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s: θ0=%g°, meanProj=%.3f, V=%.2f, p=%.3g, meanAngle=%.1f°, R̄=%.3f (N=%d)', ...
                sig.Condition(i), sig.target_deg(i), sig.meanProjection(i), ...
                sig.Vstat(i), sig.p_one_sided(i), sig.mean_angle_deg(i), sig.Rbar(i), sig.N_fish(i));
        end
    end
end

% 個体別 Rayleigh（件数サマリ）
if ~isEmptyTable(perFishRayleighTbl)
    lines{end+1} = 'Per-fish Rayleigh (counts p<alpha):';
    for j=1:numel(conds)
        cj = conds(j);
        rows = perFishRayleighTbl(perFishRayleighTbl.Condition==cj,:);
        if ~isempty(rows)
            nSig = sum(rows.p < alpha);
            lines{end+1} = sprintf('  %s: %d/%d fish significant', cj, nSig, height(rows));
        end
    end
end

% ★ fish mean pairwise（Δθ Rayleigh, 補正付）を summary に追記
if ~isEmptyTable(circPair)
    sig = circPair(circPair.p_adj < alpha, :);
    if ~isempty(sig)
        lines{end+1} = 'Circular pairwise (fish means Δθ, adjusted):';
        for i=1:height(sig)
            lines{end+1} = sprintf('  %s vs %s: μΔ=%.1f°, R̄=%.3f, Z=%.2f, p=%.3g, p_adj=%.3g (N=%d)', ...
                sig.CondA(i), sig.CondB(i), sig.MeanDelta_deg(i), sig.Rbar(i), sig.Z(i), sig.p(i), sig.p_adj(i), sig.N(i));
        end
    end
end

txt = string(lines(:));

    function T2 = pick_sig(Tin, a)
        if any(ismember(Tin.Properties.VariableNames,'p_adj'))
            T2 = Tin(Tin.p_adj < a, :);
        elseif any(ismember(Tin.Properties.VariableNames,'p'))
            T2 = Tin(Tin.p < a, :);
        else
            T2 = Tin([],:);
        end
    end
    function s = pick_padj(T,i)
        if any(ismember(T.Properties.VariableNames,'p_adj'))
            s = sprintf(', p_adj=%.3g', T.p_adj(i));
        else
            s = '';
        end
    end
    function pv = pick_p(T,i)
        if any(ismember(T.Properties.VariableNames,'p'))
            pv = T.p(i);
        else
            pv = NaN;
        end
    end
end

%% ===== 小物 =====
function pv = pick_p_from_anova(anovTbl, term)
pv = NaN;
if isempty(anovTbl), return; end
try, vnames = anovTbl.Properties.VariableNames; catch, pv=NaN; return; end
hasP = any(strcmpi(vnames,'pValue'));
hasT = any(strcmpi(vnames,'Term'));
if ~(hasP && hasT), return; end
try
    idx = strcmpi(anovTbl.Term, term);
    if any(idx), pv = anovTbl.pValue(find(idx,1)); end
catch
    pv = NaN;
end
end

function tf = isEmptyTable(T)
tf = isempty(T) || (istable(T) && height(T)==0);
end

function txt = safe_window_title(prefix, ws, wl)
txt = safe_text(sprintf('%s (Window %gs-%gs)', prefix, ws, finite_end(ws, wl)));
end

function x = finite_end(ws, wl)
if isfinite(wl), x = ws+wl; else, x = Inf; end
end

function safe_sgtitle(parentOrFig, titleStr)
titleStr = safe_text(titleStr);
try
    sgtitle(parentOrFig, titleStr, 'Interpreter','none');
catch
    try
        sgtitle(titleStr, 'Interpreter','none');
    catch
        if isempty(get(0,'CurrentFigure')), figure('Color','w'); end
        annotation(gcf,'textbox',[0 0.95 1 0.05], ...
            'String', titleStr, 'EdgeColor','none', ...
            'HorizontalAlignment','center', 'Interpreter','none');
    end
end
end

function txt = safe_text(x)
if isstring(x), x = x(1); end
if iscell(x),   x = x{1}; end
if isnumeric(x), x = num2str(x); end
txt = char(x);
txt = strrep(txt, char(8211), '-');
txt = strrep(txt, char(8212), '-');
end

function save_lines(path, lines_str)
try
    writelines(lines_str, path);
catch
    fid = fopen(path, 'w');
    if fid<0, warning('Cannot open %s', path); return; end
    for i=1:numel(lines_str)
        fprintf(fid, '%s\n', safe_text(lines_str(i)));
    end
    fclose(fid);
end
end

% --- ローカル：wrapToPi 代替（Mapping Toolbox 非依存） ---
function ang = local_wrapToPi(ang)
ang = mod(ang + pi, 2*pi) - pi;
end


function theta = normalize_target_angles(targetDirDeg, N)
if isscalar(targetDirDeg)
    theta = repmat(targetDirDeg, 1, N);
else
    theta = targetDirDeg(:).';
    if numel(theta) < N
        theta = [theta, nan(1, N-numel(theta))];
    elseif numel(theta) > N
        theta = theta(1:N);
    end
end
end

function r = inset_rect(ax, x, y, w, h)
op = get(ax,'Position'); r = [op(1)+op(3)*x, op(2)+op(4)*y, op(3)*w, op(4)*h];
end

function out = ternary(tf,a,b)
if tf, out=a; else, out=b; end
end
