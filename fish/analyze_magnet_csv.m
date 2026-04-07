function OUT = analyze_magnet_csv(csvPath, conditions, winStartSec, winLenSec)
% analyze_magnet_csv_multi
% 複数条件（>=2）を、同一 Fish（Cohort×Fish）で“対応あり”比較する一括解析。
% - 連続指標: Speed, Distance10s, TurnRate, Curvature, Tortuosity
%   * Fish対応の全ペア t 検定
%   * LME: metric ~ Condition + (1|Cohort) + (1|Fish)
% - 円統計（Mag North=0°）：meanAngle_deg, R
%   * R の全ペア t 検定
%   * Projected Normal LME: cosMag/sinMag ~ Condition + (1|Cohort) + (1|Fish)
%
% 例:
% OUT = analyze_magnet_csv_multi('speed_heading_data2025.csv',...
%     {'beringSea','okhotskSea','controlSea1'}, 9*60, inf);

if nargin < 4, winLenSec = inf; end
if ischar(conditions) || isstring(conditions), conditions = cellstr(conditions); end
assert(numel(conditions)>=2, 'conditions は2つ以上指定してください。');

% ---- 読み込み ----
T = readtable(csvPath);

% 必須列チェック
needCols = {'Level','Session','Condition','Cohort','Fish','Track','Time', ...
            'Speed','HeadingMag','cosMag','sinMag','TurnRate','Curvature'};
for k = 1:numel(needCols)
    if ~ismember(needCols{k}, T.Properties.VariableNames)
        error('CSVに必要列がありません: %s', needCols{k});
    end
end

% Level=sample のみ使用
S = T(strcmp(T.Level,'sample'), :);

% 時間窓フィルタ
if isfinite(winLenSec)
    S = S(S.Time >= winStartSec & S.Time <= winStartSec + winLenSec, :);
else
    S = S(S.Time >= winStartSec, :);
end

% 条件名の正規化（小文字照合 & 「OhotukuSea」→「okhotskSea」）
S.Condition = lower(string(S.Condition));
conditions_in = lower(string(conditions));
% conditions_in = strrep(conditions_in, "ohotukusea", "okhotsksea");

% 存在確認
conds_exist = intersect(unique(S.Condition), conditions_in);
if numel(conds_exist) < numel(conditions_in)
    missing = setdiff(conditions_in, conds_exist);
    warning('CSVに無い条件が含まれています: %s', strjoin(missing, ', '));
end
conditions_use = conditions_in(:)';

S = S(ismember(S.Condition, conditions_use), :);

% 10秒距離とトータス（10秒窓）を付与
S = add_distance10s(S);
S = add_tortuosity(S);

% ---- 条件ごとの Fish（Cohort×Fish）集合 -> 全条件共通の Fish のみ残す ----
S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
commonPairs = [];
for i = 1:numel(conditions_use)
    Pi = unique(S.Pair(S.Condition==conditions_use(i)));
    if i==1, commonPairs = Pi; else, commonPairs = intersect(commonPairs, Pi); end
end
S = S(ismember(S.Pair, commonPairs), :);

% ---- 連続指標の Fish×Condition 平均（窓内） ----
METRICS = {'Speed','Distance10s','TurnRate','Curvature','Tortuosity'};
PF = per_fish_cond_means(S, METRICS);   % columns: Pair, Condition, (metrics...)

% ---- 円統計（Fish×Condition）：meanAngle_deg, R, Rayleigh ----
PFc = per_fish_cond_circular(S);

% ---- 全ペア（条件間）の対応t検定：連続指標 / R ----
pairsCond = nchoosek(conditions_use, 2);
tt_cont = table();  % 連続指標
tt_R    = table();  % R

for p = 1:size(pairsCond,1)
    cA = pairsCond(p,1); cB = pairsCond(p,2);

    % 連続指標
    A = PF(strcmp(PF.Condition,cA),:);
    B = PF(strcmp(PF.Condition,cB),:);
    [A,B] = align_pairs(A,B);

    for m = 1:numel(METRICS)
        a = A.(METRICS{m}); b = B.(METRICS{m});
        valid = ~isnan(a) & ~isnan(b);
        if any(valid)
            [~,pval,~,st] = ttest(a(valid), b(valid));
            tt_cont = [tt_cont; table(string(cA),string(cB),string(METRICS{m}), ...
                mean(a(valid),'omitnan'), mean(b(valid),'omitnan'), ...
                mean(a(valid)-b(valid),'omitnan'), pval, st.tstat, ...
                'VariableNames', {'CondA','CondB','Metric','MeanA','MeanB','MeanDiff_AminusB','p','t'})]; %#ok<AGROW>
        end
    end

    % R
    Ac = PFc(strcmp(PFc.Condition,cA),:);
    Bc = PFc(strcmp(PFc.Condition,cB),:);
    [Ac,Bc] = align_pairs(Ac,Bc);
    validR = ~isnan(Ac.R) & ~isnan(Bc.R);
    if any(validR)
        [~,pR,~,stR] = ttest(Ac.R(validR), Bc.R(validR));
        tt_R = [tt_R; table(string(cA),string(cB), ...
            mean(Ac.R(validR),'omitnan'), mean(Bc.R(validR),'omitnan'), ...
            mean(Ac.R(validR)-Bc.R(validR),'omitnan'), pR, stR.tstat, ...
            'VariableNames', {'CondA','CondB','RA','RB','RDiff_AminusB','p','t'})]; %#ok<AGROW>
    end
end

% ---- LME（条件が2以上を一括で固定効果、Fish/Cohortはランダム効果）----
S_LME = S;  % sample粒度
S_LME.Condition = categorical(S_LME.Condition);
% 参照カテゴリは先頭
S_LME.Condition = reordercats(S_LME.Condition, conditions_use);

% 連続指標ごとに LME
lme_cont = struct();
for m = 1:numel(METRICS)
    y = METRICS{m};
    if ~all(ismember({'Cohort','Fish',y}, S_LME.Properties.VariableNames))
        continue;
    end
    mdl = fitglme(S_LME, sprintf('%s ~ 1 + Condition + (1|Cohort) + (1|Fish)', y), ...
                  'Distribution','normal','Link','identity');
    lme_cont.(y) = struct('anova', anova(mdl), 'model', mdl);
end

% Projected Normal（cos/sin）
mdl_cos = fitglme(S_LME, 'cosMag ~ 1 + Condition + (1|Cohort) + (1|Fish)', ...
                  'Distribution','normal','Link','identity');
mdl_sin = fitglme(S_LME, 'sinMag ~ 1 + Condition + (1|Cohort) + (1|Fish)', ...
                  'Distribution','normal','Link','identity');

% ---- 出力 ----
OUT = struct();
OUT.conditions          = string(conditions_use);
OUT.N_pairs_common      = numel(commonPairs);
OUT.window              = [winStartSec, winStartSec + winLenSec];
OUT.pairwise.ttest_cont = tt_cont;
OUT.pairwise.ttest_R    = tt_R;
OUT.per_fish_cond.means = PF;
OUT.per_fish_cond.circ  = PFc;
OUT.lme.continuous      = lme_cont;
OUT.lme.cos             = struct('anova', anova(mdl_cos), 'model', mdl_cos);
OUT.lme.sin             = struct('anova', anova(mdl_sin), 'model', mdl_sin);

% 画面サマリ
fprintf('Common paired Fish (Cohort×Fish): %d\n', OUT.N_pairs_common);
disp('--- Paired t-tests (continuous metrics) ---'); disp(tt_cont);
disp('--- Paired t-tests (R) ---'); disp(tt_R);
disp('--- LME anova (cos) ---'); disp(OUT.lme.cos.anova);
disp('--- LME anova (sin) ---'); disp(OUT.lme.sin.anova);

end

%% ================= ヘルパ =================

function PF = per_fish_cond_means(S, metricNames)
% Fish×Condition の窓内平均
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
% Fish×Condition の meanAngle_deg, R, Rayleigh
S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
G = findgroups(S.Pair, S.Condition);
Pair = splitapply(@(x) x(1), S.Pair, G);
Cond = splitapply(@(x) x(1), S.Condition, G);
mu = nan(size(Pair)); R = nan(size(Pair)); Z = nan(size(Pair)); p = nan(size(Pair));
uG = unique(G);
for gi = 1:numel(uG)
    g = uG(gi);
    idx = (G==g);
    ang = deg2rad(S.HeadingMag(idx));
    ang = ang(~isnan(ang));
    if ~isempty(ang)
        C = mean(cos(ang)); Ssin = mean(sin(ang));
        mu(gi) = rad2deg(atan2(Ssin, C));
        R(gi)  = hypot(C,Ssin);
        n      = numel(ang);
        Z(gi)  = n * R(gi)^2;
        p(gi)  = exp(sqrt(1+4*n+4*(n^2 - (n*R(gi))^2)) - (1+2*n));
        p(gi)  = min(max(p(gi),0),1);
    end
end
PFc = table(Pair, Cond, mu, R, Z, p, 'VariableNames', ...
    {'Pair','Condition','meanAngle_deg','R','RayleighZ','RayleighP'});
end

function [A,B] = align_pairs(A,B)
% Pair キーの共通部分で整列
[common, ia, ib] = intersect(A.Pair, B.Pair, 'stable');
A = A(ia,:); B = B(ib,:);
end

function S = add_distance10s(S)
% 10秒距離（後ろ向き窓）を各サンプル行へ
if ~ismember('Distance10s', S.Properties.VariableNames)
    S.Distance10s = nan(height(S),1);
end
keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g);
    [~, ord] = sort(S.Time(idx));
    ii = idx(ord);
    t = S.Time(ii); v = S.Speed(ii);
    if numel(t) < 2, continue; end
    dt = median(diff(t), 'omitnan');
    win = max(1, round(10 / max(dt, eps)));
    step = v .* dt;
    if numel(step) >= win
        dist10 = movsum(step, [win-1, 0], 'Endpoints','discard');
        D = nan(numel(ii),1);
        D(win:end) = dist10;
        S.Distance10s(ii) = D;
    end
end
end

function S = add_tortuosity(S)
% 速度と HeadingMag から軌跡を再構成し、10秒窓トータスを各サンプル行へ
% 代入寸法を常に一致させる（要素数不一致を根絶）
if ~ismember('Tortuosity', S.Properties.VariableNames)
    S.Tortuosity = nan(height(S),1);
elseif height(S.Tortuosity) ~= height(S)
    S.Tortuosity = nan(height(S),1);
end

keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g);
    [~, ord] = sort(S.Time(idx));
    ii = idx(ord);

    if numel(ii) < 2, continue; end

    t  = S.Time(ii);
    v  = S.Speed(ii);
    th = deg2rad(S.HeadingMag(ii));

    dt  = median(diff(t), 'omitnan');
    win = max(1, round(10 / max(dt, eps)));

    % 速度→変位累積
    dx = v .* dt .* cos(th);
    dy = v .* dt .* sin(th);
    X  = cumsum([0; dx]);               % 長さ N+1
    Y  = cumsum([0; dy]);
    step = hypot(diff(X), diff(Y));     % 長さ N (=numel(ii))

    TT = nan(numel(ii),1);              % 出力（各サンプルに1値、先頭側はNaN）

    if numel(step) >= win
        % 評価可能な終端インデックス（step基準）
        idxK = win:numel(step);         % 長さ M = N - win + 1

        % path length（窓内総行路長）
        path10 = movsum(step, [win-1, 0], 'Endpoints','discard');  % 長さ M

        % displacement（窓始点→終点の直線距離）
        disp10 = nan(size(step));
        for k = idxK
            x0 = X(k-win+1); y0 = Y(k-win+1);
            x1 = X(k+1);     y1 = Y(k+1);
            disp10(k) = hypot(x1-x0, y1-y0);
        end
        tor_k = path10 ./ disp10(idxK);      % 長さ M
        tor_k(~isfinite(tor_k)) = NaN;

        % サンプル列へ整合代入（終端サンプルは k+1）
        putIdx = idxK + 1;                    % 長さ M, 範囲  (win+1) .. (N+1) だが N+1は存在しないので注意
        putIdx(putIdx > numel(TT)) = [];      % 越えたぶんは除外
        M2 = numel(putIdx);
        if M2 > 0
            TT(putIdx) = tor_k(1:M2);
        end
    end

    % 念のため要素数を一致させてから代入
    if numel(TT) ~= numel(ii), TT = pad_or_clip(TT, numel(ii)); end
    S.Tortuosity(ii) = TT;
end
end

function v = pad_or_clip(v, N)
v = v(:);
if numel(v) > N
    v = v(1:N);
elseif numel(v) < N
    v(end+1:N) = NaN;
end
end
