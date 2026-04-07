function OUT = report_magnet_allpdf(csvPath, pdfOut)
% report_magnet_allpdf_v2
% 4条件（beringSea/okhotsksea/controlSea1/rotation90）＋強度（Intensity連続）の一括PDF。
% 時間窓は (1) 5–15分(300–900s), (2) 10–15分(600–900s)。
% NaN/単一水準による LME 落ちを徹底回避（fitglme_safe）。

if nargin < 2 || isempty(pdfOut), pdfOut = 'magnet_all_v2.pdf'; end

% ===== 入力・前処理 =====
T = readtable(csvPath);
needCols = {'Level','Session','Condition','Cohort','Fish','Track','Time', ...
            'Speed','HeadingMag','cosMag','sinMag','TurnRate','Curvature','Intensity_uT'};
for k = 1:numel(needCols)
    assert(ismember(needCols{k}, T.Properties.VariableNames), ...
        'CSVに必要列がありません: %s', needCols{k});
end

Sall = T(strcmp(T.Level,'sample'), :);
Sall.Condition = lower(string(Sall.Condition));
Sall.Condition = strrep(Sall.Condition, "ohotukusea", "okhotsksea");
Sall.Pair = strcat("C", string(Sall.Cohort), "_F", string(Sall.Fish));

% 付加指標（10s距離/10sトータス）
Sall = add_distance10s(Sall);
Sall = add_tortuosity(Sall);

% ===== 設定 =====
winList = [300 900; 600 900];                       % [start end] 秒
conds4  = ["beringsea","okhotsksea","controlsea1","rotation90"];
METRICS = {'Speed','Distance10s','TurnRate','Curvature','Tortuosity'};

% 既存PDF削除
if exist(pdfOut,'file'), delete(pdfOut); end
OUT = struct(); OUT.pdf = string(pdfOut);

for w = 1:size(winList,1)
    winStartSec = winList(w,1); winEndSec = winList(w,2);
    winTag = sprintf('[%gs–%gs]', winStartSec, winEndSec);

    S = Sall(Sall.Time >= winStartSec & Sall.Time <= winEndSec, :);

    %% 4条件：対応あり解析
    S4 = S(ismember(S.Condition, conds4), :);

    % 共通 Pair のみ
    commonPairs4 = [];
    for i = 1:numel(conds4)
        Pi = unique(S4.Pair(S4.Condition==conds4(i)));
        if i==1, commonPairs4 = Pi; else, commonPairs4 = intersect(commonPairs4, Pi); end
    end
    S4 = S4(ismember(S4.Pair, commonPairs4), :);

    PF4  = per_fish_cond_means(S4, METRICS);
    PFc4 = per_fish_cond_circular(S4);

    % 全ペア t 検定
    pairsCond4 = nchoosek(conds4, 2);
    tt4_cont = table(); tt4_R = table();
    for p = 1:size(pairsCond4,1)
        cA = pairsCond4(p,1); cB = pairsCond4(p,2);
        A = PF4(strcmp(PF4.Condition,cA),:);
        B = PF4(strcmp(PF4.Condition,cB),:);
        [A,B] = align_pairs(A,B);
        for m = 1:numel(METRICS)
            a=A.(METRICS{m}); b=B.(METRICS{m}); v=~isnan(a)&~isnan(b);
            if any(v)
                [~,pval,ci,st]=ttest(a(v),b(v));
                dz=st.tstat/sqrt(sum(v));
                tt4_cont=[tt4_cont; table(string(cA),string(cB),string(METRICS{m}), ...
                    mean(a(v),'omitnan'), mean(b(v),'omitnan'), mean(a(v)-b(v),'omitnan'), ...
                    pval, st.tstat, dz, ci(1), ci(2), ...
                    'VariableNames',{'CondA','CondB','Metric','MeanA','MeanB','MeanDiff_AminusB','p','t','dz','CI_lo','CI_hi'})]; %#ok<AGROW>
            end
        end
        Ac=PFc4(strcmp(PFc4.Condition,cA),:); Bc=PFc4(strcmp(PFc4.Condition,cB),:); [Ac,Bc]=align_pairs(Ac,Bc);
        a=Ac.R; b=Bc.R; v=~isnan(a)&~isnan(b);
        if any(v)
            [~,pval,ci,st]=ttest(a(v),b(v)); dz=st.tstat/sqrt(sum(v));
            tt4_R=[tt4_R; table(string(cA),string(cB), mean(a(v),'omitnan'), mean(b(v),'omitnan'), ...
                mean(a(v)-b(v),'omitnan'), pval, st.tstat, dz, ci(1), ci(2), ...
                'VariableNames',{'CondA','CondB','RA','RB','RDiff_AminusB','p','t','dz','CI_lo','CI_hi'})]; %#ok<AGROW>
        end
    end

    % ---- LME（連続）/ Projected Normal（cos/sin）: 安全化 ----
    S4_LME = S4;
    S4_LME.Condition = categorical(S4_LME.Condition);
    presentCats = categories(S4_LME.Condition);
    newOrder = intersect(conds4, presentCats, 'stable'); % 存在カテゴリのみ
    if ~isempty(newOrder)
        S4_LME.Condition = reordercats(S4_LME.Condition, cellstr(newOrder));
    end

    lme4_cont = struct(); lme4_cont_anova = struct();
    for m = 1:numel(METRICS)
        y = METRICS{m};
        req = {'Condition','Cohort','Fish', y};
        [mdl_ok, mdl, A] = fitglme_safe(S4_LME, sprintf('%s ~ 1 + Condition + (1|Cohort) + (1|Fish)', y), req);
        lme4_cont.(y) = mdl; lme4_cont_anova.(y) = A;
    end
    [~, mdl4_cos, anova4_cos] = fitglme_safe(S4_LME, 'cosMag ~ 1 + Condition + (1|Cohort) + (1|Fish)', {'Condition','Cohort','Fish','cosMag'});
    [~, mdl4_sin, anova4_sin] = fitglme_safe(S4_LME, 'sinMag ~ 1 + Condition + (1|Cohort) + (1|Fish)', {'Condition','Cohort','Fish','sinMag'});

    % ---- 図：ペアライン/ R / ローズ ----
    add_page_title(pdfOut, sprintf('4 Conditions %s', winTag));

    for m=1:numel(METRICS)
        f=figure('Color','w','Position',[100 100 900 550]); hold on; box on;
        pairs = unique(PF4.Pair,'stable');
        cmap = lines(max(1,numel(pairs)));
        for iP=1:numel(pairs)
            yvals=nan(1,numel(newOrder));
            for j=1:numel(newOrder)
                row=PF4(strcmp(PF4.Condition,newOrder(j)) & strcmp(PF4.Pair,pairs(iP)),:);
                if ~isempty(row), yvals(j)=row.(METRICS{m}); end
            end
            plot(1:numel(newOrder), yvals,'-o','Color',cmap(min(iP,end),:), ...
                 'MarkerFaceColor',cmap(min(iP,end),:),'LineWidth',0.8);
        end
        mu=zeros(1,numel(newOrder)); se=mu;
        for j=1:numel(newOrder)
            vv=PF4.(METRICS{m})(strcmp(PF4.Condition,newOrder(j)));
            mu(j)=mean(vv,'omitnan'); se(j)=std(vv,'omitnan')/sqrt(sum(~isnan(vv)));
        end
        if ~isempty(newOrder)
            errorbar(1:numel(newOrder), mu, se, 'k-','LineWidth',1.5,'CapSize',10);
            set(gca,'XTick',1:numel(newOrder),'XTickLabel',newOrder,'TickLabelInterpreter','none');
        else
            set(gca,'XTick',[],'XTickLabel',[]);
        end
        ylabel(METRICS{m}); title(sprintf('Paired lines + mean±SEM (%s) %s', METRICS{m}, winTag));
        exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);
    end

    f=figure('Color','w','Position',[100 100 700 500]); hold on; box on;
    if ~isempty(newOrder)
        muR=zeros(1,numel(newOrder));
        for j=1:numel(newOrder)
            Rv=PFc4.R(strcmp(PFc4.Condition,newOrder(j)));
            muR(j)=mean(Rv,'omitnan');
            scatter(j*ones(size(Rv)),Rv,16,'filled','MarkerFaceAlpha',0.4);
        end
        bar(1:numel(newOrder), muR,'FaceAlpha',0.3,'EdgeColor','none');
        set(gca,'XTick',1:numel(newOrder),'XTickLabel',newOrder,'TickLabelInterpreter','none');
    else
        text(0.5,0.5,'No conditions after window filter','HorizontalAlignment','center');
        set(gca,'XTick',[]); xlim([0 1]);
    end
    ylabel('R (集中度)'); title(sprintf('R by Condition %s', winTag));
    exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);

    f=figure('Color','w','Position',[100 100 900 320]);
    if ~isempty(newOrder)
        tiledlayout(1,numel(newOrder),'Padding','compact','TileSpacing','compact');
        for j=1:numel(newOrder)
            nexttile; pax = polaraxes; hold(pax,'on');
            ang = deg2rad(S4.HeadingMag(S4.Condition==newOrder(j))); ang = ang(~isnan(ang));
            if ~isempty(ang), polarhistogram(pax, ang, 30); else, title(pax, 'No data'); end
            title(pax, string(newOrder(j)));
        end
        sgtitle(sprintf('HeadingMag Rose %s', winTag));
    else
        ax = axes('Position',[0 0 1 1],'Visible','off'); %#ok<LAXES>
        text(0.5,0.5,'No condition left for rose plot','HorizontalAlignment','center');
    end
    exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);

    % テーブル
    figTbl(sprintf('Paired t (continuous) %s', winTag), tt4_cont, pdfOut);
    figTbl(sprintf('Paired t (R) %s', winTag), tt4_R, pdfOut);
    figTbl(sprintf('LME anova (cos) %s', winTag), tableify(anova4_cos), pdfOut);
    figTbl(sprintf('LME anova (sin) %s', winTag), tableify(anova4_sin), pdfOut);

    % 返り値
    OUT.windows(w).range = [winStartSec, winEndSec];
    OUT.windows(w).Npairs_4conds = numel(commonPairs4);
    OUT.windows(w).ttest4_cont = tt4_cont;
    OUT.windows(w).ttest4_R    = tt4_R;
    OUT.windows(w).lme4_cos    = anova4_cos;
    OUT.windows(w).lme4_sin    = anova4_sin;

    %% 強度：Intensity を連続量（2次項まで）
    Sg = S( startsWith(S.Condition,"microt") | S.Condition=="controlt1", :);

    % Intensity が無い/NaN のときのフォールバック
    if ~ismember('Intensity_uT', Sg.Properties.VariableNames)
        Sg.Intensity_uT = nan(height(Sg),1);
    end
    if all(~isfinite(Sg.Intensity_uT))
        map = containers.Map( ...
            {'microt10','microt20','microt30','microt50','microt70','microt100','controlt1'}, ...
            [10 20 30 50 70 100 50] ...
        );
        for r=1:height(Sg)
            key = string(Sg.Condition(r));
            if isKey(map,key), Sg.Intensity_uT(r) = map(key); end
        end
    end

    PFg  = per_fish_cond_means(Sg, METRICS);
    PFcg = per_fish_cond_circular(Sg);
    PFg.Intensity  = group_intensity_from_samples(Sg, PFg);
    PFcg.Intensity = group_intensity_from_samples(Sg, PFcg);

    Sg_LME = Sg;
    Sg_LME.Intensity  = double(Sg_LME.Intensity_uT);
    Sg_LME.Intensity2 = Sg_LME.Intensity.^2;

    % 連続指標の LME（安全化）
    polyLME = struct(); polyANOVA = struct();
    for m = 1:numel(METRICS)
        y = METRICS{m};
        req = {'Intensity','Intensity2','Cohort','Fish', y};
        [ok, mdl, A, msg] = fitglme_safe(Sg_LME, sprintf('%s ~ 1 + Intensity + Intensity2 + (1|Cohort) + (1|Fish)', y), req);
        polyLME.(y) = mdl; polyANOVA.(y) = A;
        if ~ok
            add_page_title(pdfOut, sprintf('WARN: LME failed (%s) %s\n%s', y, winTag, msg));
        end
    end

    % Projected Normal（cos/sin）
    [okc, mdl_cos, anova_cos, msgc] = fitglme_safe(Sg_LME, 'cosMag ~ 1 + Intensity + Intensity2 + (1|Cohort) + (1|Fish)', ...
        {'Intensity','Intensity2','Cohort','Fish','cosMag'});
    [oks, mdl_sin, anova_sin, msgs] = fitglme_safe(Sg_LME, 'sinMag ~ 1 + Intensity + Intensity2 + (1|Cohort) + (1|Fish)', ...
        {'Intensity','Intensity2','Cohort','Fish','sinMag'});
    if ~okc, add_page_title(pdfOut, sprintf('WARN: PN-cos failed %s\n%s', winTag, msgc)); end
    if ~oks, add_page_title(pdfOut, sprintf('WARN: PN-sin failed %s\n%s', winTag, msgs)); end

    % 図：各指標の強度応答（点＋2次フィット＋ p 値）
    for m=1:numel(METRICS)
        f=figure('Color','w','Position',[100 100 780 520]); hold on; box on;
        x=PFg.Intensity; y=PFg.(METRICS{m});
        scatter(x,y,18,'filled','MarkerFaceAlpha',0.35);
        [xp,yp,ciLo,ciHi]=polyfit_ci(x,y,2,200);
        plot(xp, yp,'k-','LineWidth',1.8);
        plot(xp, ciLo,'k--', xp, ciHi,'k--','LineWidth',1);
        xlabel('Intensity (µT)'); ylabel(METRICS{m});

        A = polyANOVA.(METRICS{m});
        pI  = get_p_from_anova(A, 'Intensity');
        pI2 = get_p_from_anova(A, 'Intensity2');
        title(sprintf('Intensity–Response: %s %s (p_I=%.3g, p_I^2=%.3g)', METRICS{m}, winTag, pI, pI2));
        exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);
    end

    % 図：R
    f=figure('Color','w','Position',[100 100 780 520]); hold on; box on;
    x=PFcg.Intensity; y=PFcg.R;
    scatter(x,y,18,'filled','MarkerFaceAlpha',0.35);
    [xp,yp,ciLo,ciHi]=polyfit_ci(x,y,2,200);
    plot(xp, yp,'k-','LineWidth',1.8);
    plot(xp, ciLo,'k--', xp, ciHi,'k--','LineWidth',1);
    xlabel('Intensity (µT)'); ylabel('R (集中度)');
    pIc  = get_p_from_anova(anova_cos, 'Intensity');
    pIc2 = get_p_from_anova(anova_cos, 'Intensity2');
    pIs  = get_p_from_anova(anova_sin, 'Intensity');
    pIs2 = get_p_from_anova(anova_sin, 'Intensity2');
    title(sprintf('Intensity–Response: R %s (cos: p_I=%.3g,p_I^2=%.3g | sin: p_I=%.3g,p_I^2=%.3g)', ...
          winTag, pIc, pIc2, pIs, pIs2));
    exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);

    % 表をPDFへ
    figTbl(sprintf('LME poly ANOVA (continuous metrics) %s', winTag), pack_anova_tables(polyANOVA), pdfOut);
    figTbl(sprintf('Projected Normal ANOVA (cos,sin) %s', winTag), ...
           combine_two_tables(tableify(anova_cos), tableify(anova_sin), ["cos","sin"]), pdfOut);

    % 返り値
    OUT.windows(w).poly_anova  = polyANOVA;
    OUT.windows(w).pn_anova_cos= anova_cos;
    OUT.windows(w).pn_anova_sin= anova_sin;
end

fprintf('Done. Wrote PDF: %s\n', pdfOut);
end

%% ========= 安全な LME フィッタ =========
function [ok, mdl, A, msg] = fitglme_safe(T, formula, requiredCols)
ok = false; mdl = []; A = table(); msg = '';
try
    % 必要列が存在する・数値列は有限値のみ・カテゴリは2水準以上を確保
    T2 = T;
    % 欠損除去のためのインデックス
    need = requiredCols(:)';
    miss = false(height(T2),1);
    for i=1:numel(need)
        v = T2.(need{i});
        if iscategorical(v)
            miss = miss | isundefined(v);
        else
            miss = miss | ~isfinite(double(v));
        end
    end
    T2 = T2(~miss, :);

    % 予測子カテゴリが単一水準しかない場合は中止
    % （ここでは Condition のみチェック）
    if ismember('Condition', T2.Properties.VariableNames) && iscategorical(T2.Condition)
        if numel(categories(removecats(T2.Condition))) < 2
            msg = 'Condition has <2 levels after filtering.'; return;
        end
    end

    % 応答列が実質ゼロ長の時も中止
    % 応答名は formula の左辺から抽出
    resp = strtok(formula, ' ~');
    if isempty(T2) || all(~isfinite(double(T2.(resp))))
        msg = 'No usable observations after excluding NaNs.'; return;
    end

    mdl = fitglme(T2, formula);
    A = anova(mdl);
    ok = true;
catch ME
    msg = ME.message;
    ok = false;
    mdl = [];
    A = table();
end
end

%% ========= 可視化・表出力などのヘルパ =========
function add_page_title(pdfOut, titleStr)
f=figure('Color','w','Position',[100 100 800 200]); 
annotation('textbox',[0.05 0.25 0.9 0.5],'String',titleStr, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',18,'FontWeight','bold','EdgeColor','none','BackgroundColor','w');
exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);
end

function PF = per_fish_cond_means(S, metricNames)
pairs = string(S.Pair); conds = string(S.Condition);
[G] = findgroups(pairs, conds); ug = unique(G(G>0)); nG = numel(ug);
varTypes = [{'string','string'}, repmat({'double'},1,numel(metricNames))];
varNames = [{'Pair','Condition'}, metricNames];
PF = table('Size',[nG, numel(varNames)], 'VariableTypes',varTypes, 'VariableNames',varNames);
for ii = 1:nG
    idx = (G == ug(ii));
    PF.Pair(ii)      = pairs(find(idx,1,'first'));
    PF.Condition(ii) = conds(find(idx,1,'first'));
    for m = 1:numel(metricNames)
        PF{ii, 2+m} = mean(S.(metricNames{m})(idx), 'omitnan');
    end
end
end

function PFc = per_fish_cond_circular(S)
pairs = string(S.Pair); conds = string(S.Condition);
[G] = findgroups(pairs, conds); ug = unique(G(G>0)); nG = numel(ug);
varNames = {'Pair','Condition','meanAngle_deg','R','RayleighZ','RayleighP'};
varTypes = {'string','string','double','double','double','double'};
PFc = table('Size',[nG, numel(varNames)], 'VariableTypes',varTypes, 'VariableNames',varNames);
for ii = 1:nG
    idx = (G == ug(ii)); PFc.Pair(ii)=pairs(find(idx,1,'first')); PFc.Condition(ii)=conds(find(idx,1,'first'));
    ang = deg2rad(S.HeadingMag(idx)); ang = ang(isfinite(ang));
    if isempty(ang)
        PFc{ii,3:end} = {NaN, NaN, NaN, NaN};
    else
        C = mean(cos(ang)); Ssin = mean(sin(ang)); mu = atan2(Ssin, C); Rv = hypot(C, Ssin);
        n = numel(ang); Z = n*Rv^2;
        p = exp(sqrt(1+4*n+4*(n^2 - (n*Rv)^2)) - (1+2*n)); p = min(max(p,0),1);
        PFc.meanAngle_deg(ii) = rad2deg(mu);
        PFc.R(ii)             = Rv;
        PFc.RayleighZ(ii)     = Z;
        PFc.RayleighP(ii)     = p;
    end
end
end

function [A,B] = align_pairs(A,B)
[~,ia,ib]=intersect(A.Pair,B.Pair,'stable'); A=A(ia,:); B=B(ib,:);
end

function figTbl(titleStr, tbl, pdfOut)
f=figure('Color','w','Position',[100 100 1000 600]); 
if isempty(tbl)
    uit = uitable(f,'Data', {'No data'}, 'ColumnName', {'Info'}, 'Units','normalized','Position',[0 0 1 1]);
else
    uit = uitable(f,'Data', table2cell(tbl), 'ColumnName', string(tbl.Properties.VariableNames), ...
        'Units','normalized','Position',[0 0 1 1]);
end
uicontrol('Style','text','String',titleStr,'Units','normalized','Position',[0 .94 1 .06], ...
    'BackgroundColor','w','FontWeight','bold','FontSize',12);
exportgraphics(f,pdfOut,'ContentType','vector','Append',true); close(f);
end

function TT = tableify(anv)
TT = anv;
if ~istable(TT)
    try, TT = struct2table(struct(anv)); catch, TT = table(); end
end
end

function TC = combine_two_tables(T1, T2, names)
if isempty(T1) && isempty(T2)
    TC = table(); return;
end
if isempty(T1), T1 = table(); end
if isempty(T2), T2 = table(); end
if ~isempty(T1), T1.Prop = repmat(string(names(1)), height(T1), 1); end
if ~isempty(T2), T2.Prop = repmat(string(names(2)), height(T2), 1); end
TC = [T1; T2];
if ~isempty(TC)
    vars = T1.Properties.VariableNames;
    if ~isempty(vars) && ~ismember('Prop', vars)
        TC = movevars(TC, 'Prop', 'Before', TC.Properties.VariableNames{1});
    end
end
end

function P = pack_anova_tables(polyANOVA)
metrics = fieldnames(polyANOVA);
P = table();
for i=1:numel(metrics)
    A = tableify(polyANOVA.(metrics{i}));
    if isempty(A), continue; end
    A.Prop = repmat(string(metrics{i}), height(A), 1);
    A = movevars(A, 'Prop', 'Before', A.Properties.VariableNames{1});
    P = [P; A]; %#ok<AGROW>
end
end

function p = get_p_from_anova(A, termName)
p = NaN;
try
    if istable(A) && ismember('Term', A.Properties.VariableNames)
        termCol = string(A.Term);
        if ~ismember('pValue', A.Properties.VariableNames), p = NaN; return; end
        idx = contains(lower(termCol), lower(termName));
        if any(idx), p = A.pValue(find(idx,1,'first')); end
    end
catch, p = NaN;
end
end

function [xp,yp,ciLo,ciHi] = polyfit_ci(x,y,deg,nGrid)
x = x(:); y = y(:);
v = isfinite(x) & isfinite(y);
x = x(v); y = y(v);
[x,ord] = sort(x); y = y(ord);
if numel(x) < deg+1
    if isempty(x), xp = (0:1:1)'; else, xp = linspace(min(x), max(x)+eps, max(10,nGrid))'; end
    yp = nan(size(xp)); ciLo=yp; ciHi=yp; return;
end
[p,S] = polyfit(x,y,deg); 
xp = linspace(min(x), max(x), nGrid)'; 
[yp,delta] = polyval(p, xp, S);
ciLo = yp - 1.96*delta; 
ciHi = yp + 1.96*delta;
end

function PFint = group_intensity_from_samples(Sg, PFlike)
PFint = nan(height(PFlike),1);
for r=1:height(PFlike)
    idx = string(Sg.Pair)==string(PFlike.Pair(r)) & string(Sg.Condition)==string(PFlike.Condition(r));
    PFint(r) = mean(Sg.Intensity_uT(idx), 'omitnan');
end
end

%% ========= 10s距離/トータス =========
function S = add_distance10s(S)
if ~ismember('Distance10s', S.Properties.VariableNames)
    S.Distance10s = nan(height(S),1);
elseif height(S.Distance10s) ~= height(S)
    S.Distance10s = nan(height(S),1);
end
keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g); [~,ord] = sort(S.Time(idx)); ii = idx(ord);
    if numel(ii) < 2, continue; end
    t = S.Time(ii); v = S.Speed(ii);
    dt = median(diff(t),'omitnan');
    win = max(1, round(10/max(dt,eps)));
    step = v.*dt;
    D = nan(numel(ii),1);
    if numel(step) >= win
        dist10 = movsum(step,[win-1,0],'Endpoints','discard');
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
    idx = find(G==g); [~,ord]=sort(S.Time(idx)); ii=idx(ord);
    if numel(ii) < 2, continue; end
    t=S.Time(ii); v=S.Speed(ii); th=deg2rad(S.HeadingMag(ii));
    dt = median(diff(t),'omitnan'); win=max(1, round(10/max(dt,eps)));
    dx=v.*dt.*cos(th); dy=v.*dt.*sin(th);
    X=cumsum([0;dx]); Y=cumsum([0;dy]);
    step=hypot(diff(X),diff(Y)); % N
    TT=nan(numel(ii),1);
    if numel(step)>=win
        idxK=win:numel(step);
        path10=movsum(step,[win-1,0],'Endpoints','discard');
        disp10=nan(size(step));
        for k=idxK
            x0=X(k-win+1); y0=Y(k-win+1);
            x1=X(k+1);     y1=Y(k+1);
            disp10(k)=hypot(x1-x0,y1-y0);
        end
        tor_k=path10./disp10(idxK);
        tor_k(~isfinite(tor_k))=NaN;
        putIdx=idxK+1; putIdx(putIdx>numel(TT))=[];
        if ~isempty(putIdx), TT(putIdx)=tor_k(1:numel(putIdx)); end
    end
    S.Tortuosity(ii)=TT;
end
end
