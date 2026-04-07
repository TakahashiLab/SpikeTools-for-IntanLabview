function OUT = report_magnet_stats_core(csvPath, conditions, ws, wl, opt, outdir)
if ~isfolder(outdir), mkdir(outdir); end

T = readtable(csvPath);
needCols = {'Level','Session','Condition','Cohort','Fish','Track','Time', ...
            'Speed','HeadingMag','cosMag','sinMag','TurnRate','Curvature'};
for k = 1:numel(needCols)
    assert(ismember(needCols{k}, T.Properties.VariableNames), ...
        'CSVに必要列がありません: %s', needCols{k});
end

S = T(strcmp(T.Level,'sample'), :);
if isfinite(wl)
    S = S(S.Time >= ws & S.Time <= ws + wl, :);
else
    S = S(S.Time >= ws, :);
end

S.Condition = lower(string(S.Condition));
conditions_in = lower(string(conditions));
conditions_in = strrep(conditions_in, "ohotukusea", "okhotsksea");
S = S(ismember(S.Condition, conditions_in), :);
if isfield(opt, 'cohortFilter') && ~isempty(opt.cohortFilter)
    cf = double(opt.cohortFilter);
    S = S(S.Cohort >= cf(1) & S.Cohort <= cf(2), :);
end

S = add_distance10s(S);
S = add_tortuosity(S);

S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
commonPairs = [];
for i = 1:numel(conditions_in)
    Pi = unique(S.Pair(S.Condition == conditions_in(i)));
    if i == 1
        commonPairs = Pi;
    else
        commonPairs = intersect(commonPairs, Pi);
    end
end
S = S(ismember(S.Pair, commonPairs), :);

if isfield(opt, 'rotateControlDeg') && isfield(opt, 'rotateControlNames') && ...
        ~isempty(opt.rotateControlDeg) && isfinite(opt.rotateControlDeg) && opt.rotateControlDeg ~= 0
    S = rotate_condition_cos_sin(S, opt.rotateControlDeg, opt.rotateControlNames);
end

METRICS = {'Speed','Distance10s','TurnRate','Curvature','Tortuosity'};
PF = per_fish_cond_means(S, METRICS);
PFc = per_fish_cond_circular(S);
pairsCond = nchoosek(conditions_in, 2);

tt_cont = table();
for p = 1:size(pairsCond, 1)
    cA = pairsCond(p,1);
    cB = pairsCond(p,2);
    A = PF(strcmp(PF.Condition, cA), :);
    B = PF(strcmp(PF.Condition, cB), :);
    [A, B] = align_pairs(A, B);
    for m = 1:numel(METRICS)
        a = A.(METRICS{m});
        b = B.(METRICS{m});
        valid = isfinite(a) & isfinite(b);
        if sum(valid) >= 2
            [~, pval, ci, stats] = ttest(a(valid), b(valid));
            dz = stats.tstat / sqrt(sum(valid));
            tt_cont = [tt_cont; table(string(cA), string(cB), string(METRICS{m}), ...
                mean(a(valid), 'omitnan'), mean(b(valid), 'omitnan'), ...
                mean(a(valid) - b(valid), 'omitnan'), pval, stats.tstat, dz, ci(1), ci(2), ...
                'VariableNames', {'CondA','CondB','Metric','MeanA','MeanB','MeanDiff_AminusB','p','t','dz','CI_lo','CI_hi'})]; %#ok<AGROW>
        end
    end
end

tt_R = table();
for p = 1:size(pairsCond, 1)
    cA = pairsCond(p,1);
    cB = pairsCond(p,2);
    A = PFc(strcmp(PFc.Condition, cA), :);
    B = PFc(strcmp(PFc.Condition, cB), :);
    [A, B] = align_pairs(A, B);
    a = A.R;
    b = B.R;
    valid = isfinite(a) & isfinite(b);
    if sum(valid) >= 2
        [~, pval, ci, stats] = ttest(a(valid), b(valid));
        dz = stats.tstat / sqrt(sum(valid));
        tt_R = [tt_R; table(string(cA), string(cB), ...
            mean(a(valid), 'omitnan'), mean(b(valid), 'omitnan'), ...
            mean(a(valid) - b(valid), 'omitnan'), pval, stats.tstat, dz, ci(1), ci(2), ...
            'VariableNames', {'CondA','CondB','RA','RB','RDiff_AminusB','p','t','dz','CI_lo','CI_hi'})]; %#ok<AGROW>
    end
end

tt_proj = table();
for p = 1:size(pairsCond, 1)
    cA = pairsCond(p,1);
    cB = pairsCond(p,2);
    Z = per_fish_means_cos_sin(S, cA, cB);
    if ~isempty(Z)
        [row_cos, row_sin] = paired_t_for_means(Z, cA, cB);
        tt_proj = [tt_proj; row_cos; row_sin]; %#ok<AGROW>
    end
end

switch lower(string(opt.adjust))
    case "holm"
        if ~isempty(tt_cont), tt_cont.p_adj = holm_adjust(tt_cont.p); end
        if ~isempty(tt_R), tt_R.p_adj = holm_adjust(tt_R.p); end
        if ~isempty(tt_proj), tt_proj.p_adj = holm_adjust(tt_proj.p); end
    case "fdr"
        if ~isempty(tt_cont), tt_cont.p_adj = mafdr(tt_cont.p, 'BHFDR', true); end
        if ~isempty(tt_R), tt_R.p_adj = mafdr(tt_R.p, 'BHFDR', true); end
        if ~isempty(tt_proj), tt_proj.p_adj = mafdr(tt_proj.p, 'BHFDR', true); end
end

S_LME = table();
orderUsed = false;
orderFormulaSuffix = "";
prevUsed = false;
prevFormulaSuffix = "";
groupUsed = false;
groupFormulaSuffix = "";
resolvedCohortGroupNames = {};
resolvedCohortGroupMode = "none";
lme_cont = struct();
mdl_cos = [];
mdl_sin = [];
anova_cos = [];
anova_sin = [];
modelMethodUsed = "off";
if (~isfield(opt, 'doLME') || opt.doLME) && ~strcmpi(string(opt.modelMethod), "off")
    S_LME = build_lme_table(S, METRICS);
    S_LME.Condition = categorical(S_LME.Condition);
    S_LME.Condition = reordercats(S_LME.Condition, cellstr(conditions_in));
    if isfield(opt, 'includeOrderEffect') && opt.includeOrderEffect && ismember('OrderInSession', S_LME.Properties.VariableNames)
        S_LME.OrderInSession = double(S_LME.OrderInSession);
        orderValid = isfinite(S_LME.OrderInSession);
        if sum(orderValid) >= 2 && numel(unique(S_LME.OrderInSession(orderValid))) > 1
            orderUsed = true;
            orderFormulaSuffix = " + OrderInSession";
            if isfield(opt, 'includeOrderInteraction') && opt.includeOrderInteraction
                orderFormulaSuffix = orderFormulaSuffix + " + Condition:OrderInSession";
            end
        end
    end
    if isfield(opt, 'includePrevConditionEffect') && opt.includePrevConditionEffect && ismember('PrevCondition', S_LME.Properties.VariableNames)
        S_LME.PrevCondition = categorical(S_LME.PrevCondition);
        if numel(categories(removecats(S_LME.PrevCondition))) > 1
            prevUsed = true;
            prevFormulaSuffix = " + PrevCondition";
        end
    end
    if isfield(opt, 'includeCohortGroupEffect') && opt.includeCohortGroupEffect && ~isempty(S_LME)
        [S_LME.CohortGroup, resolvedCohortGroupNames, resolvedCohortGroupMode] = build_cohort_groups(S_LME, opt);
        if numel(categories(removecats(S_LME.CohortGroup))) > 1
            groupUsed = true;
            groupFormulaSuffix = " + CohortGroup";
            if isfield(opt, 'includeCohortGroupInteraction') && opt.includeCohortGroupInteraction
                groupFormulaSuffix = groupFormulaSuffix + " + Condition:CohortGroup";
            end
        end
    end

    % In sequence-based cohort groups, OrderInSession is effectively encoded by
    % CohortGroup and often by PrevCondition as well, so the order term becomes
    % rank-deficient and disappears from ANOVA tables. Drop it proactively.
    if orderUsed && groupUsed && strcmpi(string(resolvedCohortGroupMode), "sequence")
        orderUsed = false;
        orderFormulaSuffix = "";
    end

    fixedSuffix = char(orderFormulaSuffix + prevFormulaSuffix + groupFormulaSuffix);
    tryLME = any(strcmpi(string(opt.modelMethod), ["auto","lme"]));
    lmeWorked = false;
    if tryLME
        for m = 1:numel(METRICS)
            y = METRICS{m};
            formula = char(sprintf('%s ~ 1 + Condition%s + (1|Cohort) + (1|Fish)', y, fixedSuffix));
            mdl = fitglme_safe(S_LME, formula, 'Distribution', 'normal', 'Link', 'identity');
            if isempty(mdl)
                lme_cont.(y) = struct('anova', [], 'model', []);
            else
                lmeWorked = true;
                lme_cont.(y) = struct('anova', anova(mdl), 'model', mdl);
            end
        end

        mdl_cos = fitglme_safe(S_LME, char(sprintf('cosMag ~ 1 + Condition%s + (1|Cohort) + (1|Fish)', fixedSuffix)));
        mdl_sin = fitglme_safe(S_LME, char(sprintf('sinMag ~ 1 + Condition%s + (1|Cohort) + (1|Fish)', fixedSuffix)));
        if ~isempty(mdl_cos), anova_cos = anova(mdl_cos); lmeWorked = true; end
        if ~isempty(mdl_sin), anova_sin = anova(mdl_sin); lmeWorked = true; end
    end

    if lmeWorked
        modelMethodUsed = "lme";
    elseif any(strcmpi(string(opt.modelMethod), ["auto","lm"]))
        [lme_cont, mdl_cos, mdl_sin, anova_cos, anova_sin] = fit_fixed_models(S_LME, METRICS, fixedSuffix);
        modelMethodUsed = "lm";
    end
end

Rgrp = across_fish_R(S, conditions_in);

if opt.doBootstrapDeltaR && opt.bootstrapB > 0
    dR = delta_R_bootstrap(S, conditions_in, opt.bootstrapB, opt.useParallel);
else
    dR = table();
end

uniformTbl = table();
if opt.doUniformityTest
    uniformTbl = rayleigh_uniformity_across_fish(S, conditions_in);
end

targetTbl = table();
if opt.doTargetTest
    tgt = normalize_target_angles(opt.targetDirDeg, numel(conditions_in));
    targetTbl = target_direction_vtest(S, conditions_in, tgt);
end

circPair = fishmeans_circ_pairwise_tests(S, conditions_in, opt.adjust);
draw_circ_pairwise_with_raw(S, conditions_in, outdir, ws, wl, 0.05, opt.roseWritePNG, opt.pdfRasterRes);
draw_R_visualization(S, conditions_in, outdir, ws, wl, opt.roseWritePNG, opt.pdfRasterRes);

files = struct();
files.rose_pdf = fullfile(outdir, 'rose_headingMag.pdf');
files.rose_png = ternary(opt.roseWritePNG, fullfile(outdir, 'rose_headingMag.png'), "");
draw_rose_pdf_all(S, conditions_in, files.rose_pdf, files.rose_png, ...
    ws, wl, opt.pdfRasterRes, opt.drawMeanVector, opt.drawTargetOnRose, opt.targetDirDeg, opt.roseBinN, opt.roseWritePNG);

perFishRayleighTbl = table();
perFishFiles = struct();
if opt.doPerFishPlots || opt.doPerFishRayleigh
    [perFishRayleighTbl, perFishFiles] = per_fish_rose_and_tests(S, conditions_in, outdir, ws, wl, opt);
end
files.perFish = perFishFiles;

fishMeansFiles = struct();
if opt.doFishMeansPlots
    fishMeansFiles = draw_fish_means_rose_multi(S, conditions_in, outdir, ws, wl, ...
        opt.pdfRasterRes, opt.drawMeanVector, uniformTbl, opt.alpha, opt.fishMeansTiles, opt.fishMeansWritePNG);
end
files.fishMeans = fishMeansFiles;

statRainFiles = struct();
if opt.doStatsRaincloud
    statRainFiles = draw_paired_raincloud_allmetrics(PF, conditions_in, commonPairs, METRICS, outdir, ws, wl, opt.rainTiles, opt.rainWritePNG, opt.pdfRasterRes);
end
files.statPairs = statRainFiles;

circPairFiles = struct();
if opt.doCircPairPlot
    circPairFiles = plot_circ_pairwise_overview(S, conditions_in, circPair, outdir, ws, wl, opt.alpha, opt.circPairTiles, opt.circPairWritePNG, opt.pdfRasterRes);
end
files.circPairwise = circPairFiles;

lmeFiles = draw_lme_anova_tables(lme_cont, anova_cos, anova_sin, outdir, ws, wl, opt.pdfRasterRes);
files.lmeTables = lmeFiles;

SUM = summarize_results(opt.alpha, conditions_in, ws, wl, ...
    tt_cont, tt_R, tt_proj, lme_cont, anova_cos, anova_sin, Rgrp, dR, targetTbl, uniformTbl, ...
    orderUsed, prevUsed, groupUsed, resolvedCohortGroupNames, modelMethodUsed, perFishRayleighTbl, circPair);
files.summary_txt = fullfile(outdir, 'summary.txt');
save_lines(files.summary_txt, SUM);

OUT = struct();
OUT.conditions = conditions_in;
OUT.N_pairs_common = numel(commonPairs);
OUT.window = [ws, finite_end(ws, wl)];
OUT.pairwise.ttest_cont = tt_cont;
OUT.pairwise.ttest_R = tt_R;
OUT.pairwise.ttest_proj = tt_proj;
OUT.lme.continuous = lme_cont;
OUT.lme.cos = struct('anova', anova_cos, 'model', mdl_cos);
OUT.lme.sin = struct('anova', anova_sin, 'model', mdl_sin);
OUT.modelMethodUsed = modelMethodUsed;
OUT.orderEffect = struct('used', orderUsed, 'formulaSuffix', string(orderFormulaSuffix));
OUT.prevConditionEffect = struct('used', prevUsed, 'formulaSuffix', string(prevFormulaSuffix));
OUT.cohortGroupEffect = struct('used', groupUsed, 'formulaSuffix', string(groupFormulaSuffix), ...
    'mode', resolvedCohortGroupMode, 'cutoff', opt.cohortSplitCutoff, 'names', {resolvedCohortGroupNames});
OUT.groupR = Rgrp;
OUT.deltaR = dR;
OUT.uniformityRayleigh = uniformTbl;
OUT.targetTest = targetTbl;
OUT.perFish.rayleigh = perFishRayleighTbl;
OUT.circPair = circPair;
OUT.files = files;
OUT.summary = SUM;

fprintf('Common paired Fish: %d\n', OUT.N_pairs_common);
disp('--- Paired t-tests (continuous) ---'); if ~isempty(tt_cont), disp(tt_cont); end
disp('--- Paired t-tests (R) ---'); if ~isEmptyTable(tt_R), disp(tt_R); end
disp('--- Fish-mean paired t-tests (cos/sin) ---'); if ~isEmptyTable(tt_proj), disp(tt_proj); end
disp('--- LME anova (cos) ---'); if ~isempty(anova_cos), disp(anova_cos); else, disp('N/A'); end
disp('--- LME anova (sin) ---'); if ~isempty(anova_sin), disp(anova_sin); else, disp('N/A'); end
disp('--- Across-fish R ---'); disp(Rgrp);
if ~isEmptyTable(uniformTbl), disp('--- Rayleigh across fish means ---'); disp(uniformTbl); end
if ~isEmptyTable(perFishRayleighTbl), disp('--- Per-fish Rayleigh ---'); disp(perFishRayleighTbl); end
if ~isEmptyTable(circPair), disp('--- Circular pairwise (fish means, adjusted) ---'); disp(circPair); end
end

function S_LME = build_lme_table(S, metricNames)
groupVars = {'Pair','Condition'};
if ismember('OrderInSession', S.Properties.VariableNames)
    groupVars{end+1} = 'OrderInSession';
end

G = findgroups(S(:, groupVars));
S_LME = table();
S_LME.Pair = splitapply(@(x) x(1), S.Pair, G);
S_LME.Condition = splitapply(@(x) x(1), S.Condition, G);
S_LME.Cohort = splitapply(@(x) x(1), S.Cohort, G);
S_LME.Fish = splitapply(@(x) x(1), S.Fish, G);
if ismember('OrderInSession', S.Properties.VariableNames)
    S_LME.OrderInSession = splitapply(@(x) x(1), S.OrderInSession, G);
    S_LME.PrevCondition = strings(height(S_LME), 1);
    cohorts = unique(S_LME.Cohort);
    for ci = 1:numel(cohorts)
        idx = find(S_LME.Cohort == cohorts(ci));
        [~, ord] = sort(S_LME.OrderInSession(idx));
        ii = idx(ord);
        S_LME.PrevCondition(ii) = "none";
        if numel(ii) > 1
            S_LME.PrevCondition(ii(2:end)) = string(S_LME.Condition(ii(1:end-1)));
        end
    end
end

for i = 1:numel(metricNames)
    m = metricNames{i};
    S_LME.(m) = splitapply(@(x) mean(x, 'omitnan'), S.(m), G);
end

if ismember('cosMag', S.Properties.VariableNames)
    S_LME.cosMag = splitapply(@(x) mean(x, 'omitnan'), S.cosMag, G);
end
if ismember('sinMag', S.Properties.VariableNames)
    S_LME.sinMag = splitapply(@(x) mean(x, 'omitnan'), S.sinMag, G);
end
end

function [model_cont, mdl_cos, mdl_sin, anova_cos, anova_sin] = fit_fixed_models(T, metricNames, fixedSuffix)
model_cont = struct();
mdl_cos = [];
mdl_sin = [];
anova_cos = [];
anova_sin = [];
for i = 1:numel(metricNames)
    y = metricNames{i};
    formula = char(sprintf('%s ~ 1 + Condition%s', y, fixedSuffix));
    mdl = fitlm_safe(T, formula);
    if isempty(mdl)
        model_cont.(y) = struct('anova', [], 'model', []);
    else
        model_cont.(y) = struct('anova', anova(mdl, 'component'), 'model', mdl);
    end
end
mdl_cos = fitlm_safe(T, char(sprintf('cosMag ~ 1 + Condition%s', fixedSuffix)));
mdl_sin = fitlm_safe(T, char(sprintf('sinMag ~ 1 + Condition%s', fixedSuffix)));
if ~isempty(mdl_cos), anova_cos = anova(mdl_cos, 'component'); end
if ~isempty(mdl_sin), anova_sin = anova(mdl_sin, 'component'); end
end

function [cohortGroup, groupNames, modeUsed] = build_cohort_groups(S_LME, opt)
modeUsed = "cutoff";
groupNames = {};
if isfield(opt, 'cohortGroupMode')
    modeUsed = string(opt.cohortGroupMode);
end

if strcmpi(modeUsed, "sequence") && ismember('OrderInSession', S_LME.Properties.VariableNames)
    cohortGroup = strings(height(S_LME),1);
    cohorts = unique(S_LME.Cohort);
    seqKeys = strings(numel(cohorts),1);
    seqLabels = strings(numel(cohorts),1);
    for ci = 1:numel(cohorts)
        idx = find(S_LME.Cohort == cohorts(ci));
        sub = S_LME(idx, {'OrderInSession','Condition'});
        [~, ord] = sort(sub.OrderInSession);
        sub = sub(ord,:);
        [~, firstIdx] = unique(sub.OrderInSession, 'stable');
        sub = sub(firstIdx,:);
        condSeq = string(sub.Condition);
        condSeq = condSeq(strlength(condSeq) > 0);
        seqLabels(ci) = strjoin(condSeq, " -> ");
        seqKeys(ci) = lower(strjoin(condSeq, "|"));
    end
    [~, ia, ic] = unique(seqKeys, 'stable');
    uniqueLabels = seqLabels(ia);
    for ci = 1:numel(cohorts)
        cohortGroup(S_LME.Cohort == cohorts(ci)) = uniqueLabels(ic(ci));
    end
    groupNames = cellstr(uniqueLabels);
    cohortGroup = categorical(cohortGroup, uniqueLabels);
    return;
end

cutoff = 12;
if isfield(opt, 'cohortSplitCutoff') && ~isempty(opt.cohortSplitCutoff)
    cutoff = double(opt.cohortSplitCutoff);
end
groupNames = {'cohort1to12','cohort13plus'};
if isfield(opt, 'cohortGroupNames') && numel(opt.cohortGroupNames) >= 2
    groupNames = opt.cohortGroupNames(1:2);
end
cohortGroup = strings(height(S_LME),1);
cohortGroup(S_LME.Cohort <= cutoff) = string(groupNames{1});
cohortGroup(S_LME.Cohort > cutoff) = string(groupNames{2});
cohortGroup = categorical(cohortGroup, string(groupNames));
groupNames = cellstr(string(groupNames));
end
