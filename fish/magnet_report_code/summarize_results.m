function txt = summarize_results(alpha, conds, ws, wl, tt_cont, tt_R, tt_proj, lme_cont, anova_cos, anova_sin, Rgrp, dR, targetTbl, uniformTbl, orderUsed, prevUsed, cohortGroupUsed, cohortGroupNames, modelMethodUsed, perFishRayleighTbl, circPair)
lines = {};
lines{end+1} = sprintf('Window %gs-%gs | alpha=%.3f', ws, finite_end(ws, wl), alpha);
lines{end+1} = '=== Summary Guide / 読み方 ===';
lines{end+1} = '1. Pairwise tests compare conditions directly. / ペア比較は条件どうしを直接比較しています。';
lines{end+1} = '2. Continuous metrics include Speed, Distance10s, TurnRate, Curvature, Tortuosity. / 連続指標には Speed, Distance10s, TurnRate, Curvature, Tortuosity が含まれます。';
lines{end+1} = '3. LME sections show model-based tests. / LME の欄はモデルベースの検定結果です。';
lines{end+1} = '4. Order effect asks whether presentation order matters. / 順序効果は提示順の影響を見ています。';
lines{end+1} = '5. PrevCondition asks whether the immediately previous condition matters. / PrevCondition は直前条件の影響を見ています。';
lines{end+1} = '6. Cohort-group effect asks whether cohort groups differ. / コホート群効果は群どうしの違いを見ています。';
lines{end+1} = '';

lines{end+1} = '=== Continuous Pairwise Tests / 連続指標のペア比較 ===';
if ~is_empty_table(tt_cont)
    sig = pick_sig(tt_cont, alpha);
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s vs %s | %s: delta=%.3g, p=%.3g%s', ...
                sig.CondA(i), sig.CondB(i), sig.Metric(i), sig.MeanDiff_AminusB(i), ...
                pick_p(sig,i), pick_padj(sig,i));
        end
    else
        lines{end+1} = '  No significant pairwise result for continuous metrics / 連続指標の有意なペア差は見つかりませんでした。';
    end
else
    lines{end+1} = '  Continuous pairwise test table is empty / 連続指標のペア比較表は空です。';
end

lines{end+1} = sprintf('=== Continuous Model / 連続指標のモデル解析 (%s) ===', upper(char(modelMethodUsed)));
metricNames = fieldnames(lme_cont);
if isempty(metricNames)
    lines{end+1} = '  No continuous model result / 連続指標のモデル結果はありません。';
else
    for i = 1:numel(metricNames)
        name = metricNames{i};
        A = [];
        if isfield(lme_cont.(name), 'anova')
            A = lme_cont.(name).anova;
        end
        if isempty(A)
            lines{end+1} = sprintf('  %s: no model result', name);
        else
            line = sprintf('  %s: Condition p=%.3g', name, pick_p_from_anova(A, 'Condition'));
            if orderUsed
                line = sprintf('%s | Order p=%.3g | Condition x Order p=%.3g', line, ...
                    pick_p_from_anova(A, 'OrderInSession'), pick_p_from_anova(A, 'Condition:OrderInSession'));
            end
            if prevUsed
                line = sprintf('%s | PrevCondition p=%.3g', line, pick_p_from_anova(A, 'PrevCondition'));
            end
            if cohortGroupUsed
                line = sprintf('%s | CohortGroup p=%.3g | Condition x CohortGroup p=%.3g', line, ...
                    pick_p_from_anova(A, 'CohortGroup'), pick_p_from_anova(A, 'Condition:CohortGroup'));
            end
            lines{end+1} = line;
        end
    end
end

lines{end+1} = '=== Across-Fish R Pairwise / 魚間集中度 R のペア比較 ===';
if ~is_empty_table(tt_R)
    sig = pick_sig(tt_R, alpha);
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s vs %s: deltaR=%.3g, p=%.3g%s', ...
                sig.CondA(i), sig.CondB(i), sig.RDiff_AminusB(i), pick_p(sig,i), pick_padj(sig,i));
        end
    else
        lines{end+1} = '  No significant pairwise result for R / R の有意なペア差は見つかりませんでした。';
    end
end

lines{end+1} = '=== Fish-Mean cos/sin Pairwise / 各魚平均 cos・sin のペア比較 ===';
if ~is_empty_table(tt_proj)
    sig = pick_sig(tt_proj, alpha);
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s vs %s | %s: delta=%.3g, p=%.3g%s', ...
                sig.CondA(i), sig.CondB(i), sig.Var(i), sig.MeanDiff(i), pick_p(sig,i), pick_padj(sig,i));
        end
    else
        lines{end+1} = '  No significant pairwise result for fish-mean cos/sin / 各魚平均 cos・sin の有意差は見つかりませんでした。';
    end
end

pc = pick_p_from_anova(anova_cos, 'Condition');
ps = pick_p_from_anova(anova_sin, 'Condition');
lines{end+1} = sprintf('=== Angular Model / 角度成分のモデル解析 (%s) ===', upper(char(modelMethodUsed)));
lines{end+1} = sprintf('  Condition effect / 条件主効果: cos p=%.3g | sin p=%.3g', pc, ps);
if orderUsed
    pco = pick_p_from_anova(anova_cos, 'OrderInSession');
    pso = pick_p_from_anova(anova_sin, 'OrderInSession');
    pci = pick_p_from_anova(anova_cos, 'Condition:OrderInSession');
    psi = pick_p_from_anova(anova_sin, 'Condition:OrderInSession');
    lines{end+1} = sprintf('  Order effect / 順序主効果: cos p=%.3g | sin p=%.3g', pco, pso);
    lines{end+1} = sprintf('  Condition x Order / 条件 x 順序: cos p=%.3g | sin p=%.3g', pci, psi);
    lines{end+1} = interpret_order_effect(alpha, pco, pso, pci, psi);
else
    lines{end+1} = '  Order was not modeled / 順序効果はこのモデルに入っていません。';
end

if prevUsed
    pcp = pick_p_from_anova(anova_cos, 'PrevCondition');
    psp = pick_p_from_anova(anova_sin, 'PrevCondition');
    lines{end+1} = sprintf('  PrevCondition effect / 直前条件効果: cos p=%.3g | sin p=%.3g', pcp, psp);
    lines{end+1} = interpret_prev_effect(alpha, pcp, psp);
end

if cohortGroupUsed
    pcc = pick_p_from_anova(anova_cos, 'CohortGroup');
    psc = pick_p_from_anova(anova_sin, 'CohortGroup');
    pcint = pick_p_from_anova(anova_cos, 'Condition:CohortGroup');
    psint = pick_p_from_anova(anova_sin, 'Condition:CohortGroup');
    if numel(cohortGroupNames) == 2
        groupLabelA = string(cohortGroupNames{1});
        groupLabelB = string(cohortGroupNames{2});
        lines{end+1} = sprintf('  Cohort-group effect / コホート群主効果 (%s vs %s): cos p=%.3g | sin p=%.3g', ...
            groupLabelA, groupLabelB, pcc, psc);
        lines{end+1} = sprintf('  Condition x CohortGroup / 条件 x コホート群: cos p=%.3g | sin p=%.3g', pcint, psint);
        lines{end+1} = interpret_cohort_group_effect(alpha, pcc, psc, pcint, psint, groupLabelA, groupLabelB);
    else
        lines{end+1} = sprintf('  Cohort-group effect / コホート群主効果: cos p=%.3g | sin p=%.3g', pcc, psc);
        lines{end+1} = sprintf('  Condition x CohortGroup / 条件 x コホート群: cos p=%.3g | sin p=%.3g', pcint, psint);
        lines{end+1} = sprintf('  Cohort groups / コホート群: %s', strjoin(string(cohortGroupNames), ' | '));
        lines{end+1} = interpret_multi_cohort_group_effect(alpha, pcc, psc, pcint, psint, cohortGroupNames);
    end
end

if ~isempty(Rgrp)
    lines{end+1} = '=== Group Direction Summary / 条件ごとの方向性要約 ===';
    for i = 1:height(Rgrp)
        lines{end+1} = sprintf('  %s: R=%.3g, 95%%CI %.3g-%.3g, mean=%.1f deg, N=%d', ...
            Rgrp.Condition(i), Rgrp.R_group(i), Rgrp.R_CI_lo(i), Rgrp.R_CI_hi(i), Rgrp.mu_group_deg(i), Rgrp.N_fish(i));
    end
end

if ~is_empty_table(dR)
    sig = dR(dR.CI_lo > 0 | dR.CI_hi < 0, :);
    lines{end+1} = '=== Delta R Bootstrap / Delta R ブートストラップ ===';
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s - %s: deltaR=%.3g, CI %.3g-%.3g, p_boot=%.3g', ...
                sig.CondA(i), sig.CondB(i), sig.Delta(i), sig.CI_lo(i), sig.CI_hi(i), sig.p_boot(i));
        end
    else
        lines{end+1} = '  No bootstrap-supported delta R difference / Delta R の有意差は見つかりませんでした。';
    end
end

if ~is_empty_table(uniformTbl)
    sig = uniformTbl(uniformTbl.p < alpha, :);
    lines{end+1} = '=== Rayleigh Test Across Fish / 魚平均方向の Rayleigh 検定 ===';
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s: mean=%.1f deg, Rbar=%.3f, Z=%.2f, p=%.3g, N=%d', ...
                sig.Condition(i), sig.mean_angle_deg(i), sig.Rbar_fish(i), sig.Z(i), sig.p(i), sig.N_fish(i));
        end
    else
        lines{end+1} = '  No clear directional bias / 明確な方向バイアスは見つかりませんでした。';
    end
end

if ~is_empty_table(targetTbl)
    sig = targetTbl(targetTbl.p_one_sided < alpha, :);
    lines{end+1} = '=== V-Test for Target Direction / 目標方向の V 検定 ===';
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s: target=%g deg, meanProj=%.3f, V=%.2f, p=%.3g, meanAngle=%.1f deg, Rbar=%.3f, N=%d', ...
                sig.Condition(i), sig.target_deg(i), sig.meanProjection(i), sig.Vstat(i), ...
                sig.p_one_sided(i), sig.mean_angle_deg(i), sig.Rbar(i), sig.N_fish(i));
        end
    else
        lines{end+1} = '  No significant target-direction bias / 目標方向への有意な偏りは見つかりませんでした。';
    end
end

if ~is_empty_table(perFishRayleighTbl)
    lines{end+1} = '=== Per-Fish Rayleigh / 個体ごとの Rayleigh 検定 ===';
    for j = 1:numel(conds)
        cj = conds(j);
        rows = perFishRayleighTbl(perFishRayleighTbl.Condition == cj, :);
        if ~isempty(rows)
            nSig = sum(rows.p < alpha);
            lines{end+1} = sprintf('  %s: %d/%d fish significant', cj, nSig, height(rows));
        end
    end
end

if ~is_empty_table(circPair)
    sig = circPair(circPair.p_adj < alpha, :);
    lines{end+1} = '=== Circular Pairwise Tests / 円統計の条件間比較 ===';
    if ~isempty(sig)
        for i = 1:height(sig)
            lines{end+1} = sprintf('  %s vs %s: meanDelta=%.1f deg, Rbar=%.3f, Z=%.2f, p=%.3g, p_adj=%.3g, N=%d', ...
                sig.CondA(i), sig.CondB(i), sig.MeanDelta_deg(i), sig.Rbar(i), ...
                sig.Z(i), sig.p(i), sig.p_adj(i), sig.N(i));
        end
    else
        lines{end+1} = '  No significant circular pairwise difference / 円統計で有意な条件差は見つかりませんでした。';
    end
end

txt = string(lines(:));
end

function tf = is_empty_table(T)
tf = isempty(T) || (istable(T) && height(T) == 0);
end

function T2 = pick_sig(Tin, a)
if any(ismember(Tin.Properties.VariableNames, 'p_adj'))
    T2 = Tin(Tin.p_adj < a, :);
elseif any(ismember(Tin.Properties.VariableNames, 'p'))
    T2 = Tin(Tin.p < a, :);
else
    T2 = Tin([], :);
end
end

function s = pick_padj(T, i)
if any(ismember(T.Properties.VariableNames, 'p_adj'))
    s = sprintf(', p_adj=%.3g', T.p_adj(i));
else
    s = '';
end
end

function pv = pick_p(T, i)
if any(ismember(T.Properties.VariableNames, 'p'))
    pv = T.p(i);
else
    pv = NaN;
end
end

function line = interpret_order_effect(alpha, pco, pso, pci, psi)
orderP = [pco pso];
orderP = orderP(isfinite(orderP));
interactionP = [pci psi];
interactionP = interactionP(isfinite(interactionP));
orderSig = ~isempty(orderP) && any(orderP < alpha);
interactionSig = ~isempty(interactionP) && any(interactionP < alpha);
if interactionSig
    line = '  Interpretation / 解釈: order interacts with condition, so presentation order should not be ignored. / 順序と条件の交互作用が示唆されます。';
elseif orderSig
    line = '  Interpretation / 解釈: order has a main effect, but no clear interaction. / 順序の主効果が示唆されます。';
else
    line = '  Interpretation / 解釈: no clear order effect was detected. / 明確な順序効果は見つかりませんでした。';
end
end

function line = interpret_prev_effect(alpha, pcp, psp)
prevP = [pcp psp];
prevP = prevP(isfinite(prevP));
if ~isempty(prevP) && any(prevP < alpha)
    line = '  Interpretation / 解釈: the immediately previous condition likely affects the current response. / 直前条件の影響が示唆されます。';
else
    line = '  Interpretation / 解釈: no clear previous-condition effect was detected. / 直前条件の明確な影響は見つかりませんでした。';
end
end

function line = interpret_cohort_group_effect(alpha, pcc, psc, pcint, psint, groupLabelA, groupLabelB)
groupP = [pcc psc];
groupP = groupP(isfinite(groupP));
interactionP = [pcint psint];
interactionP = interactionP(isfinite(interactionP));
groupSig = ~isempty(groupP) && any(groupP < alpha);
interactionSig = ~isempty(interactionP) && any(interactionP < alpha);
if interactionSig
    line = sprintf('  Interpretation / 解釈: cohort-group difference depends on condition. / %s と %s の差は条件によって変わる可能性があります。', groupLabelA, groupLabelB);
elseif groupSig
    line = sprintf('  Interpretation / 解釈: a cohort-group main effect was detected between %s and %s. / %s と %s の全体差が示唆されます。', ...
        groupLabelA, groupLabelB, groupLabelA, groupLabelB);
else
    line = sprintf('  Interpretation / 解釈: no clear overall difference was detected between %s and %s. / %s と %s の明確な全体差は見つかりませんでした。', ...
        groupLabelA, groupLabelB, groupLabelA, groupLabelB);
end
end

function line = interpret_multi_cohort_group_effect(alpha, pcc, psc, pcint, psint, groupNames)
groupP = [pcc psc];
groupP = groupP(isfinite(groupP));
interactionP = [pcint psint];
interactionP = interactionP(isfinite(interactionP));
groupSig = ~isempty(groupP) && any(groupP < alpha);
interactionSig = ~isempty(interactionP) && any(interactionP < alpha);
groupList = strjoin(string(groupNames), ' / ');
if interactionSig
    line = sprintf('  Interpretation / 解釈: cohort-group difference depends on condition. / コホート群 (%s) の差は条件によって変わる可能性があります。', groupList);
elseif groupSig
    line = sprintf('  Interpretation / 解釈: an overall cohort-group effect was detected. / コホート群 (%s) の全体差が示唆されます。', groupList);
else
    line = sprintf('  Interpretation / 解釈: no clear overall difference was detected across cohort groups. / コホート群 (%s) の明確な全体差は見つかりませんでした。', groupList);
end
end
