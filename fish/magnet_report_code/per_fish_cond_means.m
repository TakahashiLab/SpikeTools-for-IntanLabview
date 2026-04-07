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
