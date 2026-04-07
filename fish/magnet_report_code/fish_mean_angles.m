function [P, mu] = fish_mean_angles(Ssub)
if ~ismember('Pair', Ssub.Properties.VariableNames)
    Ssub.Pair = strcat("C", string(Ssub.Cohort), "_F", string(Ssub.Fish));
end
G = findgroups(Ssub.Pair);
cosm = splitapply(@(x) mean(x,'omitnan'), Ssub.cosMag, G);
sinm = splitapply(@(x) mean(x,'omitnan'), Ssub.sinMag, G);
ok = isfinite(cosm) & isfinite(sinm);
mu = atan2(sinm(ok), cosm(ok));
P  = splitapply(@(x) x(1), Ssub.Pair, G); P = P(ok);
end
