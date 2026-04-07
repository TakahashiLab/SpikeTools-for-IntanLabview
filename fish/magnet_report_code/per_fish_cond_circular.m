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
    ang = deg2rad(S.HeadingMag(idx)); ang = ang(isfinite(ang));
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
