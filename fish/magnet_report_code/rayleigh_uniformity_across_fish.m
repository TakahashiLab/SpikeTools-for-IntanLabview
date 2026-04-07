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
    Cname(end+1,1)   = cj;               
    Nfish(end+1,1)   = n;                
    meanAng(end+1,1) = rad2deg(mu);      
    Rbar(end+1,1)    = Rb;               
    Zs(end+1,1)      = z;                
    Ps(end+1,1)      = p;                
end
uniformTbl = table(Cname, Nfish, meanAng, Rbar, Zs, Ps, ...
    'VariableNames', {'Condition','N_fish','mean_angle_deg','Rbar_fish','Z','p'});
end
