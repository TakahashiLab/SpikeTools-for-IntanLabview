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
    meanAngDeg(end+1,1) = rad2deg(atan2(Sbar,Cbar)); 
    Rbar(end+1,1)       = hypot(Cbar,Sbar);          

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
