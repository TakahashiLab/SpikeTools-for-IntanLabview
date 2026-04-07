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
    ok = isfinite(cosm) & isfinite(sinm);
    cosm = cosm(ok); sinm = sinm(ok);
    n = numel(cosm);
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
