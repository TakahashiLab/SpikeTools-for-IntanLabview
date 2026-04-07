function Rval = one_R_core(Ssub)
if isempty(Ssub), Rval = NaN; return; end
G = findgroups(Ssub.Pair);
cosm = splitapply(@(x) mean(x,'omitnan'), Ssub.cosMag, G);
sinm = splitapply(@(x) mean(x,'omitnan'), Ssub.sinMag, G);
ok = isfinite(cosm) & isfinite(sinm);
cosm = cosm(ok); sinm = sinm(ok);
if isempty(cosm) || isempty(sinm), Rval = NaN; return; end
Rval = hypot(mean(cosm), mean(sinm));
end
