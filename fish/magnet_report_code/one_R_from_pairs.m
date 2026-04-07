function Rval = one_R_from_pairs(S, Cname, allowedPairs)
mask = (S.Condition==Cname) & ismember(S.Pair, allowedPairs);
Ssub = S(mask,:);
Rval = one_R_core(Ssub);
end
