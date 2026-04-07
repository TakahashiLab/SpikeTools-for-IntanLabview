function Rval = one_R_global(S, Cname)
Ssub = S(S.Condition==Cname,:);
Rval = one_R_core(Ssub);
end
