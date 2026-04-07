function pv = pick_p_from_anova(anovTbl, term)
pv = NaN;
if isempty(anovTbl), return; end
try, vnames = anovTbl.Properties.VariableNames; catch, pv=NaN; return; end
try
    pCol = "";
    for cand = ["pValue","PValue","p","ProbF","Pr(>F)"]
        if any(strcmpi(vnames, cand))
            pCol = string(vnames{find(strcmpi(vnames, cand), 1)});
            break;
        end
    end
    if strlength(pCol) == 0, return; end
    if any(strcmpi(vnames,'Term'))
        termNames = string(anovTbl.Term);
    elseif any(strcmpi(vnames,'Source'))
        termNames = string(anovTbl.Source);
    elseif ~isempty(anovTbl.Properties.RowNames)
        termNames = string(anovTbl.Properties.RowNames);
    else
        termNames = strings(height(anovTbl),1);
    end
    idx = strcmpi(termNames, term) | contains(lower(termNames), lower(string(term)));
    if any(idx), pv = anovTbl.(pCol)(find(idx,1,'first')); end
catch
    pv = NaN;
end
end
