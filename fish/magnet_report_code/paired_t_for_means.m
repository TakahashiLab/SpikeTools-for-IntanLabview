function [row_cos, row_sin] = paired_t_for_means(Z, cA, cB)
row_cos = local_paired_row(Z.cosA, Z.cosB, cA, cB, "cos");
row_sin = local_paired_row(Z.sinA, Z.sinB, cA, cB, "sin");
end

function row = local_paired_row(a, b, cA, cB, varName)
valid = isfinite(a) & isfinite(b);
n = sum(valid);
if n < 2
    row = table(string(cA), string(cB), string(varName), mean(a(valid),'omitnan'), mean(b(valid),'omitnan'), ...
        mean(a(valid)-b(valid),'omitnan'), NaN, NaN, NaN, NaN, NaN, n, ...
        'VariableNames', {'CondA','CondB','Var','MeanA','MeanB','MeanDiff','p','t','dz','CI_lo','CI_hi','N'});
    return;
end
[~,p,ci,st] = ttest(a(valid), b(valid));
dz = st.tstat / sqrt(n);
row = table(string(cA), string(cB), string(varName), mean(a(valid),'omitnan'), mean(b(valid),'omitnan'), ...
      mean(a(valid)-b(valid),'omitnan'), p, st.tstat, dz, ci(1), ci(2), n, ...
      'VariableNames', {'CondA','CondB','Var','MeanA','MeanB','MeanDiff','p','t','dz','CI_lo','CI_hi','N'});
end
