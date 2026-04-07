function circPair = fishmeans_circ_pairwise_tests(S, conditions_in, adjustMethod)
conds = string(conditions_in(:));
if ~ismember('Pair', S.Properties.VariableNames)
    S.Pair = strcat("C", string(S.Cohort), "_F", string(S.Fish));
end
condCol = string(S.Condition);

A_all = strings(0,1); B_all = strings(0,1);
N_all = []; Z_all = []; P_all = []; MeanDeltaDeg_all = []; Rbar_all = [];

pairs = nchoosek(1:numel(conds),2);
for k = 1:size(pairs,1)
    cA = conds(pairs(k,1)); cB = conds(pairs(k,2));
    d = compute_deltas_rad(S, condCol, cA, cB);
    d = d(~isnan(d));
    if isempty(d)
        A_all(end+1,1)=cA; B_all(end+1,1)=cB;
        N_all(end+1,1)=0; Z_all(end+1,1)=NaN; P_all(end+1,1)=NaN;
        MeanDeltaDeg_all(end+1,1)=NaN; Rbar_all(end+1,1)=NaN; %#ok<AGROW>
        continue;
    end
    n = numel(d);
    C = mean(cos(d)); Ssin = mean(sin(d));
    Rbar = hypot(C,Ssin); Z = n * Rbar^2;
    if exist('circ_rtest','file')==2
        [pval, ~] = circ_rtest(d);
    else
        pval = exp(sqrt(1+4*n+4*(n^2 - (n*Rbar)^2)) - (1+2*n));
        pval = min(max(pval,0),1);
    end
    mu = atan2(Ssin, C);
    A_all(end+1,1)=cA; B_all(end+1,1)=cB;
    N_all(end+1,1)=n; Z_all(end+1,1)=Z; P_all(end+1,1)=pval;
    MeanDeltaDeg_all(end+1,1)=rad2deg(mu); Rbar_all(end+1,1)=Rbar; %#ok<AGROW>
end

circPair = table(A_all, B_all, N_all, Z_all, P_all, MeanDeltaDeg_all, Rbar_all, ...
    'VariableNames', {'CondA','CondB','N','Z','p','MeanDelta_deg','Rbar'});

switch lower(string(adjustMethod))
    case "holm"
        circPair.p_adj = holm_adjust(circPair.p);
    case "fdr"
        circPair.p_adj = mafdr(circPair.p,'BHFDR',true);
    otherwise
        circPair.p_adj = circPair.p;
end
end
