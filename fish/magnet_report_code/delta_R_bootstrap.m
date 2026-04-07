function dR = delta_R_bootstrap(S, conditions_in, B, useParallel)
pairs = nchoosek(conditions_in,2);
dR = table();
for p=1:size(pairs,1)
    A  = pairs(p,1); Bc = pairs(p,2);
    Ra = one_R_global(S, A);
    Rb = one_R_global(S, Bc);
    if isnan(Ra) || isnan(Rb)
        dR = [dR; table(string(A),string(Bc),NaN,NaN,NaN,NaN, ...
            'VariableNames',{'CondA','CondB','Delta','CI_lo','CI_hi','p_boot'})];
        continue
    end
    fish = unique(S.Pair(ismember(S.Condition,[A;Bc])));
    n = numel(fish);
    D = zeros(B,1);
    if useParallel && license('test','Distrib_Computing_Toolbox')
        parfor b=1:B
            pick = fish(randi(n,n,1));
            Ra_b = one_R_from_pairs(S, A, pick);
            Rb_b = one_R_from_pairs(S, Bc, pick);
            D(b)  = Ra_b - Rb_b;
        end
    else
        for b=1:B
            pick = fish(randi(n,n,1));
            Ra_b = one_R_from_pairs(S, A, pick);
            Rb_b = one_R_from_pairs(S, Bc, pick);
            D(b)  = Ra_b - Rb_b;
        end
    end
    dR = [dR; table(string(A),string(Bc), Ra-Rb, quantile(D,0.025), quantile(D,0.975), ...
        mean((D - (Ra-Rb))>=0), 'VariableNames',{'CondA','CondB','Delta','CI_lo','CI_hi','p_boot'})]; %#ok<AGROW>
end
end
