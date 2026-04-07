function Z = per_fish_means_cos_sin(S, cA, cB)
Pa = unique(S.Pair(S.Condition==cA));
Pb = unique(S.Pair(S.Condition==cB));
P  = intersect(Pa,Pb,'stable');
if isempty(P), Z = table(); return; end
Z = table('Size',[0 5], ...
          'VariableTypes', {'string','double','double','double','double'}, ...
          'VariableNames', {'Pair','cosA','cosB','sinA','sinB'});
for k=1:numel(P)
    p = P(k);
    Sa = S(S.Pair==p & S.Condition==cA,:);
    Sb = S(S.Pair==p & S.Condition==cB,:);
    row = table( string(p), ...
                 mean(Sa.cosMag,'omitnan'), mean(Sb.cosMag,'omitnan'), ...
                 mean(Sa.sinMag,'omitnan'), mean(Sb.sinMag,'omitnan'), ...
                 'VariableNames', {'Pair','cosA','cosB','sinA','sinB'});
    Z = [Z; row]; %#ok<AGROW>
end
end
