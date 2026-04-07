function d = compute_deltas_rad(S, condCol, cA, cB)
SjA = S(condCol==cA,:); SjB = S(condCol==cB,:);
[Pa, muA] = fish_mean_angles(SjA);
[Pb, muB] = fish_mean_angles(SjB);
[~, ia, ib] = intersect(Pa, Pb, 'stable');
if isempty(ia), d = []; return; end
d = local_wrapToPi(muA(ia) - muB(ib));
end
