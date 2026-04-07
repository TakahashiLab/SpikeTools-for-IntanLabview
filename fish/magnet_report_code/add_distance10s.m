function S = add_distance10s(S)
if ~ismember('Distance10s', S.Properties.VariableNames)
    S.Distance10s = nan(height(S),1);
elseif height(S.Distance10s) ~= height(S)
    S.Distance10s = nan(height(S),1);
end
keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g);
    [~, ord] = sort(S.Time(idx)); ii = idx(ord);
    if numel(ii) < 2, continue; end
    t = S.Time(ii); v = S.Speed(ii);
    valid = isfinite(t) & isfinite(v);
    dt = median(diff(t), 'omitnan');
    if ~isfinite(dt) || dt <= 0
        continue;
    end
    win = max(1, round(10 / max(dt, eps)));
    D = nan(numel(ii),1);
    edges = diff([false; valid; false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;
    for k = 1:numel(starts)
        seg = starts(k):stops(k);
        if numel(seg) < win
            continue;
        end
        step = v(seg) .* dt;
        dist10 = movsum(step, [win-1, 0], 'Endpoints','discard');
        D(seg(win:end)) = dist10;
    end
    S.Distance10s(ii) = D;
end
end
