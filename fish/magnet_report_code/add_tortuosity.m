function S = add_tortuosity(S)
if ~ismember('Tortuosity', S.Properties.VariableNames)
    S.Tortuosity = nan(height(S),1);
elseif height(S.Tortuosity) ~= height(S)
    S.Tortuosity = nan(height(S),1);
end
keys = strcat("C", string(S.Cohort), "_F", string(S.Fish), "_", string(S.Session), "_", string(S.Condition));
G = findgroups(keys);
for g = 1:max(G)
    idx = find(G==g);
    [~, ord] = sort(S.Time(idx)); ii = idx(ord);
    if numel(ii) < 2, continue; end
    t  = S.Time(ii);
    v  = S.Speed(ii);
    th = deg2rad(S.HeadingMag(ii));
    valid = isfinite(t) & isfinite(v) & isfinite(th);
    dt  = median(diff(t), 'omitnan');
    if ~isfinite(dt) || dt <= 0
        continue;
    end
    win = max(1, round(10 / max(dt, eps)));
    TT = nan(numel(ii),1);
    edges = diff([false; valid; false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;
    for s = 1:numel(starts)
        seg = starts(s):stops(s);
        if numel(seg) < win
            continue;
        end
        dx = v(seg) .* dt .* cos(th(seg));
        dy = v(seg) .* dt .* sin(th(seg));
        X  = cumsum([0; dx]);
        Y  = cumsum([0; dy]);
        step = hypot(diff(X), diff(Y));
        idxK = win:numel(step);
        path10 = movsum(step, [win-1, 0], 'Endpoints','discard');
        disp10 = nan(size(step));
        for k = idxK
            x0 = X(k-win+1); y0 = Y(k-win+1);
            x1 = X(k+1);     y1 = Y(k+1);
            disp10(k) = hypot(x1-x0, y1-y0);
        end
        tor_k = path10 ./ disp10(idxK);
        tor_k(~isfinite(tor_k)) = NaN;
        putIdx = seg(idxK);
        TT(putIdx) = tor_k;
    end
    S.Tortuosity(ii) = TT;
end
end
