function behav = magnetBehav2(fn, opts)
% magnetBehav2.m
% 幾何北基準で座標を回転し、速度・ヘディングを算出。
% さらに磁場北=0°で再基準化したヘディング、回頭率、曲率、
% 時間窓ごとの円統計（平均角・ベクトル長・Rayleigh検定）、
% トータス／直進性を出力する。
%
% 出力フィールド:
%   behav.speed, heading, headingMag, turnRate, curvature
%   behav.tortuosityWin, behav.straightnessWin
%   behav.winStats (table; 高さ = num_fish * nWins)
%   behav.meta

%% --- オプション設定 ---
if nargin < 2, opts = struct(); end
fps              = getfielddef(opts, 'fps', 30);
smooth_win       = getfielddef(opts, 'smooth_win', fps);
magNorthDeg      = getfielddef(opts, 'magNorthDeg', 0);
window_starts    = getfielddef(opts, 'window_starts_sec', 5*60);
win_len_sec      = getfielddef(opts, 'win_len_sec', inf);
hover_threshold  = getfielddef(opts, 'hover_threshold', []);
posture_dispersion_threshold = getfielddef(opts, 'posture_dispersion_threshold', []);
posture_dispersion_mad_scale = getfielddef(opts, 'posture_dispersion_mad_scale', 20);
heading_jump_threshold_deg = getfielddef(opts, 'heading_jump_threshold_deg', []);

%% --- データ読込 ---
data = load(fn);
if isfield(data, 'xydims')
    xydim = data.xydims;
elseif isfield(data, 'coordinates')
    xydim = data.coordinates';
else
    error('Neither xydims nor coordinates found in file.');
end

%% --- 寸法確認 ---
if size(xydim,1) == 132
    num_points = 7;
    idx8 = [129,130]; idx9 = [131,132];
elseif size(xydim,1) == 78
    num_points = 4;
    idx8 = [75,76]; idx9 = [77,78];
else
    error('Unexpected data shape: rows=%d', size(xydim,1));
end
num_fish = 9;
num_samples = size(xydim,2);
num_rows = size(xydim,1);

%% --- 幾何北基準への回転 ---
pos8 = [xydim(idx8(1),1); xydim(idx8(2),1)];
pos9 = [xydim(idx9(1),1); xydim(idx9(2),1)];
vec9_8 = pos8 - pos9;
angle_to_north = atan2d(vec9_8(2), vec9_8(1));
rotation_angle = 90 - angle_to_north;
R = [cosd(rotation_angle) -sind(rotation_angle);
     sind(rotation_angle)  cosd(rotation_angle)];
xydim_rot = xydim;
for s = 1:num_samples
    for r = 1:2:num_rows
        p = [xydim(r,s); xydim(r+1,s)];
        pr = R*p;
        xydim_rot(r,s)   = pr(1);
        xydim_rot(r+1,s) = pr(2);
    end
end

%% --- 速度・ヘディング（幾何北基準） ---
speed = nan(num_fish, num_samples-1);
heading = nan(num_fish, num_samples-1);
step_dist = nan(num_fish, num_samples-1);
posture_dispersion = nan(num_fish, num_samples);
for f = 1:num_fish
    base = 2*(f-1)*num_points;
    for t = 1:num_samples
        XYf = reshape(xydim_rot(base+1:base+2*num_points, t), 2, [])';
        if any(~isfinite(XYf), 'all')
            continue;
        end
        ctr = mean(XYf, 1);
        posture_dispersion(f,t) = sqrt(mean(sum((XYf - ctr).^2, 2)));
    end
    for t = 1:(num_samples-1)
        head_curr = [xydim_rot(base+1,t), xydim_rot(base+2,t)];
        head_next = [xydim_rot(base+1,t+1), xydim_rot(base+2,t+1)];
        if any(~isfinite(head_curr)) || any(~isfinite(head_next))
            continue;
        end
        dist = hypot(head_next(1)-head_curr(1), head_next(2)-head_curr(2));
        step_dist(f,t) = dist;
        speed(f,t) = dist * fps;
        p2 = [xydim_rot(base+3,t), xydim_rot(base+4,t)];
        if any(~isfinite(p2))
            continue;
        end
        v1 = p2 - head_curr;
        heading(f,t) = atan2d(v1(2), v1(1));
    end
end

posture_thresholds = resolve_posture_thresholds(posture_dispersion, posture_dispersion_threshold, posture_dispersion_mad_scale);
posture_bad_frames = posture_dispersion > posture_thresholds;
heading_jump_bad_idx = false(num_fish, size(heading,2));
if ~isempty(heading_jump_threshold_deg) && isfinite(heading_jump_threshold_deg) && heading_jump_threshold_deg > 0
    heading_jump = abs(wrapTo180(diff(heading, 1, 2)));
    jump_cols = heading_jump > heading_jump_threshold_deg;
    if ~isempty(jump_cols)
        heading_jump_bad_idx(:,1:end-1) = heading_jump_bad_idx(:,1:end-1) | jump_cols;
        heading_jump_bad_idx(:,2:end) = heading_jump_bad_idx(:,2:end) | jump_cols;
    end
end
bad_idx = posture_bad_frames(:,1:end-1) | posture_bad_frames(:,2:end) | heading_jump_bad_idx | ...
    ~isfinite(step_dist) | ~isfinite(heading) | ~isfinite(speed);
speed(bad_idx) = NaN;
heading(bad_idx) = NaN;

for i = 1:num_fish
    if smooth_win > 1
        speed(i,:) = smooth_preserve_nan(speed(i,:), smooth_win);
    end
end

time = (1:size(speed,2)) / fps;

%% --- 磁場北 = 0° 再基準化 ---
T = size(speed,2);
if isscalar(magNorthDeg)
    magNorthDeg = repmat(magNorthDeg, 1, T);
elseif numel(magNorthDeg) ~= T
    error('magNorthDeg must be scalar or match time length.');
end
headingMag = wrapTo180(heading - magNorthDeg);

%% --- 回頭率・曲率 ---
dt = 1/fps;
if T >= 2
    dtheta = diff(headingMag,1,2);
    dtheta = wrapTo180(dtheta);
    turnRate = dtheta / dt;
    spd_mid = (speed(:,1:end-1) + speed(:,2:end)) / 2;
    curvature = abs(turnRate) ./ max(spd_mid, eps);
else
    turnRate = zeros(num_fish,0);
    curvature = zeros(num_fish,0);
end

%% --- 時間窓統計 ---
if ~iscell(window_starts), window_starts = num2cell(window_starts); end
nWins = numel(window_starts);
tortuosityWin = nan(num_fish, nWins);
straightnessWin = nan(num_fish, nWins);
meanAngle = nan(num_fish, nWins);
Rlen = nan(num_fish, nWins);
rayZ = nan(num_fish, nWins);
rayP = nan(num_fish, nWins);

for w = 1:nWins
    t0 = window_starts{w};
    if isfinite(win_len_sec)
        t1 = t0 + win_len_sec;
    else
        t1 = time(end);
    end
    idx = (time >= t0) & (time <= t1);
    if ~any(idx), continue; end

    for f = 1:num_fish
        if any(bad_idx(f,idx))
            continue;
        end
        base = 2*(f-1)*num_points;
        XY = [xydim_rot(base+1,idx)', xydim_rot(base+2,idx)'];
        if size(XY,1) >= 2
            seg = hypot(diff(XY(:,1)), diff(XY(:,2)));
            pathLen = sum(seg, 'omitnan');
            dispLen = hypot(XY(end,1)-XY(1,1), XY(end,2)-XY(1,2));
            if dispLen > 0
                tor = pathLen / dispLen;
                tortuosityWin(f,w) = tor;
                straightnessWin(f,w) = 1/tor;
            end
        end
        ang = deg2rad(headingMag(f,idx));
        ang = ang(~isnan(ang));
        if ~isempty(ang)
            [mu,R] = circ_mean_R(ang);
            meanAngle(f,w) = rad2deg(mu);
            Rlen(f,w) = R;
            n = numel(ang);
            Z = n * R^2;
            p = exp(sqrt(1+4*n+4*(n^2 - (n*R)^2)) - (1+2*n));
            rayZ(f,w) = Z;
            rayP(f,w) = min(max(p,0),1);
        end
    end
end

%% --- hover閾値処理 ---
hover_idx = [];
if ~isempty(hover_threshold)
    hover_idx = speed < hover_threshold;
end

%% --- 出力パッケージ ---
behav = struct();
behav.speed = speed;
behav.heading = heading;
behav.headingMag = headingMag;
behav.turnRate = turnRate;
behav.curvature = curvature;
behav.tortuosityWin = tortuosityWin;
behav.straightnessWin = straightnessWin;

% winStats table（必ず num_fish * nWins 行）
Nrows = num_fish * nWins;
meanAngle_col = reshape(meanAngle, [], 1);
Rlen_col = reshape(Rlen, [], 1);
rayZ_col = reshape(rayZ, [], 1);
rayP_col = reshape(rayP, [], 1);
MetaWinIndex = repelem((1:nWins)', num_fish);
FishIndex = repmat((1:num_fish)', nWins, 1);

% すべて同じ長さに調整
meanAngle_col = fixlen(meanAngle_col, Nrows);
Rlen_col      = fixlen(Rlen_col, Nrows);
rayZ_col      = fixlen(rayZ_col, Nrows);
rayP_col      = fixlen(rayP_col, Nrows);
MetaWinIndex  = fixlen(MetaWinIndex, Nrows);
FishIndex     = fixlen(FishIndex, Nrows);

behav.winStats = table(meanAngle_col, Rlen_col, rayZ_col, rayP_col, ...
                       MetaWinIndex, FishIndex, ...
    'VariableNames', {'meanAngle_deg','R','RayleighZ','RayleighP','MetaWinIndex','Fish'});

behav.meta = struct('fps', fps, 'time', time, ...
    'window_starts_sec', cell2mat_safe(window_starts), ...
    'win_len_sec', win_len_sec, ...
    'magNorthDeg', magNorthDeg, ...
    'hover_threshold', hover_threshold, ...
    'heading_jump_threshold_deg', heading_jump_threshold_deg, ...
    'posture_dispersion_threshold', posture_thresholds, ...
    'posture_dispersion_mad_scale', posture_dispersion_mad_scale, ...
    'num_points', num_points);

behav.hover_idx = hover_idx;
behav.bad_idx = bad_idx;
behav.step_dist = step_dist;
behav.posture_dispersion = posture_dispersion;
behav.posture_bad_frames = posture_bad_frames;
behav.heading_jump_bad_idx = heading_jump_bad_idx;
behav.posture_bad_idx = bad_idx;

end

%% --- 補助関数群 ---
function v = getfielddef(S, name, def)
if isfield(S, name), v = S.(name); else, v = def; end
end

function ang = wrapTo180(ang)
ang = mod(ang + 180, 360) - 180;
end

function [mu, R] = circ_mean_R(theta)
theta = theta(~isnan(theta));
if isempty(theta)
    mu = NaN; R = NaN; return;
end
C = mean(cos(theta));
S = mean(sin(theta));
mu = atan2(S, C);
R = hypot(C, S);
end

function y = smooth_preserve_nan(x, win)
x = x(:)';
y = nan(size(x));
valid = ~isnan(x);
if ~any(valid)
    return;
end
edges = diff([false, valid, false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
for k = 1:numel(starts)
    seg_idx = starts(k):stops(k);
    if numel(seg_idx) == 1 || win <= 1
        y(seg_idx) = x(seg_idx);
    else
        y(seg_idx) = smooth(x(seg_idx), win);
    end
end
end

function th = resolve_posture_thresholds(posture_dispersion, posture_dispersion_threshold, posture_dispersion_mad_scale)
num_fish = size(posture_dispersion, 1);
if isempty(posture_dispersion_threshold)
    th = nan(num_fish, 1);
    for f = 1:num_fish
        vals = posture_dispersion(f, :);
        vals = vals(isfinite(vals));
        if isempty(vals)
            th(f) = inf;
            continue;
        end
        medv = median(vals);
        madv = median(abs(vals - medv));
        if madv == 0
            th(f) = max(vals) * 3;
        else
            th(f) = medv + posture_dispersion_mad_scale * madv;
        end
        if ~isfinite(th(f)) || th(f) <= 0
            th(f) = inf;
        end
    end
elseif isscalar(posture_dispersion_threshold)
    th = repmat(posture_dispersion_threshold, num_fish, 1);
elseif numel(posture_dispersion_threshold) == num_fish
    th = posture_dispersion_threshold(:);
else
    error('posture_dispersion_threshold must be empty, scalar, or length num_fish.');
end
end

function v = fixlen(v, N)
v = v(:);
nv = numel(v);
if nv < N
    v(end+1:N) = NaN;
elseif nv > N
    v = v(1:N);
end
end

function arr = cell2mat_safe(ca)
if isempty(ca)
    arr = [];
    return;
end
arr = nan(1,numel(ca));
for k = 1:numel(ca)
    val = ca{k};
    if isnumeric(val) && isscalar(val)
        arr(k) = val;
    else
        arr(k) = NaN;
    end
end
end
