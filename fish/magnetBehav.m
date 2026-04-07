function [speed,heading]=magnetBehav(fn)
% Step 1: Load data from the MAT file
data = load(fn);

if isfield(data, 'xydims')
    xydims = data.xydims;
    % Assuming the variable name in the file
elseif isfield(data, 'coordinates')
    xydims = data.coordinates';
else
    error('Neither xydims nor coordinates found in the loaded file.');
end
% Step 3: Analyze the data and plot results
fps = 30;
hover_threshold = 5;  % Example threshold for hovering speed in units per second

% Load data and calculate speed and heading
[speed, heading] = calculate_speed_and_heading(xydims, fps);

% Time in seconds
time = (1:size(speed, 2)) / fps;

% Time after 5 minutes
time_after_5min = time > 9 * 60;

if 0
figure;
hold on;

for fish = 1:9
    % Find indices where speed is below threshold after 5 minutes
    hover_indices = find(speed(fish, :) < hover_threshold & time_after_5min);

    if ~isempty(hover_indices)
        % Plot the heading directions during hovering
        %       plot(time(hover_indices), heading(fish, hover_indices), 'DisplayName', sprintf('Fish %d', fish));
        subplot(3,3,fish)
        polarhistogram(deg2rad(heading(fish,hover_indices)), 'BinEdges', linspace(-pi, pi, 30));
    end
end
end
%xlabel('Time (seconds)');
%ylabel('Heading Direction (degrees)');
%title('Heading Direction During Hovering');
%legend;
%hold off;
return;
%%%%%%%%%%%%%%%
function [speed, heading] = calculate_speed_and_heading(xydim, fps)
    % Settings
   

    num_fish   = 9;   % 魚の数
    if size(xydim,1)==132
        num_points = 7;   % 1匹あたりのトラック点の数
    else
        num_points = 4;
    end
    
    num_samples = size(xydim, 2);
    num_rows    = size(xydim, 1);

    % 出力
    speed   = zeros(num_fish, num_samples-1);
    heading = zeros(num_fish, num_samples-1);

    % --- フレームごとに北向きに回す（8番→9番のベクトルを北(=+90°)に合わせる）---
    % 8番: 行129,130（[x;y]）
    % 9番: 行131,132（[x;y]）
    xydim_rot = xydim; % 回転後を書き込むバッファ

   if size(xydim,1)==132
        pos8 = [xydim(129, 1); xydim(130, 1)];
        pos9 = [xydim(131, 1); xydim(132, 1)];
   else
        pos8 = [xydim(75, 1); xydim(76, 1)];
        pos9 = [xydim(77, 1); xydim(78, 1)];
   end
        vec9_8 = pos8 - pos9;

        angle_to_north = atan2d(vec9_8(2), vec9_8(1));  % x軸基準の角度
        rotation_angle = 90 - angle_to_north;           % 北(+Y)に合わせるための回転角

        c = cosd(rotation_angle);
        sgn = sind(rotation_angle);
        R = [c, -sgn; sgn, c];
    
    for s=1:num_samples
        % このフレーム内の全座標を回転
        for r = 1:2:num_rows
            p = [xydim(r, s); xydim(r+1, s)];
            pr = R * p;
            xydim_rot(r,   s) = pr(1);
            xydim_rot(r+1, s) = pr(2);
        end
    end

    % --- 速度とヘディング計算（回転後の座標で）---
    for fish = 1:num_fish
        base = 2*(fish-1)*num_points;

        for t = 1:(num_samples-1)
            % head: point 1 の現在位置・次時刻位置
            head_curr = [xydim_rot(base+1, t),   xydim_rot(base+2, t)];
            head_next = [xydim_rot(base+1, t+1), xydim_rot(base+2, t+1)];

            % 速度（フレーム間距離 × fps）
            dist = hypot(head_next(1)-head_curr(1), head_next(2)-head_curr(2));
            speed(fish, t) = dist * fps;

            % ヘディング（point2 − point1 の向き）
            p2 = [xydim_rot(base+3, t), xydim_rot(base+4, t)];
            v1 = p2 - head_curr;
            heading(fish, t) = atan2d(v1(2), v1(1));  % [-180,180] deg
        end
    end

    % 平滑化（各行ごとに窓長 ≈ fps）
    for i = 1:size(speed, 1)
        speed(i, :) = smooth(speed(i, :), fps);
    end

return;
