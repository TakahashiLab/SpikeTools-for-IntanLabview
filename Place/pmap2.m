%%%%%%%%%%%%%%%%%
%Traj/Pos=Posistion(n x posture)
%msT/PosT=time of position data
%posture: Index of head position in Traj array
%fstart: file start time (for logger), not necessary if thethered recording
%maxdist: the maximum distance of arena/maze in a centimeter scale
%khz:sampling rate (intan(takahashilab):25.0, mouselog:31.25)
%spon: speed thresholding
%speed: cm/s
%binwidth: place field binning in a centimeter scale
%grid: kernel size
%verbose: output
%%%%%%%%%%%%%%%%
%%For bird,
%[ratemap,~,~,ocmap]=pmap2(ensemble{i,3},Traj,msT,'fstart',fstart,'animal','any','maxdist',150,'khz',31.25,'spon',1,'speed',5,'binwidth',2.5,'grid',[10 10],'verbose',1);
%%%%%%%%%%%%%%%%
%%For mouse,%
%[ratemap,~,~,ocmap]=pmap2(ensemble{i,3},Traj,msT,'animal','rodent','maxdist',150,'khz',25,'spon',1,'speed',2.5,'binwidth',2.5,'grid',[10
%10]);
%
%
%%for display
%%imagePmap(ratemap,ocmap);

function [rate_map, spatial_scale, SpksShuffle, oc_map, pp, rate_mapB, FrIndex, xdim, ydim, instValues] = pmap2(Spks, Traj, msT, varargin)
    corrFlag = 0;
    year = 2018;
    shuffleN = 100;
    shuffleLag = 20000; %20sec
    SpksShuffle = [];
    clockWise = 1;
    xl = 3;
    yl = 4;

    if size(Traj, 1) == 1 | size(Traj, 2) == 1
        Len = max(Traj);
        Linear = 1;
        xl = 1;
        yl = 1;
    else
        yLen = max(Traj(:, 2)) - min(Traj(:, 2));
        xLen = max(Traj(:, 1)) - min(Traj(:, 1));
        Len = yLen;

        if yLen < xLen
            Len = xLen;
        end

        %     Len=min([xLen yLen])

    end

    p = inputParser;
    p.addParamValue('fstart', 0, @isnumeric);
    p.addParamValue('linear', 0, @isnumeric);
    p.addParamValue('maxdist', 400, @isnumeric);
    p.addParamValue('animal', 'any', @ischar);
    p.addParamValue('khz', 31.25, @isnumeric);
    p.addParamValue('spon', 0, @isnumeric);
    p.addParamValue('speed', 2.5, @isnumeric);
    p.addParamValue('shuffle', 0, @isnumeric);
    p.addParamValue('shufflen', 100, @isnumeric);
    p.addParamValue('shuffletype', 'shift', @ischar);

    p.addParamValue('verbose', 0, @isnumeric);
    p.addParamValue('posture', [3 4 5 6], @isvector);
    p.addParamValue('binwidth', 2.5, @isnumeric);
    p.addParamValue('std', 3, @isnumeric);
    p.addParamValue('xdim', [], @isvector);
    p.addParamValue('ydim', [], @isvector);
    p.addParamValue('grid', [], @isvector);
    p.addParamValue('spatialscale', -1, @isnumeric);
    p.addParamValue('mode', 'spatial', @ischar);
    p.addParamValue('dirtype', 'head', @ischar);

    p.parse(varargin{:});

    fstart = p.Results.fstart;
    Linear = p.Results.linear;
    maxdist = p.Results.maxdist;
    animal = p.Results.animal;
    kHz = p.Results.khz;
    spON = p.Results.spon;
    ThS = p.Results.speed;
    shuffle = p.Results.shuffle;
    shuffleN = p.Results.shufflen;
    shuffleType = p.Results.shuffletype;
    verbose = p.Results.verbose;
    posture = p.Results.posture;
    BinWidthCm = p.Results.binwidth;
    binside = p.Results.std;
    spatial_scale = p.Results.spatialscale;
    calcMode = p.Results.mode;
    dirType = p.Results.dirtype;

    xdim = p.Results.xdim;
    ydim = p.Results.ydim;
    GRID = p.Results.grid;

    cmPerPixel = maxdist / Len; %

    switch (calcMode)
        case 'hdseg',
            calcMode = 1;
        case 'velocity-vector',
            calcMode = 2;
        case 'spatialview',
            calcMode = 3;
        case 'spatial',
            calcMode = 0;
        otherwise
            calcMode = 0;
    end

    %%%%%%%%%%%%%%
    %%%For bird, large arena
    %%'animal','bird','maxdist',120
    %%
    %%if spikegadgets
    %%%'kHz',30
    %%if small arena
    %%'maxdist',80
    %%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%For rat,
    %%%if Linear, then 'maxdist',400
    %%%else 'maxdist',160

    if size(Traj, 1) ~= size(msT, 1)

        msFPS = median(diff(msT));
        seq = 1:length(Traj);
        seq = seq - 1;
        msT = seq * msFPS + msT(1);
        msT = msT';
    end

    switch lower(animal)
        case 'any'
            cmPerPixel = maxdist / Len; %120cm circle;
        case 'rodent'
            cmPerPixel = maxdist / Len; %
            msT = msT / kHz;
    end

    if spatial_scale < 0
        spatial_scale = cmPerPixel / BinWidthCm;
    end

    FPS = floor(1 / (median(diff(msT)) / 1000));
    fs_video = FPS;
    msFPS = floor(1 / FPS * 1000);

    Spks = ceil(Spks / kHz) + fstart; %msec

    Traj = Traj';
    x = Traj(posture(1), :);
    y = Traj(posture(2), :);
    x2 = Traj(posture(3), :);
    y2 = Traj(posture(4), :);
    Traj = Traj';

    if Linear
        dTraj = diff(Traj);
        ind = find(abs(dTraj) > max(Traj) / 2);
        dTraj(ind) = dTraj(ind) - sign(dTraj(ind)) * max(Traj);
        movement = sqrt(sum(dTraj .^ 2, 2)) * cmPerPixel; %
    else
        % movement = sqrt(sum(diff(Traj(:, posture(1:2))) .^ 2, 2)) * cmPerPixel; % 1:2 -> 5:6

        % First, smoothed position was obtained by Matlab’s smooth function using a width of 0.5 s.
        movX = diff(movmean(Traj(:, posture(1)), 500 / msFPS)) * cmPerPixel;
        movY = diff(movmean(Traj(:, posture(2)), 500 / msFPS)) * cmPerPixel;
    end
    
    headdir = atan2d(y - y2, x - x2)';
    
    % Second, speed was calculated independently in the x and y directions, smoothed by Matlab’s smoothing function with a width of 0.8 s.
    speedX = movmean(movX .* 1000 / msFPS, 800 / msFPS)'; % cm/sec
    speedY = movmean(movY .* 1000 / msFPS, 800 / msFPS)'; % cm/sec
    % Lastly, the animal’s running speed was calculated as the combination of speed in the x and y directions.
    speed = sqrt(speedX .^ 2 + speedY .^ 2);

    speedX2=[speedX(msFPS:end) speedX(1:msFPS-1)];
    speedY2=[speedY(msFPS:end) speedY(1:msFPS-1)];

    movedir = atan2d(speedY2 - speedY, speedX2 - speedX)';
    hovering = find(speed < 5)'; %5cm/s

    if hovering(1) == 1
        hovering(1) = [];
    end

    %hovering
    movedir(hovering) = headdir(hovering);

    deltadir = movedir - headdir(1:end-1);
    deltadir(deltadir > 180) = deltadir(deltadir > 180) - 360;
    deltadir(deltadir < -180) = deltadir(deltadir <- 180) + 360;

    switch (dirType)
        case 'move',
            headdir = movedir;

        case 'diff',
            headdir = deltadir;
    end

    %analyzed if the speed is more than 5cm/s
    if spON
        Good = find(speed >= ThS);
        Traj = Traj(Good, :);
        msT = msT(Good);
        speed = speed(Good);
        headdir=headdir(Good);
    end

    if Linear
        x = Traj;
        y = ones(size(Traj));
    else
        x = ceil(Traj(:, xl));
        y = ceil(Traj(:, yl));
    end

    StartTraj = msT(1);
    EndTraj = msT(end);
    Spks = Spks(find(Spks > StartTraj & Spks < EndTraj));
   
    
    if isempty(Spks)
        binside = 2.5;

        if Linear
            dim = binside;
            xdim = min(x):spatial_scale ^ -1 * dim:max(x); %edges of x and y
            dummy = histcounts(x, xdim);
        else
            dummy = Occupancy(x, y, spatial_scale, binside, fs_video);
        end

        rate_map = zeros(size(dummy));
        spatial_scale = [];
        SpksShuffle = [];
        oc_map = ones(size(dummy));
        rate_mapB = [];
        FrIndex = [];
        return;
    end

    %head direction segmentation
    if calcMode == 1
        % 進行方向を8つのセグメントに分割
        direction_segments = -180:45:180;
        num_segments = length(direction_segments) - 1;
        % 方向セグメントの定義
        % 22.5度から開始し、45度刻みでセグメントを定義

        direction_segments = direction_segments + 22.5;

        % ヘッドディレクションデータの範囲調整
        % -157.5° 以下を 202.5° 以上に変換
        headdir(headdir < -157.5) = headdir(headdir < -157.5) + 360;

        % 各方向セグメントに対応するスパイクデータを格納するためのセル配列を初期化
        spk_x_segments = cell(1, num_segments);
        spk_y_segments = cell(1, num_segments);

        % 全位置データを方向セグメントに分類
        % 各位置データポイントのセグメントインデックスを計算
        [~, ~, segment_indices] = histcounts(headdir, direction_segments);

        % 各方向セグメントごとの座標を格納するためのセル配列
        x_segments = cell(1, num_segments);
        y_segments = cell(1, num_segments);

        % 各セグメントに対して x, y 座標を分類
        for i = 1:num_segments
            segment_mask = segment_indices == i;
            x_segments{i} = Traj(segment_mask, xl);
            y_segments{i} = Traj(segment_mask, yl);
        end

        % 各スパイクに対して方向セグメントを決定し、対応するセル配列に追加
        for i = 1:length(Spks)
            idx = find(msT <= Spks(i), 1, 'last');

            if ~isempty(idx)
                direction = headdir(idx);
                % 方向をセグメントに分類
                for j = 1:num_segments

                    if direction >= direction_segments(j) && direction < direction_segments(j + 1)
                        spk_x_segments{j} = [spk_x_segments{j}; Traj(idx, xl)];
                        spk_y_segments{j} = [spk_y_segments{j}; Traj(idx, yl)];
                        break;
                    end

                end

            end

        end

    end

    %Firing rate contrast Index
    windowSize = 1000; %1s
    Th = nanmedian(speed);

    FrIndex = frInd(Spks, StartTraj, windowSize, EndTraj, msT, Th, speed, deltadir);
    spk_x = [];
    spk_y = [];

    fieldShuffleFlag = 0;

    if shuffle

        beginSpks = msT(1);
        endSpks = msT(end);
        entireLength = ceil(endSpks - beginSpks);
        rp = randperm(entireLength - shuffleLag * 2 + 1);
        rp = rp(1:shuffleN);
        orgRp = shuffleLag:(entireLength - shuffleLag);
        rp = orgRp(rp);
        lenOrgRp = length(orgRp);
        lenSpks = length(Spks);
        SpksShuffle = zeros(shuffleN, lenSpks);

        switch (shuffleType)
            case 'shift',

                for i = 1:shuffleN
                    spkss = Spks + rp(i);
                    topSpks = find(spkss >= endSpks);
                    spkss(topSpks) = spkss(topSpks) - endSpks + beginSpks;
                    SpksShuffle(i, :) = sort(spkss);
                end

            case 'bootstrap',

                for i = 1:shuffleN
                    bootsamp = sort(randsample(length(Spks), length(Spks), true));
                    SpksShuffle(i, :) = Spks(bootsamp);
                end

        end

        if ~fieldShuffleFlag

            parfor i = 1:shuffleN
                Spks = SpksShuffle(i, :);
                spk_x = [];
                spk_y = [];

                for j = 1:(size(Traj, 1) - 1)
                    SpkCnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);

                    for k = 1:SpkCnt
                        spk_x = [spk_x; Traj(j, xl)];
                        spk_y = [spk_y; Traj(j, yl)];
                    end

                end

                if Linear
                    [rate_map(i, :, :)] = rateLinearMap(spk_x, x, spatial_scale, fs_video, binside);
                elseif calcMode == 2 %velocity-vector
                    SpkTimes = [];

                    for j = 1:(size(Traj, 1) - 1)
                        spkcnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);
                        SpkTimes = [SpkTimes spkcnt];
                    end

                    hd = headdir;
                    hd(end) = [];

                    [rate_map(i, :, :)] = velocityVector(SpkTimes, hd, speed, 100, fs_video);
                elseif calcMode == 3 %spatial view

                    [~, ~, ~, oc_map] = ratemap(spk_x, spk_y, x, y, spatial_scale, fs_video, 1);

                    [svp] = svmap(spk_x, spk_y, x, y, spatial_scale, fs_video, oc_map, ...
                        headdir, Traj);

                    sv_x = [];
                    sv_y = [];

                    for j = 1:(size(Traj, 1) - 1)
                        SpkCnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);
                        sv_x = [sv_x; repmat(svp(j, 1), SpkCnt, 1)];
                        sv_y = [sv_y; repmat(svp(j, 2), SpkCnt, 1)];
                    end

                    [occupancy, xdim, ydim] = Occupancy(svp(:, 1), svp(:, 2), spatial_scale, binside, fs_video);

                    spikes = hist3([sv_x, sv_y], 'Edges', {xdim, ydim});
                    rate_map(i, :, :) = SmoothMat(spikes, GRID, binside);
                    oc_map = occupancy;

                else %spatial map
                    [rate_map(i, :, :)] = ratemap(spk_x, spk_y, x, y, ...
                        spatial_scale, fs_video, binside);
                end

                %Firing rate index
                FrIndex(i, :) = frInd(Spks, StartTraj, windowSize, EndTraj, msT, Th, speed, deltadir);
            end

            oc_map = [];
            pp = [];
            rate_mapB = [];

        else %fieldShuffle

            SpkCnt = [];
            spk_x = [];
            spk_y = [];

            for j = 1:(size(Traj, 1) - 1)
                spkcnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);

                if spkcnt

                    for k = 1:spkcnt
                        spk_x = [spk_x; Traj(j, xl)];
                        spk_y = [spk_y; Traj(j, yl)];
                    end

                end

            end

            [rm, ~, ~, oc_map] = ratemap(spk_x, spk_y, x, y, spatial_scale, fs_video, binside);

            for i = 1:shuffleN
                [rate_map(i, :, :)] = fieldShuffleEE(rm, oc_map);
            end

            %for Border cell
            for i = 1:shuffleN
                Spks = SpksShuffle(i, :);
                spk_x = [];
                spk_y = [];

                for j = 1:(size(Traj, 1) - 1)
                    SpkCnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);

                    for k = 1:SpkCnt
                        spk_x = [spk_x; Traj(j, xl)];
                        spk_y = [spk_y; Traj(j, yl)];
                    end

                end

                [rate_mapB(i, :, :)] = ratemap(spk_x, spk_y, x, y, spatial_scale, fs_video, binside);
            end

        end

    else
        SpkCnt = [];
        spk_x = [];
        spk_y = [];
        spk_x_disp = [];
        spk_y_disp = [];
        jitter = [-5:1:5];

        for j = 1:(size(Traj, 1) - 1)
            spkcnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);

            if spkcnt

                for k = 1:spkcnt
                    rp = randperm(11);
                    jitterS = jitter(rp(1));
                    spk_x_disp = [spk_x_disp; Traj(j, xl) + jitterS];
                    spk_y_disp = [spk_y_disp; Traj(j, yl) + jitterS];
                    spk_x = [spk_x; Traj(j, xl)];
                    spk_y = [spk_y; Traj(j, yl)];
                end

            end

        end

        if verbose
            line(Traj(:, xl), Traj(:, yl), 'Color', 'k');
            hold on;
            scatter(spk_x, spk_y, '.r', 'SizeData', 25);

            % 進行方向ベクトルを格納するための配列を初期化
            U = zeros(size(spk_x));
            V = zeros(size(spk_y));

            % 各スパイクの進行方向ベクトルを計算
            IDX = [];

            for i = 1:length(spk_x)
                idx = find(msT <= Spks(i), 1, 'last');
                IDX = [IDX idx];

                if ~isempty(idx)
                    U(i) = cosd(headdir(idx));
                    V(i) = sind(headdir(idx));
                end

            end

            %outl=find(headdir(IDX)>45 | headdir(IDX)<-45 );

            % すべての矢印を一度に描画
            quiver(spk_x, spk_y, U, V, 'AutoScale', 'on', 'Color', 'b');
            %quiver(spk_x(outl), spk_y(outl), U(outl), V(outl), 'AutoScale', 'on', 'Color', 'b');

            axis equal off;
            set(gca, 'xdir', 'normal');
            set(gca, 'ydir', 'reverse');
            maxX = max(spk_x);
            maxY = max(spk_y);
        end

        if Linear
            [rate_map, xdim, oc_map] = rateLinearMap(spk_x, x, spatial_scale, fs_video, binside);

        elseif calcMode == 1
            % 各方向セグメントに対してプレースマップを計算
            pmaps = cell(1, num_segments);

            for i = 1:num_segments

                if ~isempty(spk_x_segments{i})
                    %inactivate omitIsLands function
                    [pmaps{i}, ~, ~, ocmaps{i}] = ratemap(spk_x_segments{i}, spk_y_segments{i}, x_segments{i}, y_segments{i}, spatial_scale, fs_video, binside, -1);
                else
                    pmaps{i} = []; % 空のセグメントの場合
                    ocmaps{i} = [];
                end

            end

            rate_map = pmaps;
            oc_map = ocmaps;

        elseif calcMode == 2 %velocity-vector

            SpkTimes = [];

            for j = 1:(size(Traj, 1) - 1)
                spkcnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);
                SpkTimes = [SpkTimes spkcnt];
            end
            
            hd = headdir';
        
            [rate_map, oc_map] = velocityVector(SpkTimes, hd, speed, 50, fs_video);
            xdim = [];
            ydim = [];

        elseif calcMode == 3 % spatial view
            [~, ~, ~, oc_map] = ratemap(spk_x, spk_y, x, y, spatial_scale, fs_video, 1);

            [svp] = svmap(spk_x, spk_y, x, y, spatial_scale, fs_video, oc_map, ...
                headdir, Traj);

            sv_x = [];
            sv_y = [];

            for j = 1:(size(Traj, 1) - 1)
                SpkCnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);
                sv_x = [sv_x; repmat(svp(j, 1), SpkCnt, 1)];
                sv_y = [sv_y; repmat(svp(j, 2), SpkCnt, 1)];
            end

            [occupancy, xdim, ydim] = Occupancy(svp(:, 1), svp(:, 2), spatial_scale, binside, fs_video);
            %    occupancy = OmitIslands(occupancy);

            spikes = hist3([sv_x, sv_y], 'Edges', {xdim, ydim});
            rate_map = SmoothMat(spikes, GRID, binside);
            oc_map = occupancy;
        else

            if ~isempty(xdim) & ~isempty(GRID)
                [rate_map, xdim, ydim, oc_map] = ratemap(spk_x, spk_y, x, y, ...
                    spatial_scale, ...
                    fs_video, binside, GRID, xdim, ydim);
            elseif ~isempty(GRID)
                [rate_map, xdim, ydim, oc_map] = ratemap(spk_x, spk_y, x, y, ...
                    spatial_scale, ...
                    fs_video, binside, GRID);
            else
                [rate_map, xdim, ydim, oc_map] = ratemap(spk_x, spk_y, x, y, ...
                    spatial_scale, ...
                    fs_video, binside);
            end

        end

        if ~iscell(rate_map)
            pp = pfIdent(rate_map, oc_map, Traj, spatial_scale, binside, 0);
        else
            pp = [];
        end

        %Firing rate index
        [FrIndex, instValues] = frInd(Spks, StartTraj, windowSize, EndTraj, msT, Th, speed, deltadir);

        rate_mapB = [];

    end

    %%%%
    function [FrIndex, instValues] = frInd(Spks, StartTraj, windowSize, EndTraj, msT, Th, speed, direction)
        debug = 0;

        FrHigh = [];
        FrLow = [];
        instSpeeds = [];
        instDegs = [];
        Frs = [];

        dir = circ_ang2rad(direction);
        msFPS = ceil(1 ./ median(diff(msT)) * 1000);

        % windowSize=FPS;
        % for k = StartTraj:windowSize:(EndTraj - windowSize)
        %     fr = sum(Spks > k & Spks < k + windowSize) ./ windowSize .* 1000;
        %     instSpeed = nanmedian(speed(msT > k & msT < k + windowSize));
        %     instSpeeds = [instSpeeds instSpeed];
        %     Frs = [Frs fr];

        %     instDeg = circ_mean(dir(msT > k & msT < k + windowSize));
        %     instDegs = [instDegs instDeg];

        %     if instSpeed > Th
        %         FrHigh = [FrHigh fr];
        %     else
        %         FrLow = [FrLow fr];
        %     end
        % end

        % Instantaneous firing rate was smoothed with a 400 ms-wide Gaussian filter.
        for j = 1:(size(speed, 2))
            spkcnt = sum(Spks >= msT(j) & Spks < msT(j) + msFPS);
            Frs = [Frs spkcnt * 1000 / msFPS];
        end

        instSpeeds = speed;
        Frs = movmean(Frs, 400);

        % FrHigh = sum(FrHigh) / (length(FrHigh) * (windowSize / 1000));
        % FrLow = sum(FrLow) / (length(FrLow) * (windowSize / 1000));

        % FrIndex(1, 1) = (FrHigh - FrLow) / (FrHigh + FrLow);
        % FrIndex(1, 2) = FrHigh;
        % FrIndex(1, 3) = FrLow;

        instSpeeds = instSpeeds';
        Frs = Frs';

        ThSpeed = 5;

        indx = find(instSpeeds < ThSpeed);
        instValues = [instSpeeds Frs];
        instValues(indx, :) = [];

        runningSpeed = instValues(:, 1);
        % Running speed の範囲に基づいてビンエッジを計算
        minSpeed = floor(min(runningSpeed));
        maxSpeed = ceil(max(runningSpeed)) + 1;
        binEdges = minSpeed:2:maxSpeed; % 2 cm/s のビン幅
        binIndices = discretize(runningSpeed, binEdges);
        % inhomogenous coverage,> 1% bin is considered
        c = histcounts(runningSpeed, binEdges);
        idx = find(~ismember(binIndices, find(c ./ sum(c) * 100 > .5)));
        instValues(idx, :) = [];

        FrIndex(1, 4) = corr(instValues(:, 1), instValues(:, 2));

        if debug
            figure;
            plot(instSpeeds)
            hold on;
            plot(Frs);
        end

        return;
