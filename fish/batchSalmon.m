% MATLABスクリプト
% 引数 dirName に指定されたディレクトリ以下のすべてのサブディレクトリを検索
function batchSalmon(dirName, varargin)
    posture4d = [1:4]; %tracking points to calculate head directions
    orgP4d = posture4d;

    p = inputParser;
    p.addParamValue('mode', 'spatial', @ischar);
    p.addParamValue('spon', 5, @isnumeric);
    p.addParamValue('dirtype', 'head', @ischar);
    p.addParamValue('shufflen', 100, @isnumeric);
    p.addParamValue('shuffleon', 0, @isnumeric);
    p.addParamValue('shuffletype', 'shift', @ischar);
    p.addParamValue('percentile', 0.95, @isnumeric);
    p.addParamValue('sftype', 1, @isnumeric);

    p.parse(varargin{:});

    switch (lower(p.Results.mode))
        case 'spatial',
            spVsdir = 1;
        case 'direction',
            spVsdir = 0;
        case 'both',
            spVsdir = 2;
        case 'raw',
            spVsdir = 3;
        case 'hdseg',
            spVsdir = 4;
        case 'velocity-vector',
            spVsdir = 5;
        case 'spatialview',
            spVsdir = 6;
        otherwise
            spVsdir = 1;
    end

    speed = p.Results.spon;
    dirType = p.Results.dirtype;
    shuffleN = p.Results.shufflen;
    percentile = p.Results.percentile;
    shuffleON = p.Results.shuffleon;
    shuffleType = p.Results.shuffletype;
    SFtype = p.Results.sftype;
    
    spon = 0;

    gTh1 = 0.4;
    gTh2 = 0.4;
    dTh = 0.2;
    pTh = 0.4; %place cell 0.6
    sTh = 0.3;
    fTh1 = 0.5;
    fTh2 = -0.5;

    if speed > 0
        spon = 1;
    end

    % ディレクトリ内のすべてのファイルとサブディレクトリを取得
    files = dir(fullfile(dirName, '**', 'subject*.mat'));

    cnt = 1;
    cntSp = 0; %spatial cell type
    % 各.matファイルをループ処理
    for file = 1:length(files)
        %for file=1:2
        % ファイル名の取得
        fileName = fullfile(files(file).folder, files(file).name);

        % ファイルから変数を読み込み
        data = load(fileName, 'ensemble', 'msT', 'fstart', 'Traj');

        % 変数が存在するかをチェック
        if isfield(data, 'ensemble') && isfield(data, 'msT') && ...
                isfield(data, 'fstart') && isfield(data, 'Traj')

            % ここでカスタム関数を実行
            for i = 1:size(data.ensemble, 1)

                if size(data.Traj, 2) <= 4 && posture4d(1)~=1
                    posture4d = orgP4d - 2;
                else
                    posture4d = orgP4d;
                end

                if spVsdir == 2 %both
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType);
                    gs = gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 0);
                    [~, ~, mr] = plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', [1:4], 'animal', 'fish');

                    [~, IPS, IS] = calcInfoNew(ratemap, ocmap);

                    if gs > gTh || mr > dTh || IPS > pTh || IS < sTh
                        figure;
                        cntSp = cntSp + 1;
                        fprintf('%d::%d::%s::gs=%1.2f,mr=%1.2f,IPS=%1.2f, IS=%1.2f\n', cntSp, cnt, fileName, gs, mr, IPS, IS);
                        subplot(2, 2, 1);
                        imagePmap(ratemap, ocmap);
                        subplot(2, 2, 2);
                        gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 1);
                        subplot(2, 2, 3);
                        plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', posture4d, 'animal', 'fish', 'verbose', 2);

                    end

                elseif spVsdir == 1 %spatial
                    figure;
                    [ratemap, spatial_scale, ~, ocmap, pp, ~, FrIndex] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType);

                    subplot(2, 2, 1);
                    imagePmap(ratemap, ocmap, max(ratemap(:)), pp);
                    subplot(2, 2, 2);
                    gs1 = gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 1);
                    subplot(2, 2, 3);
                    gs2 = gridSCore(moserac(ratemap), 'allen', 1);

                    [~, IPS, IS] = calcInfoNew(ratemap, ocmap);
                    bn = Borderness(ratemap, spatial_scale, ocmap, 0);

                    if shuffleON
                        %shuffle
                        [ratemap, ~, ~, ~, ~, ~, Frs] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                            'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                            'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType, 'shuffle', 1, 'shuffleN', shuffleN);

                        for k = 1:size(ratemap, 1)
                            [~, IPSs(k), ISs(k)] = calcInfoNew(ratemap(k, :, :), ocmap);
                            GS1(k) = gridSCore(moserac(imresize(squeeze(ratemap(k, :, :)), [58, 58])), 'allen', 0);
                            GS2(k) = gridSCore(moserac(squeeze(ratemap(k, :, :))), 'allen', 0);
                        end

                        IPSs = sort(IPSs);
                        ISs = sort(ISs);
                        GS1 = sort(GS1);
                        GS2 = sort(GS2);
                        Frs = sort(Frs(:, SFtype));
                        pTh = IPSs(floor(shuffleN * percentile));
                        sTh = ISs(floor(shuffleN * (1 - percentile)));
                        gTh1 = GS1(floor(shuffleN * percentile));
                        gTh2 = GS2(floor(shuffleN * percentile));
                        fTh1 = Frs(floor(shuffleN * (1 - (1 - percentile) / 2)));
                        fTh2 = Frs(floor(shuffleN * (1 - percentile) / 2));

                        if gTh1 < 0 || gTh2 < 0
                            gTh1 = 10;
                            gTh2 = 10;
                        end

                    end

                    fprintf('%d::%d::%s::gs1=%1.2f,gs2=%1.2f. IPS=%1.2f, IS=%1.2f, FrIndex=%1.2f,Border=%1.2f\n', cntSp, cnt, fileName, gs1, gs2, IPS, IS, FrIndex(SFtype), bn);

                    if gs1 > gTh1 || gs2 > gTh2 || IPS > pTh || IS < sTh || FrIndex(1) > fTh1 || FrIndex(SFtype) < fTh2
                        cntSp = cntSp + 1;

                        if gs1 > gTh1 || gs2 > gTh2
                            fprintf('*Grid*');
                        end

                        if IPS > pTh || IS < sTh
                            fprintf('*Place*');
                        end

                        if FrIndex(SFtype) > fTh1
                            fprintf('*Speed-pos*');
                        elseif FrIndex(SFtype) < fTh2
                            fprintf('*Speed-neg*');
                        end

                    end

                    fprintf('FrIndex=%1.2f,fTh1=%1.2f,fTh2=%1.2f', FrIndex(SFtype), fTh1, fTh2);

                    fprintf('pTh=%1.2f,sTh=%1.2f,gTh1=%1.2f, gTh2=%1.2f\n', pTh, sTh, gTh1, gTh2);

                elseif spVsdir == 3 %raw
                    figure;
                    subplot(2, 2, 1);
                    [ratemap, ~, ~, ocmap, pp] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 1, 'posture', [1:4], 'dirtype', dirType);

                    subplot(2, 2, 2);
                    imagePmap(ratemap, ocmap, max(ratemap(:)), pp);
                    subplot(2, 2, 3);
                    gs1 = gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 1);
                    subplot(2, 2, 4);
                    gs2 = gridSCore(moserac(ratemap), 'allen', 1);

                    [~, IPS, IS] = calcInfoNew(ratemap, ocmap);

                    if gs1 > gTh1 || gs2 > gTh2 || IPS > pTh || IS < sTh
                        cntSp = cntSp + 1;
                        fprintf('%d::%d::%s::gs1=%1.2f,gs2=%1.2f. IPS=%1.2f, IS=%1.2f\n', cntSp, cnt, fileName, gs1, gs2, IPS, IS);
                    end

                elseif spVsdir == 4 %hdseg
                    figure;
                    [ratemap, ~, ~, ocmap, pp] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'hdseg', 1, 'dirtype', dirType);

                    k = [9 6 3 2 1 4 7 8];
                    % 各セルをチェックし、空の場合は0で置き換える
                    max_values = cellfun(@(x) isempty(x) * 0 + ~isempty(x) * max(x(:)), ratemap, 'UniformOutput', false);
                    max_rate = zeros(1, length(max_values));

                    for z = 1:length(max_values)

                        if isempty(max_values{z})
                            max_rate(z) = 0;
                        else
                            max_rate(z) = max(max_values{z});
                        end

                    end

                    max_rate = max(max_rate);

                    % 結果の表示
                    for j = [1:8]
                        subplot(3, 3, k(j));

                        if ~isempty(ratemap{j}) && ~isempty(ocmap{j})
                            imagePmap(ratemap{j}, ocmap{j}, max_rate, pp);
                        end

                    end

                    subplot(3, 3, 5);
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4]);
                    imagePmap(ratemap, ocmap);

                elseif spVsdir == 5 %velocity-vector
                    figure;
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType, 'mode', 'velocity-vector');

                    imagePmap(ratemap, ocmap);

                    [~, IPS, IS] = calcInfoNew(ratemap, ocmap);

                    if shuffleON
                        %shuffle
                        [ratemap, ~, ~, ~, ~, ~, Frs] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                            'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                            'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType, 'shuffle', 1, 'shuffleN', shuffleN, 'mode', 'velocity-vector');

                        for k = 1:size(ratemap, 1)
                            [~, IPSs(k), ISs(k)] = calcInfoNew(ratemap(k, :, :), ocmap);
                        end

                        IPSs = sort(IPSs);
                        ISs = sort(ISs);
                        pTh = IPSs(floor(shuffleN * percentile));
                        sTh = ISs(floor(shuffleN * (1 - percentile)));
                    end

                    fprintf('%d::%d::%s::IPS=%1.2f, IS=%1.2f\n', cntSp, cnt, fileName, IPS, IS);

                    if IPS > pTh || IS < sTh
                        cntSp = cntSp + 1;

                        if IPS > pTh || IS < sTh
                            fprintf('*Velocity-vector*');
                        end

                    end

                    fprintf('pTh=%1.2f,sTh=%1.2f\n', pTh, sTh);

                elseif spVsdir == 6 %spatial view
                    figure;
                    [ratemap, spatial_scale, ~, ocmap, pp, ~, FrIndex] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType, 'mode', 'spatialview');

                    subplot(2, 2, 1);
                   
                    imagePmap(ratemap, ocmap);

                    [~, IPS, IS] = calcInfoNew(ratemap, ocmap);

                    if shuffleON
                        %shuffle
                        [ratemap, ~, ~, ~, ~, ~, Frs] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                            'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                            'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'dirtype', dirType, 'shuffle', 1, 'shuffleN', shuffleN, 'mode', 'spatialview');
                       
                        for k = 1:size(ratemap, 1)
                            [~, IPSs(k), ISs(k)] = calcInfoNew(ratemap(k, :, :), ocmap);
                        end

                        IPSs = sort(IPSs);
                        ISs = sort(ISs);
                        Frs = sort(Frs(:, SFtype));
                        pTh = IPSs(floor(shuffleN * percentile));
                        sTh = ISs(floor(shuffleN * (1 - percentile)));

                        fTh1 = Frs(floor(shuffleN * (1 - (1 - percentile) / 2)));
                        fTh2 = Frs(floor(shuffleN * (1 - percentile) / 2));

                    end

                    fprintf('%d::%d::%s::IPS=%1.2f, IS=%1.2f, FrIndex=%1.2f\n', cntSp, cnt, fileName, IPS, IS, FrIndex(SFtype));

                    if IPS > pTh || IS < sTh || FrIndex(1) > fTh1 || FrIndex(SFtype) < fTh2
                        cntSp = cntSp + 1;

                        if IPS > pTh || IS < sTh
                            fprintf('*View*');
                        end

                        if FrIndex(SFtype) > fTh1
                            fprintf('*Speed-pos*');
                        elseif FrIndex(SFtype) < fTh2
                            fprintf('*Speed-neg*');
                        end

                    end

                    fprintf('pTh=%1.2f,sTh=%1.2f\n', pTh, sTh);

                else %direction
                    figure;
                    
                    [~, ~, mr] = plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', posture4d, 'animal', 'fish', 'verbose', 1, 'dirtype', dirType);

                    if shuffleON
                        %shuffle
                        [~, ~, mrs] = plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', posture4d, 'animal', 'fish', 'dirtype', dirType, 'shuffle', 1, 'shufflen', shuffleN, 'shuffleType', shuffleType);
                        mrs = sort(mrs);
                        dTh = mrs(floor(shuffleN * percentile));
                    end

                    fprintf('mr=%1.2f', dTh);

                    if mr > dTh
                        cntSp = cntSp + 1;
                        fprintf('*Dir*');
                    end

                    fprintf('%d::%d::%s::mr=%1.2f\n', cntSp, cnt, fileName, mr);

                end

                cnt = cnt + 1;

            end

        end

    end

    fprintf('spatial cell types=%d / %d\n', cntSp, cnt - 1);
end
