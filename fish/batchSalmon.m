% MATLABスクリプト
% 引数 dirName に指定されたディレクトリ以下のすべてのサブディレクトリを検索
function batchSalmon(dirName, varargin)

    p = inputParser;
    p.addParamValue('mode', 'spatial', @ischar);
    p.addParamValue('spon', 5, @isnumeric);

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
        otherwise
            spVsdir = 1;
    end

    speed = p.Results.spon;

    spon = 0;

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

                if spVsdir == 2 %both
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4]);
                    gs = gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 0);
                    [~, ~, mr] = plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', [1:4], 'animal', 'fish');

                    if gs > 0.3 || mr > 0.2
                        figure;
                        cntSp = cntSp + 1;
                        fprintf('%d::%d::%s::gs=%1.2f,mr=%1.2f\n', cntSp, cnt, fileName, gs, mr);
                        subplot(2, 2, 1);
                        imagePmap(ratemap, ocmap);
                        subplot(2, 2, 2);
                        gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 1);
                        subplot(2, 2, 3);
                        plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', [1:4], 'animal', 'fish', 'verbose', 2);

                    end

                elseif spVsdir == 1 %spatial
                    figure;
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4]);

                    subplot(2, 2, 1);
                    imagePmap(ratemap, ocmap);
                    subplot(2, 2, 2);
                    gs1 = gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 1);
                    subplot(2, 2, 3);
                    gs2 = gridSCore(moserac(ratemap), 'allen', 1);

                    if gs1 > 0.3 || gs2 > 0.3
                        cntSp = cntSp + 1;
                        fprintf('%d::%d::%s::gs1=%1.2f,gs2=%1.2f\n', cntSp, cnt, fileName, gs1, gs2);
                    end

                elseif spVsdir == 3 %raw
                    figure;
                    subplot(2, 2, 1);
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 1, 'posture', [1:4]);

                    subplot(2, 2, 2);
                    imagePmap(ratemap, ocmap);
                    subplot(2, 2, 3);
                    gs1 = gridSCore(moserac(imresize(ratemap, [58, 58])), 'allen', 1);
                    subplot(2, 2, 4);
                    gs2 = gridSCore(moserac(ratemap), 'allen', 1);

                    if gs1 > 0.3 || gs2 > 0.3
                        cntSp = cntSp + 1;
                        fprintf('%d::%d::%s::gs1=%1.2f,gs2=%1.2f\n', cntSp, cnt, fileName, gs1, gs2);
                    end

                elseif spVsdir == 4 %hdseg
                    figure;
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4], 'hdseg', 1);

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

                            imagePmap(ratemap{j}, ocmap{j}, max_rate);
                        end

                    end

                    subplot(3, 3, 5);
                    [ratemap, ~, ~, ocmap] = pmap2(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, ...
                        'maxDist', 150, 'khz', 31.25, 'spon', spon, 'speed', speed, ...
                        'binwidth', 2.5, 'grid', [10 10], 'verbose', 0, 'posture', [1:4]);
                    imagePmap(ratemap, ocmap);

                else %direction
                    figure;
                    [~, ~, mr] = plot_polar_rate_map(data.ensemble{i, 3}, data.Traj, data.msT, 'fstart', data.fstart, 'khz', 31.25, 'posture', [1:4], 'animal', 'fish', 'verbose', 1);

                    if mr > 0.2
                        cntSp = cntSp + 1;
                        fprintf('%d::%d::%s:mr=%1.2f\n', cntSp, cnt, fileName, mr);

                    end

                end

                cnt = cnt + 1;

            end

        end

    end

    fprintf('spatial cell types=%d / %d\n', cntSp, cnt-1);
end
