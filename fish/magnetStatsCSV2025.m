function results_stats = magnetStatsCSV()

    % mat_files: {filename, condition, cohortID}
    % 例：最初の9匹の群(cohort=1)が複数ファイルで繰り返し登場（同一個体群）、
    %     次のファイル群は別個体群(cohort=2), ... という指定。
    mat_files = {
        '2025-05-28 10-01-55_vmc_checked.mat', 'beringSea',   1;
        '2025-05-28 10-26-38_vmc_checked.mat', 'okhotskSea',  1;
        '2025-05-28 10-47-04_vmc_checked.mat', 'controlSea1', 1;
        '2025-05-28 11-05-32_vmc_checked.mat', 'beringSea',   2;
        '2025-05-28 11-30-09_vmc_checked.mat', 'okhotskSea',  2;
        '2025-05-28 11-49-48_vmc_checked.mat', 'controlSea1', 2;
        '2025-05-28 12-10-14_vmc_checked.mat', 'beringSea',   3;
        '2025-05-28 12-35-16_vmc_checked.mat', 'okhotskSea',  3;
        '2025-05-28 12-54-39_vmc_checked.mat', 'controlSea1', 3;
        '2025-05-28 13-14-15_vmc_checked.mat', 'beringSea',   4;
        '2025-05-28 13-40-20_vmc_checked.mat', 'okhotskSea',  4;
        '2025-05-28 14-00-01_vmc_checked.mat', 'controlSea1', 4;
        '2025-05-28 14-20-03_vmc_checked.mat', 'controlT1',   5;
        '2025-05-28 14-45-15_vmc_checked.mat', 'microT100',   5;
        '2025-05-28 15-05-29_vmc_checked.mat', 'microT10',    5;
        '2025-05-28 15-25-08_vmc_checked.mat', 'rotation90',  5;
        '2025-05-28 15-43-24_vmc_checked.mat', 'controlT1',   6;
        '2025-05-28 16-08-08_vmc_checked.mat', 'microT100',   6;
        '2025-05-28 16-27-27_vmc_checked.mat', 'microT10',    6;
        '2025-05-28 16-47-22_vmc_checked.mat', 'rotation90',  6;
        '2025-05-28 17-06-11_vmc_checked.mat', 'controlT1',   7;
        '2025-05-28 17-25-21_vmc_checked.mat', 'microT100',   7;
        '2025-05-28 17-44-35_vmc_checked.mat', 'microT10',    7;
        '2025-05-28 18-03-16_vmc_checked.mat', 'rotation90',  7
    };

    % --- 初期化 ---
    results = struct();
    data_table = table();

    % --- パラメータ ---
    fps = 30;
    hover_threshold = 0;
    num_tracks = 9;  % ファイル内のトラック数（=魚数）

    % --- 各ファイル処理 ---
    for i = 1:size(mat_files, 1)
        filename  = mat_files{i,1};
        condition = mat_files{i,2};
        cohort_id = mat_files{i,3};   % ★群ID（同じ個体群を示す）

        % 速度とヘディング
        [speed, heading] = magnetBehav(filename);

        % 時間ベクトル
        time = (1:size(speed, 2)) / fps;
        time_after_5min = time > (5 * 60);

        % 結果構造体の初期化
        if ~isfield(results, condition)
            results.(condition) = struct();
            results.(condition).sessions = [];
        end

        % セッション構造体
        session.file = filename;
        session.cohort_id = cohort_id;      % ★群IDを保持
        session.speed = speed;
        session.heading = heading;
        session.time = time;
        session.time_after_5min = time_after_5min;
        session.hover_indices = cell(1, num_tracks);
        session.hover_headings = NaN(1, num_tracks);
        session.fish_global_ids = (cohort_id-1)*num_tracks + (1:num_tracks); % ★このセッションの各トラックのグローバル個体ID

        for trk = 1:num_tracks
            % ★ 要件どおりの通し番号（1..9, 10..18, 19..27, ...）
            fish_global = (cohort_id - 1) * num_tracks + trk;

            hover_idx = find(speed(trk, :) < hover_threshold & time_after_5min);
            session.hover_indices{trk} = hover_idx;
            if ~isempty(hover_idx)
                session.hover_headings(trk) = circ_mean(heading(trk, hover_idx));
            end

            % CSV用テーブル（1サンプル=1行）
            tmp = table();
            tmp.Session   = repmat(string(filename), numel(time), 1);
            tmp.Condition = repmat(string(condition), numel(time), 1);
            tmp.Cohort    = repmat(cohort_id,         numel(time), 1);
            tmp.Fish      = repmat(fish_global,       numel(time), 1);   % ★ここを通し番号に
            tmp.Track     = repmat(trk,               numel(time), 1);   % ファイル内トラック(1..9)
            tmp.Time      = time(:);
            tmp.Speed     = speed(trk, :)';
            tmp.Heading   = heading(trk, :)';
            data_table = [data_table; tmp]; %#ok<AGROW>
        end

        % セッション追加
        results.(condition).sessions = [results.(condition).sessions; session]; %#ok<AGROW>
    end

    % CSV出力
    writetable(data_table, 'speed_heading_data2025.csv');

    results_stats = results;
end
