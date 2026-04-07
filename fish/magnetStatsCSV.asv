function results_stats = magnetStatsCSV()
    % Ensure Circular Statistics Toolbox is installed and add its path
    % addpath('path_to_circular_statistics_toolbox');

    % List of MAT files and their corresponding conditions
    % mat_files = {
    %     '2024-04-23-14-31.mat', 'controlSea1';
    %     '2024-04-23-15-19.mat', 'okhotskSea';
    %     '2024-04-23-16-08.mat', 'beringSea';
    %     '2024-04-23-16-56.mat', 'controlSea2';
    %     '2024-04-24-09-02.mat', 'controlRot1';
    %     '2024-04-24-09-49.mat', 'rotation90';
    %     '2024-04-24-10-38.mat', 'rotation180';
    %     '2024-04-24-11-28.mat', 'rotation270';
    %     '2024-04-24-13-04.mat', 'controlRot2';%duplicate
    %     '2024-04-24-13-04.mat', 'controlT1';
    %     '2024-04-24-13-53.mat', 'microT200';
    %     '2024-04-24-14-41.mat', 'microT100';
    %     '2024-04-24-15-31.mat', 'microT150';
    %     '2024-04-24-16-21.mat', 'controlT2'
    % };

    mat_files = {
        '2025-05-28 10-01-55_vmc_checked.mat', 'beringSea',1;
        '2025-05-28 10-26-38_vmc_checked.mat', 'okhotskSea',1;
        '2025-05-28 10-47-04_vmc_checked.mat', 'controlSea1',1;
        '2025-05-28 11-05-32_vmc_checked.mat', 'beringSea',2;
        '2025-05-28 11-30-09_vmc_checked.mat', 'okhotskSea',2;
        '2025-05-28 11-49-48_vmc_checked.mat', 'controlSea1',2;
        '2025-05-28 12-10-14_vmc_checked.mat', 'beringSea',3;
        '2025-05-28 12-35-16_vmc_checked.mat', 'okhotskSea',3;
        '2025-05-28 12-54-39_vmc_checked.mat', 'controlSea1',3;
        '2025-05-28 13-14-15_vmc_checked.mat', 'beringSea',4;
        '2025-05-28 13-40-20_vmc_checked.mat', 'okhotskSea',4;
        '2025-05-28 14-00-01_vmc_checked.mat', 'controlSea1',4;
        '2025-05-28 14-20-03_vmc_checked.mat', 'controlT1',5;
        '2025-05-28 14-45-15_vmc_checked.mat', 'microT100',5;
        '2025-05-28 15-05-29_vmc_checked.mat', 'microT10',5;
        '2025-05-28 15-25-08_vmc_checked.mat', 'rotation90',5;
        '2025-05-28 15-43-24_vmc_checked.mat', 'controlT1',6;
        '2025-05-28 16-08-08_vmc_checked.mat', 'microT100',6;
        '2025-05-28 16-27-27_vmc_checked.mat', 'microT10',6;
        '2025-05-28 16-47-22_vmc_checked.mat', 'rotation90',6;
        '2025-05-28 17-06-11_vmc_checked.mat', 'controlT1',7;
        '2025-05-28 17-25-21_vmc_checked.mat', 'microT100',7;
        '2025-05-28 17-44-35_vmc_checked.mat', 'microT10',7;
        '2025-05-28 18-03-16_vmc_checked.mat', 'rotation90',7
    };

    % Initialize results storage
    results = struct();

    % Initialize table for storing speed and heading data
    data_table = table();

    % Parameters
    fps = 30;
    hover_threshold = 0;  %5 Example threshold for hovering speed in units per second

    % Include the path to the directory containing the magnetBehav.m file
    % addpath('/mnt/data');

    % Process each MAT file
    for i = 1:size(mat_files, 1)
        % Calculate speed and heading using magnetBehav function
        [speed, heading] = magnetBehav(mat_files{i, 1});

        % Time in seconds
        time = (1:size(speed, 2)) / fps;

        % Time after 5 minutes
        time_after_5min = time >  5* 60;  % 5 minutes

        % Store results for each condition
        results.(mat_files{i, 2}).speed = speed;
        results.(mat_files{i, 2}).heading = heading;
        results.(mat_files{i, 2}).time = time;
        results.(mat_files{i, 2}).time_after_5min = time_after_5min;

        % Find indices where speed is below threshold after 5 minutes
        for fish = 1:9
            hover_indices = find(speed(fish, :) < hover_threshold & time_after_5min);
            results.(mat_files{i, 2}).hover_indices{fish} = hover_indices;
            results.(mat_files{i, 2}).hover_headings{fish} = circ_mean(heading(fish, hover_indices));

            % Create a temporary table for the current fish and condition
            temp_table = table();
            % temp_table.Condition = repmat({mat_files{i, 2}}, length(hover_indices), 1);
            % temp_table.Fish = repmat(fish, length(hover_indices), 1);
            % temp_table.Speed = speed(fish, hover_indices)';
            % temp_table.Heading = heading(fish, hover_indices)';
            temp_table.Condition = repmat({mat_files{i, 2}}, length(speed), 1);
            temp_table.Fish = repmat(fish, length(speed), 1);
            temp_table.Speed = speed(fish, :)';
            temp_table.Heading = heading(fish, :)';

            % Append to the main data table
            data_table = [data_table; temp_table];
        end
    end

    % Save the data table for later use
    writetable(data_table, 'speed_heading_data.csv');

    % Return the results struct
    results_stats = results;
end
