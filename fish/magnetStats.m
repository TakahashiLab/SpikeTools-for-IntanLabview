function results_stats=magnetStats()
% Ensure Circular Statistics Toolbox is installed and add its path
% addpath('path_to_circular_statistics_toolbox');

% List of MAT files and their corresponding conditions
mat_files = {
    '2024-04-23-14-31.mat', 'controlSea1';
    '2024-04-23-15-19.mat', 'ohotsukuSea';
    '2024-04-23-16-08.mat', 'beringSea';
    '2024-04-23-16-56.mat', 'controlSea2';
    '2024-04-24-09-02.mat', 'controlRot';
    '2024-04-24-09-49.mat', 'rotation90';
    '2024-04-24-10-38.mat', 'rotation180';
    '2024-04-24-11-28.mat', 'rotation270';
    '2024-04-24-13-04.mat', 'controlT1';
    '2024-04-24-13-53.mat', 'microT200';
    '2024-04-24-14-41.mat', 'microT100';
    '2024-04-24-15-31.mat', 'microT150';
    '2024-04-24-16-21.mat', 'controlT2'
};

% Initialize results storage
results = struct();

% Parameters
fps = 30;
hover_threshold = 5;  % Example threshold for hovering speed in units per second

% Include the path to the directory containing the magnetBehav.m file
%addpath('/mnt/data');

% Process each MAT file
for i = 1:size(mat_files, 1)
    % Calculate speed and heading using magnetBehav function
    [speed, heading] = magnetBehav(mat_files{i, 1});

    % Time in seconds
    time = (1:size(speed, 2)) / fps;

    % Time after 5 minutes
    time_after_5min = time > 9 * 60;%9 min

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
    end
end

% Define conditions for comparison
control_conditions = {'control'};
sea_conditions = {'ohotsukuSea', 'beringSea'};
rotation_conditions = {'rotation90', 'rotation180', 'rotation270'};
magnetic_conditions = {'microT100', 'microT150', 'microT200'};

% Perform circular statistical tests and apply Bonferroni correction
results_stats = struct();

all_conditions = [sea_conditions, rotation_conditions, magnetic_conditions];
num_tests = length(all_conditions) * 9; % Total number of tests

for i = 1:length(all_conditions)
    condition = all_conditions{i};
    control_data = deg2rad(cell2mat(results.controlRot.hover_headings));
    test_data = deg2rad(cell2mat(results.(condition).hover_headings));
    if ~isempty(control_data) && ~isempty(test_data)
        [pval, table] = circular_test(control_data, test_data);
       
        % Apply Bonferroni correction
        pval
        pval_corrected = min(pval * num_tests, 1);
        results_stats.(condition).pval = pval;
        results_stats.(condition).table = table;
    else
        results_stats.(condition).pval = NaN;
    end
end

% Display the p-values for each condition
disp('Statistical analysis results (Bonferroni corrected p-values):');
disp(results_stats);
return;


%%%%
% Helper function for circular test
function [pval, table] = circular_test(control_headings, test_headings)
    % Perform Watson-Williams test
    [pval, table] = circ_wwtest(control_headings, test_headings);
return;

%%%
% Plotting the results
figure;
hold on;

conditions = fieldnames(results);
for i = 1:length(conditions)
    condition = conditions{i};
    headings = [];
    for fish = 1:9
        hover_indices = results.(condition).hover_indices{fish};
        if ~isempty(hover_indices)
            headings = [headings, results.(condition).heading(fish, hover_indices)];
        end
    end
    if ~isempty(headings)
   %     polarhistogram(deg2rad(headings), 30, 'DisplayName', condition);
    end
end

title('Heading Direction During Hovering for Different Conditions');
legend;
hold off;

return;
