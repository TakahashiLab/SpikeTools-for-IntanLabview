function analyze_fish_data(fn)
debug=0;
    % Load the data
    data = readtable(fn);

    % Define constants
    fps = 30;
    time_limit = 20 * 60 * fps; % 20 minutes in frames

    % Get unique conditions and fish IDs
    conditions = unique(data.Condition);
    fish_ids = unique(data.Fish);

    % Initialize arrays to store averaged data
    averaged_data = table();

    % Loop through each condition and fish to filter and calculate mean values
    for i = 1:length(conditions)
        for j = 1:length(fish_ids)
            % Filter data for the current condition and fish
            condition_fish_data = data(strcmp(data.Condition, conditions{i}) & data.Fish == fish_ids(j), :);

            % Apply the time limit to speed and heading data
            if height(condition_fish_data) > time_limit
                condition_fish_data = condition_fish_data(time_limit+1:end, :);
            else
                continue;
            end

            % Calculate mean values
            mean_speed = mean(condition_fish_data.Speed);
            mean_acceleration = mean([0; diff(condition_fish_data.Speed)]);
            %mean_activity_frequency = mean(condition_fish_data.Speed > 0.1); % Example threshold for activity frequency
            mean_activity_frequency = calculateMomentaryActivityIndex(condition_fish_data.Speed);
            mean_heading = circ_mean(deg2rad(condition_fish_data.Heading));
            mean_heading_revolutions = mean(cumsum([0; diff(condition_fish_data.Heading)]));

            mean_values = table(mean_speed, mean_acceleration, mean_activity_frequency, mean_heading, mean_heading_revolutions, ...
                                fish_ids(j), conditions(i), ...
                                'VariableNames', {'Speed', 'Acceleration', 'ActivityFrequency', 'Heading', 'HeadingRevolutions', 'Fish', 'Condition'});

            % Append the mean values to the averaged data table
            averaged_data = [averaged_data; mean_values];
        end
 
    end
 
    % Define control conditions for different types of conditions
    control_sea = 'controlSea2';
    control_rot = 'controlRot1';
    control_t1 = 'controlT1';

    % Separate data into different experimental groups, excluding specified controls
    %sea_conditions = setdiff(conditions(contains(conditions, 'Sea') & ~contains(conditions, 'control')), {control_sea});
    sea_conditions = setdiff(conditions(contains(conditions, 'Sea') ), {control_sea});
    %rot_conditions = setdiff(conditions(contains(conditions, 'rotation') & ~contains(conditions, 'control')), {control_rot});
    rot_conditions = setdiff(conditions(contains(conditions, 'rotation') | contains(conditions, 'controlRot2')), {control_rot});
    %microt_conditions = setdiff(conditions(contains(conditions, 'microT') & ~contains(conditions, 'control')), {control_t1});
    microt_conditions = setdiff(conditions(contains(conditions, 'microT') | contains(conditions, 'controlT')), {control_t1});

    % Perform ANOVA and post-hoc tests for Sea conditions
    group_label = 'sea';
    perform_anova_and_posthoc(averaged_data, control_sea, sea_conditions, group_label);

    % Perform ANOVA and post-hoc tests for Rot conditions
    group_label = 'rotation';
    perform_anova_and_posthoc(averaged_data, control_rot, rot_conditions, group_label);

    % Perform ANOVA and post-hoc tests for microT conditions
    group_label = 'microT';
    perform_anova_and_posthoc(averaged_data, control_t1, microt_conditions, group_label);

return

function perform_anova_and_posthoc(data, control_label, condition_group, group_label)
    % Combine control data with experimental data for ANOVA
    anova_data = table();
    anova_conditions = {};
    anova_groups = {};
    anova_fish = [];

    % Filter out specified control and condition data
    control_data = data(strcmp(data.Condition, control_label), :);
    for i = 1:length(condition_group)
        condition = condition_group{i};
        condition_data = data(strcmp(data.Condition, condition), :);

        % Combine control and condition data
        combined_data = [control_data; condition_data];
        combined_conditions = [repmat({control_label}, size(control_data, 1), 1); repmat({condition}, size(condition_data, 1), 1)];
        combined_groups = repmat({group_label}, size(combined_data, 1), 1);
        combined_fish = [control_data.Fish; condition_data.Fish];

        % Append to ANOVA data
        anova_data = [anova_data; combined_data];
        anova_conditions = [anova_conditions; combined_conditions];
        anova_groups = [anova_groups; combined_groups];
        anova_fish = [anova_fish; combined_fish];
    end
    
    
      %  unique(anova_conditions)
    
    % Convert data for ANOVA
    anova_table = table(anova_conditions, anova_groups, anova_fish, ...
                        anova_data.Speed, anova_data.Acceleration, anova_data.ActivityFrequency, ...
                        'VariableNames', {'Condition', 'Group', 'Fish', 'Speed', 'Acceleration', 'ActivityFrequency'});

    % Create a repeated measures model
    rm = fitrm(anova_table, 'Speed-ActivityFrequency~Condition', 'WithinDesign', table([1 2 3]', 'VariableNames', {'Measurements'}));

    % Perform two-way repeated measures ANOVA
    ranova_results = ranova(rm);

    % Display results
    fprintf('***---***\n');
    disp(['Results for ', group_label, ' conditions:']);
    
    % Display ANOVA table
    disp(ranova_results);

    % Determine which variables are significantly different
    sig_measures = find(ranova_results.pValue < 0.05);
    sig_measure_names = {'Speed', 'Acceleration', 'ActivityFrequency'};
    if ~isempty(sig_measures)
        disp('Significant differences found in:');
        for i = 1:length(sig_measures)
            disp(sig_measure_names{sig_measures(i)});
        end
    else
        disp('No significant differences found.');
    end

    % Perform post-hoc tests if significant
    if any(ranova_results.pValue < 0.05)
        posthoc_results = multcompare(rm, 'Condition');
        sig_posthoc = posthoc_results(posthoc_results.pValue < 0.05, :);
        disp('Significant post-hoc comparisons:');
        for i = 1:height(sig_posthoc)
            disp([sig_posthoc.Condition_1{i}, ' vs ', sig_posthoc.Condition_2{i}, ': p-value = ', num2str(sig_posthoc.pValue(i))]);
        end
    end

   
    % Suppress detailed warnings
    warning('off', 'all');
    anova_conditions2=unique(anova_conditions);%%
      
     % Perform Watson-Williams test for Heading
    heading_data = anova_data.Heading;
    heading_conditions = grp2idx(anova_conditions2); % Convert conditions to numeric
    

    % Perform pairwise Watson-Williams tests for heading
    significant_pairs = {};
    unique_conditions = unique(heading_conditions);
    num_comparisons = 0;
    for i = 1:length(unique_conditions)
        for j = i+1:length(unique_conditions)
            cond1 = unique_conditions(i);
            cond2 = unique_conditions(j);
            data1 = heading_data(heading_conditions == cond1);
            data2 = heading_data(heading_conditions == cond2);
            
            try
                [p_pairwise,~] = circ_wwtest(data1, data2);
                num_comparisons = num_comparisons + 1;
                if p_pairwise < 0.05 / num_comparisons % Bonferroni correction
                    significant_pairs{end+1} = [anova_conditions2{cond1}, ' vs ', anova_conditions2{cond2}, ': p-value = ', num2str(p_pairwise)];
                end
            catch
                % Ignore errors
            end
        end
    end
    
    % Re-enable warnings
    warning('on', 'all');

    % Display results for Heading
    disp('Significant differences found in Heading between conditions:');
    for k = 1:length(significant_pairs)
        disp(significant_pairs{k});
    end

return

function activity_index = calculateMomentaryActivityIndex(speed)
% Initialize activity index array


% Calculate the standard deviation of speed and heading variance for each fish
speed_std = std(speed);  % Standard deviation across time points for each fish

% Define the threshold as 3 times the standard deviation
speed_threshold = 3 * speed_std;


% Calculate the activity index as a binary indicator if speed is above the threshold
is_active_speed = speed > speed_threshold;



% Combine both conditions for activity index (either condition being true)
activity_index = sum(is_active_speed)/length(speed);

return;

