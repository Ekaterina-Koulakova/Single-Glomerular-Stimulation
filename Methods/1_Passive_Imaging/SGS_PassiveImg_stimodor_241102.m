clear all
close all
clc

set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultTextFontSize', 10);
set(0,'defaultAxesFontSize',10);
set(0,'DefaultFigureColormap',jet)
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultAxesTitleFontWeight', 'normal')
set(0,'DefaultTextInterpreter','none')

%% FONT SIZES

axis_size = 14;
title_size = 20;
legend_size = 12;

%% DEFINE FILE PATHS
% 
% FilePaths = {'mouse0773/240806_Glom_StimOdorMapping/aligned/240806_0773_Glom_StimOdorMapping_Z0_1X_S_v73.mat'};
% ROIindex = [10];
% 
% stim_y_lim = [-0.5, 1] % y-limits of produced plots
% odor_y_lim = [-0.2, 1]
% stim_odor_y_lim = [-0.2, 1]
% stim_time_avg_y_lim = [-0.4, 1]

% FilePaths = {'mouse0773/240806_MC_StimOdorMapping/aligned/240806_MC_StimOdorMapping_Z186_3X_S_v73.mat'};
% ROIindex = [3,7,9];
% 
% stim_y_lim = [-1, 1] % y-limits of produced plots
% odor_y_lim = [-1.2, 1]
% stim_odor_y_lim = [-1.2, 1]
% stim_time_avg_y_lim = [-1, 1]

FilePaths = {'mouse0781/240806_Glom_StimOdorMapping/aligned/240806_0781_Glom_StimOdorMapping_Z0_1X_S_v73.mat'};
ROIindex = [15]; 

stim_y_lim = [-0.5, 1] % y-limits of produced plots
odor_y_lim = [-0.5, 1.2]
stim_odor_y_lim = [-0.5, 1.2]
stim_time_avg_y_lim = [-0.2, 1]
% 
% FilePaths = {'mouse0781/240806_MC_StimOdorMapping/aligned/240806_0781_MC_StimOdorMapping_Z105_3X_S_v73.mat'};
% ROIindex = [7,9,18,25];
% 
% stim_y_lim = [-1, 1] % y-limits of produced plots
% odor_y_lim = [-1, 1]
% stim_odor_y_lim = [-1.5, 1]
% stim_time_avg_y_lim = [-0.8, 0.8]

%Extract the parent directory for mouse0773
commonDirectory = fileparts(fileparts(fileparts(FilePaths{1})));

% Create the common "plots" directory
plot_dir = fullfile(commonDirectory, 'plots');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

for FilePathIndex = 1:numel(FilePaths)
    % Determine the title based on the FilePath
    if contains(FilePaths{FilePathIndex}, 'Glom')
        layer_name = 'Glomerular Layer';
        layer_acr = 'Glom';
    elseif contains(FilePaths{FilePathIndex}, 'MC')
        layer_name = 'Mitral Cell Layer';
        layer_acr = 'MC';
    end

end

%% DEFINE SCALAR VALUES

% Load the data from the first file in FilePaths
data = load(FilePaths{1}, 'Session');

n_ROI = size(data.Session.OdorResponse{1,1}, 2);
n_trials = size(data.Session.OdorResponse{1}, 3);
n_conds = numel(data.Session.OdorResponse);

% TIME VARIABLES
pre_inh = data.Session.Infos.pre_inh;
post_inh = data.Session.Infos.post_inh;
n_frames = pre_inh + post_inh;
n_fps = data.Session.Infos.fps;
tt = linspace(-(pre_inh/n_fps), (post_inh/n_fps), n_frames);

%% SCALAR VALUES USED TO AVERAGE TIME INCREMENTS
% x-axis = latency values
% legend:
% EMPTY
% STIM-ONLY
% ODOR-ONLY (low and high)
% ODOR+STIM (low and high)

stim_time = [10, 30, 60, 120, 180, 240];

% stim_time to frame_time
stim_frame = (stim_time * 30 / 1000);
stim_frame = pre_inh + round(stim_frame);
stim_frame_end = stim_frame + 15;
stim_frame_end = round(stim_frame_end);

%% ODOR RESPONSE

ODOR_RESPONSE = data.Session.OdorResponse;


%% LINEAR INTERPOLATION TO GUISE PMT GATING

% Loop through each cell in ODOR_RESPONSE
for cellIdx = 1:numel(ODOR_RESPONSE)
    % Extract the current 180x27x10 array from the cell
    datta = ODOR_RESPONSE{cellIdx};
    
    % Create a copy of the original data for comparison  
    original_datta = datta;
    
    % Loop through each column in the 180x27x10 array
    for colIdx = 1:size(datta, 2)
        for sliceIdx = 1:size(datta, 3) % Loop through each "slice" (3rd dimension)

            local_avg = mean(datta(40:80, colIdx, sliceIdx));
            
            % Define the range of rows to search
            % 10 ms latency
            if cellIdx == 2 || cellIdx == 9 || cellIdx == 16
                range = 49:63;
            % 30 ms latency
            elseif cellIdx == 3 || cellIdx == 10 || cellIdx == 17
                range = 60:65;
            % 60 ms latency
            elseif cellIdx == 4 || cellIdx == 11 || cellIdx == 18
                range = 60:65;
            % 120 ms latency
            elseif cellIdx == 5 || cellIdx == 12 || cellIdx == 19
                range = 58:67;
            % 180 ms latency
            elseif cellIdx == 6 || cellIdx == 13 || cellIdx == 20
                range = 61:70;
            % 240 ms latency
            elseif cellIdx == 7 || cellIdx == 14 || cellIdx == 21
                range = 61:70;
            else
                range = []; % No range if cellIdx doesn't match any case
            end
            
            if ~isempty(range)
                % Extract data within the specified range
                columnData = datta(range, colIdx, sliceIdx);
                
                % Condition: Proceed only if the criteria are met
                if (local_avg < 0 && any(columnData < -0.04)) || (local_avg > 0 && any(columnData < -0.02))
                    
                    % Step 1: Calculate the 3-point running average within the range
                    runningAvg = movmean(columnData, 3, 'Endpoints', 'discard');

                    % Step 2: Find the indices of the 3-point segment with the minimum running average
                    [~, minStartIdx] = min(runningAvg); % Starting index of the minimum 3-point segment
                    minIdxs = minStartIdx:minStartIdx+2; % Indices of the 3-row segment

                    % Adjust indices relative to the entire data array
                    minRowIdxs = range(minIdxs);

                    % Step 3: Linearize each of the identified 3 rows using interpolation
                    if ~isempty(minRowIdxs)
                        % Get the indices before and after the segment for interpolation
                        startIdx = minRowIdxs(1) - 1; % row 60 in the example
                        endIdx = minRowIdxs(end) + 1; % row 64 in the example

                        % Ensure bounds to avoid errors at the edges
                        if startIdx > 0 && endIdx <= size(datta, 1)
                            % Values at the start and end indices
                            startValue = datta(startIdx, colIdx, sliceIdx);
                            endValue = datta(endIdx, colIdx, sliceIdx);

                            % Linearly interpolate values between startValue and endValue
                            interpolatedValues = linspace(startValue, endValue, numel(minRowIdxs) + 2);

                            % Replace values in minRowIdxs with interpolated values (ignoring the first and last)
                            datta(minRowIdxs, colIdx, sliceIdx) = interpolatedValues(2:end-1);
                        end
                    end
                end
            end
            
            % VISUAL DEBUGGING AID
%             if (cellIdx == 21) && (colIdx == 15)
%                 figure(100)
%                 clf
%                 plot(range, original_datta(range, colIdx, sliceIdx), 'r.-', 'DisplayName', 'Data Before Correction'), hold on
%                 plot(1:180, datta(:, colIdx, sliceIdx), 'bx-', 'DisplayName', 'Data After PMT Correction')
%                 legend
%                 title(sprintf('Debugging: Cell %d, Column %d, Slice %d', cellIdx, colIdx, sliceIdx))
%                 pause
%             end
        end
    end

    % Store the modified data back in the cell array
    ODOR_RESPONSE{cellIdx} = datta;
end


%% LINEAR INTERPOLATION TO GUISE STIM ARTIFACT

% Loop through each cell in ODOR_RESPONSE
for cellIdx = 1:numel(ODOR_RESPONSE)
    % Extract the current 180x27x10 array from the cell
    datta = ODOR_RESPONSE{cellIdx};
    
    % Create a copy of the original data for comparison  
    original_datta = datta;
    
    % Loop through each column in the 180x27x10 array
    for colIdx = 1:size(datta, 2)
        for sliceIdx = 1:size(datta, 3) % Loop through each "slice" (3rd dimension)

            local_avg = mean(datta(40:80, colIdx, sliceIdx));
            
            % Define the range of rows to search
            % 10 ms latency
            if cellIdx == 2 || cellIdx == 9 || cellIdx == 16
                range = 49:75;
            % 30 ms latency
            elseif cellIdx == 3 || cellIdx == 10 || cellIdx == 17
                range = 60:75;
            % 60 ms latency
            elseif cellIdx == 4 || cellIdx == 11 || cellIdx == 18
                range = 60:75;
            % 120 ms latency
            elseif cellIdx == 5 || cellIdx == 12 || cellIdx == 19
                range = 58:75;
            % 180 ms latency
            elseif cellIdx == 6 || cellIdx == 13 || cellIdx == 20
                range = 61:75;
            % 240 ms latency
            elseif cellIdx == 7 || cellIdx == 14 || cellIdx == 21
                range = 61:75;
            else
                range = []; % No range if cellIdx doesn't match any case
            end
            
            if ~isempty(range)
                % Extract data within the specified range
                columnData = datta(range, colIdx, sliceIdx);
                
                % Condition: Proceed only if the criteria are met
                if (local_avg < 0 && any(columnData > 0.04)) || (local_avg > 0 && any(columnData > 0.06))
                    
                    % Step 1: Calculate the 3-point running average within the range
                    runningAvg = movmean(columnData, 5, 'Endpoints', 'discard');

                    % Step 2: Find the indices of the 3-point segment with the minimum running average
                    [~, maxStartIdx] = max(runningAvg); % Starting index of the minimum 3-point segment
                    maxIdxs = maxStartIdx:maxStartIdx+4; % Indices of the 3-row segment

                    % Adjust indices relative to the entire data array
                    maxRowIdxs = range(maxIdxs);

                    % Step 3: Linearize each of the identified 3 rows using interpolation
                    if ~isempty(maxRowIdxs)
                        % Get the indices before and after the segment for interpolation
                        startIdx = maxRowIdxs(1) - 1; % row 60 in the example
                        endIdx = maxRowIdxs(end) + 1; % row 64 in the example

                        % Ensure bounds to avoid errors at the edges
                        if startIdx > 0 && endIdx <= size(datta, 1)
                            % Values at the start and end indices
                            startValue = datta(startIdx, colIdx, sliceIdx);
                            endValue = datta(endIdx, colIdx, sliceIdx);

                            % Linearly interpolate values between startValue and endValue
                            interpolatedValues = linspace(startValue, endValue, numel(maxRowIdxs) + 2);

                            % Replace values in minRowIdxs with interpolated values (ignoring the first and last)
                            datta(maxRowIdxs, colIdx, sliceIdx) = interpolatedValues(2:end-1);
                        end
                    end
                end
            end
            
            % VISUAL DEBUGGING AID
            % VISUAL DEBUGGING AID
%             if (cellIdx == 3) && (colIdx == 15)
%                 figure(100)
%                 clf
%                 plot(range, original_datta(range, colIdx, sliceIdx), 'r.-', 'DisplayName', 'Data Before Correction'), hold on
%                 plot(1:180, datta(:, colIdx, sliceIdx), 'bx-', 'DisplayName', 'Data After Stim Correction')
%                 legend
%                 title(sprintf('Debugging: Cell %d, Column %d, Slice %d', cellIdx, colIdx, sliceIdx))
%                 pause
%             end
        end
    end

    % Store the modified data back in the cell array
    ODOR_RESPONSE{cellIdx} = datta;
end


%% Average ODOR_RESPONSE accross trials (3rd dimention in double)

ODOR_RESPONSE_avg = cell(size(ODOR_RESPONSE));
ODOR_RESPONSE_avg = cellfun(@(x) mean(x, 3), ODOR_RESPONSE, 'UniformOutput', false);


%% NORMALIZE ODOR RESPONSE

ODOR_RESPONSE_gmax = -Inf;

% Loop through the averaged cells to find the global maximum
for idx_ROI = 1:length(ROIindex)
    r = ROIindex(idx_ROI);
    for i = 1:numel(ODOR_RESPONSE_avg)
        currentMax = max(ODOR_RESPONSE_avg{i}(45:90,r)); % Find the max in the current matrix
        if currentMax > ODOR_RESPONSE_gmax
            ODOR_RESPONSE_gmax = currentMax;
        end
    end
end

% Normalize ODOR_RESPONSE_avg
ODOR_RESPONSE_norm = cell(size(ODOR_RESPONSE_avg));
for i = 1:numel(ODOR_RESPONSE_avg)
    ODOR_RESPONSE_norm{i} = ODOR_RESPONSE_avg{i} / ODOR_RESPONSE_gmax;
end

%% CHECK WHERE GLOBAL MAX IS LOCATED

% Loop through each condition
for cond_idx = 1:n_conds
    % Get the current 180x27 matrix
    response_matrix = ODOR_RESPONSE_avg{cond_idx};
    
    % Find the location of the value equal to ODOR_RESPONSE_gmax
    [row, col] = find(response_matrix == ODOR_RESPONSE_gmax);
    
    % Check if the value is found
    if ~isempty(row)
        % Output the location (frame, ROI, condition)
        fprintf('Found ODOR_RESPONSE_gmax at:\n');
        for idx = 1:length(row)
            condition = data.Session.UniqueConds{4}
            fprintf('The maximum value occurs at: \n');
            fprintf('Condition: %d, Frame: %d, ROI: %d\n', cond_idx, row(idx), col(idx));
            fprintf('%s', condition);
        end
    end
end


%% ORGANIZING UNIQUE CONDITIONS
% 
% Number of conditions
n_conds = length(data.Session.UniqueConds);

% Initialize cell arrays for odors, concentrations, and latencies
cond_odors = cell(n_conds, 1);
cond_concentration = cell(n_conds, 1);
cond_latencies = cell(n_conds, 1);

% Loop through each element in UniqueConds
for i = 1:n_conds
    cond = data.Session.UniqueConds{i};
    
    % Extract odor (string before the ':')
    % Extract substring prior to ':'
    before_colon = strtok(cond, ':');
    
    % Split the substring by '/'
    odor_parts = strsplit(before_colon, '/');
    
    % Initialize variable to store valid odor
    valid_odors = {};
    
    % Check each part and store those longer than 2 characters
    for j = 1:length(odor_parts)
        if length(odor_parts{j}) > 2
            valid_odors{end + 1} = odor_parts{j}; % Add valid part to list
        end
    end
    
    % Store the first valid odor found, or leave empty if none
    if ~isempty(valid_odors)
        cond_odors{i} = valid_odors{1}; % Store the first valid odor
    else
        cond_odors{i} = ''; % If no valid odor, leave empty
    end
    
    % Extract concentration
    % Extract substring after ':'
    after_colon = strtrim(cond(length(before_colon) + 1:end)); % Get the part after ':'
    
    % Split the substring by '/'
    concentration_parts = strsplit(after_colon, '/');
    
    % Initialize variable to store valid concentrations
    valid_concentrations = {};
    
    % Check each part and store those that are numeric and not '0.000'
    for j = 1:length(concentration_parts)
        is_numeric = ~isempty(regexp(concentration_parts{j}, '^\d+(\.\d+)?$', 'once'));
        is_not_zero = ~strcmp(concentration_parts{j}, '0.000');
        
        if is_numeric && is_not_zero
            % Convert to number, divide by 10, and store as string
            conc_value = str2double(concentration_parts{j}) / 10;
            valid_concentrations{end + 1} = num2str(conc_value, '%.10g'); % Store as string without trailing zeros
        end
    end
    
    % Store the first valid concentration found, or leave empty if none
    if ~isempty(valid_concentrations)
        cond_concentration{i} = valid_concentrations{1}; % Store the first valid concentration
    end
    
    % Extract latency
    latency_match = regexp(cond, '-SL\s*(\d+)', 'tokens');
    if ~isempty(latency_match)
        cond_latencies{i} = latency_match{1}{1};
    else
        cond_latencies{i} = ''; % If no match, leave empty
    end
end

% Remove empty cells (if necessary) before finding unique values
cond_concentration_ne = cond_concentration(~cellfun('isempty', cond_concentration));
cond_latencies_ne = cond_latencies(~cellfun('isempty', cond_latencies));

unique_odors = unique(cond_odors, 'stable');
unique_concentration = unique(cond_concentration_ne, 'stable');
unique_latencies = num2cell(unique(cond_latencies_ne, 'stable'));

n_latencies = numel(unique_latencies);
n_odors = numel(unique_odors)-1;
n_conc = numel(unique_concentration);

%% ORGANIZING ODOR RESPONSE BASED ON UNIQUE CONDITIONS

% If-then statements mock-up:
% ODOR_ONLY[1]: cond_odors = unique_odors[1], latencies = [empty]
% ODOR_ONLY[2]: cond_odors = unique_odors[2], latencies = [empty]
% ODOR_STIM: cond_odors = unique_odors[1], cond_latencies = [10,30,60,120,180,240]
% STIM_ONLY: cond_odors = 'empty', cond_latencies = [10,30,60,120,180,240]
% EMPTY: cond_odors = 'empty', cond_latencies = [empty]

% Initialize the output cell arrays
STIM_ONLY = cell(1, n_latencies);  % 6 
ODOR_ONLY = cell(n_conc, n_odors);  % 2
ODOR_STIM = cell(n_conc, n_odors);  % 12
for i = 1:n_odors   
    ODOR_STIM{i} = cell(1, n_latencies);
end
EMPTY = []; % 1

% ? unique_latencies = unique_latencies'

unique_latencies = {'10', '30', '60', '120', '180', '240'};

% Iterate through each element in ODOR_RESPONSE_norm
for i = 1:n_conds
    % Get the index of cond_latencies in unique_latencies
    latency_idx = find(strcmp(unique_latencies, cond_latencies{i}));
    
    % Get the index of cond_concentration in unique_concentration
    concentration_idx = find(strcmp(unique_concentration, cond_concentration{i}));
    
    % Check for EMPTY condition
    if strcmp(cond_odors{i}, 'empty') && isempty(cond_latencies{i})
        EMPTY = ODOR_RESPONSE_norm{i};
        
    % Check for STIM_ONLY condition
    elseif strcmp(cond_odors{i}, 'empty') && ~isempty(latency_idx)
        STIM_ONLY{latency_idx} = ODOR_RESPONSE_norm{i};
    
    % Check for ODOR_ONLY condition
    elseif ~strcmp(cond_odors{i}, 'empty') && isempty(cond_latencies{i})
        odor_idx = find(strcmp(unique_odors, cond_odors{i}));
        if ~isempty(odor_idx) && concentration_idx > 0
            ODOR_ONLY{odor_idx, concentration_idx} = ODOR_RESPONSE_norm{i};
        end
    
    % Check for ODOR_STIM condition
    elseif ~strcmp(cond_odors{i}, 'empty') && ~isempty(latency_idx)
        odor_idx = find(strcmp(unique_odors, cond_odors{i}));
        if ~isempty(odor_idx) && concentration_idx > 0
            ODOR_STIM{odor_idx, concentration_idx}{latency_idx} = ODOR_RESPONSE_norm{i};
        end
    end
end

%% ___________________________ STIM - ONLY ________________________________
%% EMPTY [EMPTY_frame_avg]

% Extract the columns from EMPTY corresponding to ROIindex
selected_columns = EMPTY(:, ROIindex);

% Average the selected columns along the rows
EMPTY_roi = mean(selected_columns, 2); % Average over rows, resulting in a 180x1 column vector

EMPTY_frame_avg = zeros(1, numel(stim_frame)); % Preallocate for 5 averages

% Loop through each pair of indices to compute the average
for idx = 1:numel(stim_frame)
    % Get the start and end indices for the current range
    start_idx = stim_frame(idx); % Convert to 1-based index
    end_idx = stim_frame_end(idx);
    
    % Average the values in EMPTY_roi from start_idx to end_idx
    EMPTY_frame_avg(idx) = mean(EMPTY_roi(start_idx:end_idx));
end

%% STIM-ONLY PLOTS

for ROI_idx = ROIindex %n_ROI

    % Create a new figure for each ROI
    figure_handle = figure('Name', ['ROI ', num2str(ROI_idx)]);

    % Set the figure dimensions: [left, bottom, width, height]
    set(figure_handle, 'Position', [10, 100, 600, 300]);

    tiledlayout(1, 1, 'Padding', 'compact'); % Single tile layout
        
    hold on;

    % Loop through each latency and plot the current ROI column of each matrix in OdorResponse_mean
    for latency_idx = 1:n_latencies
        
        currentLatency = STIM_ONLY{1,latency_idx}(:, ROI_idx);           
        
        legend_label = sprintf('%s ms', unique_latencies{latency_idx});
        
        % Set the alpha value based on the latency index
        alpha_value = 1 - (latency_idx - 1) * 0.15;
        
        % Create a plot with transparency
        h = line(tt, currentLatency, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', 'r');
        h.Color(4) = alpha_value;  % Set the alpha value

        hold on;

        xlim([-0.5, 1]);
        ylim(stim_y_lim); % Adjusted to [0, 1] since the data is normalized

        % Dashed vertical line at x = 0
        line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');

        line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

        % Add labels, title, etc. as needed
        xlabel('Time [s]','FontSize', axis_size);
        ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}','Interpreter','latex','FontName','Arial','FontSize', axis_size);

    end
    
    % Plot the no_stim data for the current ROI
    no_stim_ROI = EMPTY(:, ROI_idx); 
    
    % Plot the no_stim data
    plot(tt, no_stim_ROI, 'LineWidth', 1.5, 'DisplayName', 'no stim', 'Color', 'k');
        
    % Add legend
    legend_handle = legend('Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none', 'FontSize', legend_size);
    
    % Set the title for the entire figure
    sgtitle(['Single Glomerular Stimulation | ', layer_acr, ' ', num2str(ROI_idx)], 'FontSize', title_size);
    
    saveas(gcf, fullfile(plot_dir, sprintf('STIM-ONLY_%s_%d_plot.png', layer_acr, ROI_idx)));

end


%% STIM-ONLY [STIM_ONLY_avg]

for idx_latencies = 1:n_latencies
    selected_columns = STIM_ONLY{1,idx_latencies}(:, ROIindex);
    STIM_ONLY_roi_avg{1,idx_latencies} = mean(selected_columns, 2); 

    for idx = 1:numel(stim_frame)
    % Get the start and end indices for the current range
    start_idx = stim_frame(idx); % Convert to 1-based index
    end_idx = stim_frame_end(idx);
    
    % Average the values in EMPTY_roi from start_idx to end_idx
    STIM_ONLY_avg(idx) = mean(STIM_ONLY_roi_avg{1,idx_latencies}(start_idx:end_idx));
    end
end

%% STIM-ONLY PLOTS - AVERAGES ROI's
% Here, we will plot STIM_ONLY_roi_avg
% This plot averages the fluorescent curves for all of our ROI's for each
% latency
% In other words, all MC curves will be plotted together

% Create a new figure for each ROI
figure_handle = figure('Name', ['STIM-ONLY_', layer_acr, '_avg_roi']);

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 300]);

tiledlayout(1, 1, 'Padding', 'compact'); % Single tile layout

hold on;

% Loop through each latency and plot the current ROI column of each matrix in OdorResponse_mean
for latency_idx = 1:n_latencies

    currentLatency = STIM_ONLY_roi_avg{1,latency_idx};           

    legend_label = sprintf('%s ms', unique_latencies{latency_idx});

    % Set the alpha value based on the latency index
    alpha_value = 1 - (latency_idx - 1) * 0.15;

    % Create a plot with transparency
    h = line(tt, currentLatency, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', 'r');
    h.Color(4) = alpha_value;  % Set the alpha value

    hold on;

    xlim([-0.5, 1]);
    ylim(stim_y_lim); % Adjusted to [0, 1] since the data is normalized

    % Dashed vertical line at x = 0
    line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');

    line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

    % Add labels, title, etc. as needed

    if strcmp(layer_acr, 'Glom')
        % no x-label for Glom layer
    elseif strcmp(layer_acr, 'MC')
        xlabel('Time [s]','FontSize', axis_size);  
    end
    ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}','Interpreter','latex','FontName','Arial','FontSize', axis_size);

end

% Plot the no_stim data for the current ROI
no_stim_ROI = EMPTY_roi; 

% Plot the no_stim data
plot(tt, no_stim_ROI, 'LineWidth', 1.5, 'DisplayName', 'no stim', 'Color', 'k');

% Add legend
legend_handle = legend('Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none', 'FontSize', legend_size);


if strcmp(layer_acr, 'Glom')
    sgtitle('Selected Glomerulus');
elseif strcmp(layer_acr, 'MC')
    sgtitle('Averaged Daughter Mitral Cells');
end
% Set the title for the entire figure


saveas(gcf, fullfile(plot_dir, sprintf('%s_STIM-ONLY_avg_plot.png', layer_acr)));



%% ___________________________ ODOR - ONLY ________________________________
%% ODOR-ONLY Fluorescence Cross-Trace Plot
% (add low and high concentrations)

tt = tt(:); % turn tt into a column vector

% Prep no-odor
EMPTY_ROIindex = mean(EMPTY(:, ROIindex), 2);

% Define colors with alpha values
high_concentration_color = [0, 0, 1, 1.0];  % Blue color for higher concentration (RGBA)
low_concentration_color = [0, 0.5, 1, 0.5];   % Same blue for lower concentration (RGBA)

% Initialize legend entries
legend_entries = {};

% Loop through the chosen odor indices
for idx_odor = 1:n_odors
    
    % Get current odor name
    chosen_odor = unique_odors{idx_odor};

    % Create a new figure for each chosen odor index
    figure_handle = figure('Name', chosen_odor);
    % Set the figure dimensions: [left, bottom, width, height]
    set(figure_handle, 'Position', [10, 50, 600, 300]);

    hold on;
    
    % Dummy handles for the legend entries
    hLow = plot(nan, nan, 'Color', low_concentration_color(1:3), 'LineWidth', 1.5); % Dummy for low concentration
    hHigh = plot(nan, nan, 'Color', high_concentration_color(1:3), 'LineWidth', 1.5); % Dummy for high concentration
    hNoOdor = plot(nan, nan, 'Color', 'k', 'LineWidth', 1.5); % Dummy for no-odor

    for idx_conc = 1:n_conc

        odor_matrix = ODOR_ONLY{idx_odor, idx_conc};

        % Loop through each ROI index and plot the corresponding columns
        for j = 1:numel(ROIindex)
            roi_index = ROIindex(j);
    
            % Get the column corresponding to the current ROI
            ROI_odor_column = odor_matrix(:, roi_index);
    
            % Set the color based on the concentration index
            if idx_conc == 1 % Assuming this is the lower concentration
                color = low_concentration_color; % Set to lower concentration color
            else
                color = high_concentration_color; % Set to higher concentration color
            end
    
            % Plot the column for the current ROI
            hLine = plot(tt, ROI_odor_column, 'LineWidth', 1.5, 'Color', color(1:3)); % Plot without alpha first
            
            % Set the color with alpha using the Color property
            set(hLine, 'Color', color(1:3)); % Set RGB color
            alpha_value = color(4); % Set alpha value
            set(hLine, 'Color', [color(1:3) alpha_value]); % Set RGBA
        end
    end

    % Plot the no-odor response
    plot(tt, EMPTY_ROIindex, 'LineWidth', 1.5, 'DisplayName', 'no-odor', 'Color', 'k');

    xlim([-0.5, 1.5]);
    ylim(odor_y_lim);

    % Dashed vertical line at x = 0
    line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
    % Horizontal line at y = 0
    line(xlim, [0, 0], 'Color', [0, 0, 0, 0.2], 'LineStyle', '-', 'HandleVisibility', 'off');

    % Add labels, title, etc. as needed
    if strcmp(layer_acr, 'Glom')
        % no x-label for Glom layer
    elseif strcmp(layer_acr, 'MC')
        xlabel('Time [s]', 'FontSize', axis_size);  
    end
    
    ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', axis_size);

    % Add legend with the defined entries
    legend([hLow, hHigh, hNoOdor], {'low conc.', 'high conc.', 'no-odor'}, 'Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none');
    
    hold off;

    if strcmp(layer_acr, 'Glom')
        sgtitle('Selected Glomerulus');
    elseif strcmp(layer_acr, 'MC')
        sgtitle('Daughter Mitral Cells');
    end

    % Save the figure
    saveas(gcf, fullfile(plot_dir, sprintf('%s_ODOR-ONLY_plot.png', layer_acr)));

end

%% ODOR-ONLY Fluorescence Cross-Trace Plot (average ROIs)
% (add low and high concentrations)

tt = tt(:); % turn tt into a column vector

% Prep no-odor
EMPTY_ROIindex = mean(EMPTY(:, ROIindex), 2);

% Define colors with alpha values
high_concentration_color = [0, 0, 1, 1.0];  % Blue color for higher concentration (RGBA)
low_concentration_color = [0, 0.5, 1, 0.5]; % Lighter blue for lower concentration (RGBA)

% Loop through the chosen odor indices
for idx_odor = 1:n_odors
    
    % Get current odor name
    chosen_odor = unique_odors{idx_odor};

    % Create a new figure for each chosen odor index
    figure_handle = figure('Name', chosen_odor);
    % Set the figure dimensions: [left, bottom, width, height]
    set(figure_handle, 'Position', [10, 50, 600, 300]);

    hold on;

    % Initialize a list to store handles for the plotted lines
    plot_handles = [];
    
    for idx_conc = 1:n_conc

        odor_matrix = ODOR_ONLY{idx_odor, idx_conc};
        odor_matrix_avg_roi = mean(odor_matrix(:, ROIindex), 2);
    
        % Set the color based on the concentration index
        if idx_conc == 1 % Assuming this is the lower concentration
            color = low_concentration_color; % Set to lower concentration color
            display_name = sprintf('low %s conc.', chosen_odor);
        else
            color = high_concentration_color; % Set to higher concentration color
            display_name = sprintf('high %s conc.', chosen_odor);
        end

        % Plot the current concentration trace
        hLine = plot(tt, odor_matrix_avg_roi, 'LineWidth', 1.5, 'Color', color(1:3), 'DisplayName', display_name);
        
        % Apply transparency with the alpha function
        alpha(color(4));

        % Store the plot handle for the legend
        plot_handles = [plot_handles, hLine];
    end

    % Plot the no-odor response and store the handle
    hNoOdor = plot(tt, EMPTY_ROIindex, 'LineWidth', 1.5, 'DisplayName', 'no-odor', 'Color', 'k');
    plot_handles = [plot_handles, hNoOdor];

    % Set plot limits
    xlim([-0.5, 1.5]);
    ylim(odor_y_lim);

    % Add a vertical dashed line at x = 0 and horizontal line at y = 0
    line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
    line(xlim, [0, 0], 'Color', [0, 0, 0, 0.2], 'LineStyle', '-', 'HandleVisibility', 'off');

    % Add labels and titles
    if strcmp(layer_acr, 'Glom')
        % no x-label for Glom layer
    elseif strcmp(layer_acr, 'MC')
        xlabel('Time [s]', 'FontSize', axis_size);  
    end
    
    ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', axis_size);

    % Create the legend using the plot handles
    legend(plot_handles, 'Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none', 'FontSize', legend_size);
    
    hold off;

    % Add title depending on the layer
    if strcmp(layer_acr, 'Glom')
        sgtitle('Selected Glomerulus');
    elseif strcmp(layer_acr, 'MC')
        sgtitle('Averaged Daughter Mitral Cells');
    end

    % Save the figure
    saveas(gcf, fullfile(plot_dir, sprintf('%s_ODOR-ONLY_plot.png', layer_acr)));

end

%% ODOR-ONLY [ODOR_ONLY_avg]

for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
         selected_columns = ODOR_ONLY{idx_odors, idx_conc}(:, ROIindex);
         ODOR_ONLY_roi_avg{idx_odors, idx_conc} = mean(selected_columns, 2); 

        % Loop through each pair of indices to compute the average
        for idx = 1:numel(stim_frame)
            % Get the start and end indices for the current range
            start_idx = stim_frame(idx); % Convert to 1-based index
            end_idx = stim_frame_end(idx);
            
            % Average the values in EMPTY_roi from start_idx to end_idx
            ODOR_ONLY_avg{idx_odors, idx_conc}(idx) = mean(ODOR_ONLY_roi_avg{idx_odors, idx_conc}(start_idx:end_idx));
        end
    end
end

%% ___________________________ ODOR - STIM ________________________________
%% ODOR + STIM [ODOR_STIM_avg]

ODOR_STIM_roi = cell(n_odors, n_conc); % 1x2 cell
ODOR_STIM_avg = cell(n_odors, n_conc); % 1x2 cell

for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
        % Check if the current cell in ODOR_STIM is empty
        if isempty(ODOR_STIM{idx_odors, idx_conc})
            continue; % Skip to the next iteration if the cell is empty
        end
        
        % Initialize a 1x6 cell for each concentration
        ODOR_STIM_roi{idx_odors, idx_conc} = cell(1, n_latencies);
        ODOR_STIM_avg{idx_odors, idx_conc} = zeros(1, n_latencies); % Use a numeric array to store averages directly

        for idx_latency = 1:n_latencies
            
            data_values = ODOR_STIM{idx_odors, idx_conc}{idx_latency}; % Accessing as a 180x27 double
            
            % Select the columns specified by ROIindex and compute the mean
            selected_columns = data_values(:, ROIindex); % Select specified columns
            ODOR_STIM_roi{idx_odors, idx_conc}{idx_latency} = mean(selected_columns, 2); % Average across rows

            % Get the start and end indices for the current range
            start_idx = stim_frame(idx_latency); % Convert to 1-based index
            end_idx = stim_frame_end(idx_latency);

            ODOR_STIM_avg{idx_odors, idx_conc}(idx_latency) = mean(ODOR_STIM_roi{idx_odors, idx_conc}{1, idx_latency}(start_idx:end_idx));

        end
    end
end

%% ODOR+STIM FLUORESCENCE TRACE PLOTS

% Define colors with alpha values
high_concentration_color = [0, 0, 1, 1.0];  % Blue color for higher concentration (RGBA)
low_concentration_color = [0, 0.5, 1, 1.0];   % Same blue for lower concentration (RGBA)

% Create a new figure for all ROIs and concentrations
figure_handle = figure('Name', 'All ODOR+STIM Fluorescence Trace Plots');
% Set the figure dimensions: [left, bottom, width, height]
if strcmp(layer_acr, 'MC')
    set(figure_handle, 'Position', [10, 100, 1000, 300]);
elseif strcmp(layer_acr, 'Glom')
    set(figure_handle, 'Position', [10, 100, 1000, 300]);
end

% Create a tiled layout with n_ROIs rows and 2 columns for the concentrations
tiledlayout(1, n_conc, 'Padding', 'compact');

for idx_odor = 1:n_odors
                  
    for idx_conc = 1:n_conc
        % Get the current odor stimulation data
        current_ODOR_STIM = ODOR_STIM_roi{idx_odor, idx_conc}; % Use () since it's a double, not a cell

        % Create the next tile for the current ROI and concentration
        nexttile;

        hold on;

        for idx_latency = 1:n_latencies
            % Extract the current latency's data using parentheses for double array indexing
            current_ODOR_STIM_lat = current_ODOR_STIM{idx_latency}(:,1);  % Use () since it's a double, not a cell

            legend_label = sprintf('%s ms', unique_latencies{idx_latency});

            % Set the alpha value based on the latency index
            alpha_value = 1 - (idx_latency - 1) * 0.15;

            % Set color based on concentration
            if idx_conc == 1 % Low concentration
                color = low_concentration_color;
            else % High concentration
                color = high_concentration_color;
            end

            % Create a plot with transparency
            h = line(tt, current_ODOR_STIM_lat, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', color);
            h.Color(4) = alpha_value;  % Set the alpha value
        end

        % Plot the no_stim data with fixed color
        plot(tt, EMPTY_roi, 'LineWidth', 1.5, 'DisplayName', 'no stim', 'Color', 'k');

        % Set limits for x and y axes
        xlim([-0.5, 1]);
        ylim(stim_odor_y_lim); % Adjusted to [0, 1] since the data is normalized

        % Dashed vertical line at x = 0
        line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
        % Horizontal line at y = 0
        line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

        % Add labels
        if strcmp(layer_acr, 'MC')
            xlabel('Time [s]','FontSize', axis_size);
        elseif strcmp(layer_acr, 'Glom')
            % No x-axis label if layer_acr is 'Glom'
        end
        
        if idx_conc == 1
            ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', axis_size);
        elseif idx_conc == 2
            % No x-axis label if idx_coc is associated with the high concentration
        end
        
        % Conditionally set the title based on layer_acr
        if strcmp(layer_acr, 'Glom')
            title(sprintf('%s Conc. %s%%', unique_odors{idx_odor}, unique_concentration{idx_conc}), 'FontSize', title_size);
        elseif strcmp(layer_acr, 'MC')
            % No title if layer_acr is 'MC'
        end
        
        % Add the legend conditionally
        lgd = legend('show', 'Location', 'northeastoutside', 'FontSize', legend_size); 
        lgd.EdgeColor = 'none';  % Remove the outline from the legend

        % For plots not in the first row, set legend text to white to appear "blank"
        if ROI_idx ~= ROIindex(1)
            lgd.Visible = 'off';

        end
    end
end

% Save the figure
saveas(gcf, fullfile(plot_dir, sprintf('%s_ODOR-STIM_dF_plot.png', layer_acr)));


if 1
%% ________________________ STIM-TIME AVG. PLOTS __________________________
%% COMBINED PLOTS - ORGANIZE DATA

% EMPTY_frame_avg - 1x6 double
% STIM_ONLY_avg - 1x6 double
% ODOR_ONLY_avg - 1x2 cell - 1x6 double
%%% Odor_1_Conc_1
%%% Odor_1_Conc_2
% ODOR_STIM_avg - 1x2 cel - 1x6 double
%%% Odor_1_Conc_1
%%% Odor_1_Conc_2

% Create the struct for storing plot data
ODOR_STIM_data = struct();

% Assign EMPTY_frame_avg and STIM_ONLY_avg to the respective fields
ODOR_STIM_data.Empty = EMPTY_frame_avg;
ODOR_STIM_data.StimOnly = STIM_ONLY_avg;

% Initialize fields for storing odor and stimulus data
for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
        % Create field names dynamically based on indices
        odor_field_name = sprintf('Odor%dConc%d', idx_odors, idx_conc);
        
        % Assign the respective averages to the fields for Odor
        ODOR_STIM_data.(odor_field_name) = ODOR_ONLY_avg{idx_odors, idx_conc};
    end
end

% Now assign the OdorStim fields after the Odor fields
for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
        % Create field names dynamically based on indices
        odor_stim_field_name = sprintf('OdorStim%dConc%d', idx_odors, idx_conc);
        
        % Assign the respective averages to the fields for OdorStim
        ODOR_STIM_data.(odor_stim_field_name) = ODOR_STIM_avg{idx_odors, idx_conc};
    end
end

%% UNCERTAINTIES

EMPTY_std = zeros(size(EMPTY_frame_avg));
EMPTY_sem = zeros(size(EMPTY_frame_avg));

% Compute SEM for EMPTY
for i = 1:length(ROIindex)
    for j = 1:n_latencies
        r = ROIindex(i);  
        % Calculate the standard deviation across the specified time points for each ROI
        EMPTY_std(i,j) = std(EMPTY(stim_frame(j):stim_frame_end(j), r), 0, 1);  % Standard deviation for the specified time range and ROI
    end
end

% Calculate the SEM for STIM_ONLY_numeric
EMPTY_sem = EMPTY_std / sqrt(n_trials);
EMPTY_sem_avg = mean(EMPTY_sem, 1);

% Initialize arrays to hold standard deviations and SEMs
STIM_ONLY_std = zeros(size(STIM_ONLY_avg));
STIM_ONLY_sem = zeros(size(STIM_ONLY_avg));

% Calculate the mean and SEM across the second dimension of STIM_ONLY_numeric
% The SEM for STIM_ONLY_numeric needs to incorporate the averages conducted previously
% In this case, the values in STIM_ONLY_numeric are averaged over 20 trials for STIM_RESPONSE_MEAN and 15 time points for POST_STIM_MEAN.

for i = 1:length(ROIindex)
    for j = 1:n_latencies
        % Get the specific ROI index
        r = ROIindex(i);  
        % Calculate the standard deviation across the specified time points for each ROI
        STIM_ONLY_std(i,j) = std(STIM_ONLY{1,j}(stim_frame(j):stim_frame_end(j), r), 0, 1);  % Standard deviation for the specified time range and ROI
    end
end

% Calculate the SEM for STIM_ONLY_numeric
STIM_ONLY_sem = STIM_ONLY_std / sqrt(n_trials);
STIM_ONLY_sem_avg = mean(STIM_ONLY_sem, 1);

% Compute SEM for ODOR_ONLY
% Initialize arrays to hold standard deviations and SEMs
ODOR_ONLY_std = cell(1, n_conc);
ODOR_ONLY_sem = cell(1, n_conc);
ODOR_ONLY_sem_avg = cell(1, n_conc);

for idx_conc = 1:n_conc
    % Initialize arrays for this concentration
    ODOR_ONLY_std{idx_conc} = zeros(length(ROIindex), n_latencies);
    ODOR_ONLY_sem{idx_conc} = zeros(length(ROIindex), n_latencies);
    
    for i = 1:length(ROIindex)
        for j = 1:n_latencies
            % Get the specific ROI index
            r = ROIindex(i);  
            % Calculate the standard deviation across the specified time points for each ROI
            % Extract the data for the specified range and ROI
            response_data = ODOR_ONLY{1,idx_conc}(stim_frame(j):stim_frame_end(j), r);
            
            % Calculate the standard deviation for this data
            ODOR_ONLY_std{idx_conc}(i, j) = std(response_data, 0, 1);  % 0 for normalization by (n-1)
        end
    end
    
    % Calculate the SEM for this concentration by dividing the std by sqrt(n_trials)
    ODOR_ONLY_sem{1,idx_conc} = ODOR_ONLY_std{idx_conc} / sqrt(n_trials);
    
    % Average the SEM across ROIs to get a single value for each latency
    ODOR_ONLY_sem_avg{idx_conc} = mean(ODOR_ONLY_sem{idx_conc}, 1);
end


% Compute SEM for ODOR_STIM
% Initialize arrays to hold standard deviations and SEMs
ODOR_STIM_std = cell(1, n_conc);
ODOR_STIM_sem = cell(1, n_conc);
ODOR_STIM_sem_avg = cell(1, n_conc);

for idx_conc = 1:n_conc
    % Initialize arrays for this concentration
    ODOR_STIM_std{1, idx_conc} = zeros(length(ROIindex), n_latencies);
    ODOR_STIM_sem{1, idx_conc} = zeros(length(ROIindex), n_latencies);
    
    for i = 1:length(ROIindex)
        for j = 1:n_latencies
            % Get the specific ROI index
            r = ROIindex(i);  
            % Calculate the standard deviation across the specified time points for each ROI
            % Extract the data for the specified range and ROI
            response_data = ODOR_STIM{1,idx_conc}{j}(stim_frame(j):stim_frame_end(j), r);
            
            % Calculate the standard deviation for this data
            ODOR_STIM_std{1,idx_conc}(i, j) = std(response_data, 0, 1);  % 0 for normalization by (n-1)
        end
    end
    
    % Calculate the SEM for this concentration by dividing the std by sqrt(n_trials)
    ODOR_STIM_sem{1, idx_conc} = ODOR_STIM_std{idx_conc} / sqrt(n_trials);
    
    % Average the SEM across ROIs to get a single value for each latency
    ODOR_STIM_sem_avg{1,idx_conc} = mean(ODOR_STIM_sem{1, idx_conc}, 1);
end


%% ODOR_STIM data SEM for fill function in following plots

% Create the struct for storing plot data
ODOR_STIM_data_sem = struct();

% Assign EMPTY_frame_avg and STIM_ONLY_avg to the respective fields
ODOR_STIM_data_sem.Empty = EMPTY_sem_avg;
ODOR_STIM_data_sem.StimOnly = STIM_ONLY_sem_avg;

% Initialize fields for storing odor and stimulus data
for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
        % Create field names dynamically based on indices
        odor_field_name = sprintf('Odor%dConc%d', idx_odors, idx_conc);
        
        % Assign the respective averages to the fields for Odor
        ODOR_STIM_data_sem.(odor_field_name) = ODOR_ONLY_sem_avg{idx_odors, idx_conc};
    end
end

% Now assign the OdorStim fields after the Odor fields
for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
        % Create field names dynamically based on indices
        odor_stim_field_name = sprintf('OdorStim%dConc%d', idx_odors, idx_conc);
        
        % Assign the respective averages to the fields for OdorStim
        ODOR_STIM_data_sem.(odor_stim_field_name) = ODOR_STIM_sem_avg{idx_odors, idx_conc};
    end
end

%% STIM-ONLY STIM-TIME AVG 

% no stim or odor (black)
% stim (red)

figure_handle = figure('Name', ['STIM-ONLY STIM-TIME AVG']);

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 500]);

% Hold on to plot multiple lines
hold on;

% Define custom legend labels and colors
custom_labels = {
    'no stim nor odor',      [0 0 0 1];           % ODOR_STIM_data.Empty (black)
    'stim',                 [1 0 0 1];           % ODOR_STIM_data.StimOnly (red)
};

% Get the number of fields in ODOR_STIM_data (excluding StimLatencies)
field_names = fieldnames(ODOR_STIM_data);
num_fields = length(field_names); 

% Store scatter handles for legend
legend_handles = 2; % Preallocate for legend handles

% Plot each field against StimLatencies with larger dots
for i = 1:num_fields
    % Get the field name
    field_name = field_names{i};
    
    % Determine the corresponding label and color
    switch field_name
        case 'Empty'
            label = custom_labels{1, 1};
            color = custom_labels{1, 2}; % Black
        case 'StimOnly'
            label = custom_labels{2, 1};
            color = custom_labels{2, 2}; % Red
        otherwise
            continue; % Skip any fields not explicitly listed
    end
    
    %FILL FUNCTION
    % Calculate the SEM bounds
    mean_values = ODOR_STIM_data.(field_name);
    sem_values = ODOR_STIM_data_sem.(field_name);
    upper_bound = mean_values + sem_values;
    lower_bound = mean_values - sem_values;
    
    % Plot SEM shaded region using fill
    fill([stim_time, fliplr(stim_time)], [upper_bound, fliplr(lower_bound)], ...
         color(1:3), 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 20% transparency for SEM
    
    % Plot lines with alpha transparency
    plot(stim_time, mean_values, 'Color', [color(1:3), color(4)], ...
         'LineWidth', 1, 'LineStyle', '-');

    % Plot markers with alpha transparency
    scatter_handle = scatter(stim_time, mean_values, 25, ...
            'MarkerEdgeColor', color(1:3), ...
            'MarkerFaceColor', color(1:3), ...
            'MarkerFaceAlpha', color(4)); % Set marker face alpha

    % Store the scatter handle for the legend
    legend_handles(i) = scatter_handle;
    scatter_handle.DisplayName = label; % Assign label for legend
end

% Dashed vertical line at x = 0
line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

ylim(stim_time_avg_y_lim); % Adjusted to [-0.5, 0.5] since the data is normalized

xticks(stim_time);

% Add labels and title

if strcmp(layer_acr, 'MC')
    xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
end
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex', 'FontSize', axis_size); % Updated y-label

% Add legend without an outline, using only the scatter handles
lgd = legend(legend_handles, 'Location', 'northeast', 'FontSize', legend_size);
lgd.EdgeColor = 'none'; % Remove the outline from the legend

% Add a main title for the entire figure
if strcmp(layer_acr, 'MC')
    sgtitle('Daughter Mitral Cells', 'FontSize', axis_size);
elseif strcmp(layer_acr, 'Glom')
    sgtitle('Selected Glomerulus', 'FontSize', axis_size);
end

% Save the figure
saveas(gcf, fullfile(plot_dir, sprintf('%s_STIM-ONLY_stim-time_avg_plot.png', layer_acr)));

% Hold off to finish plotting
hold off;



%% STIM-ONLY & ODOR-ONLY STIM-TIME AVG.

% no stim or odor (black)
% stim (red)
% low BENZ (blue) (0.5 alpha)
% high BENZ (blue) (1.0 alpha)
% stim + low BENZ (purple) (0.5 alpha)
% stim + high BENZ (purple) (1.0 alpha)

figure_handle = figure('Name', ['STIM-ONLY & ODOR-ONLY STIM-TIME AVG.']);

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 500]);

% Hold on to plot multiple lines
hold on;

% Define custom legend labels and colors
custom_labels = {
    'no stim nor odor',           [0 0 0 1];           % ODOR_STIM_data.Empty (black)
    'stim',                       [1 0 0 1];            % ODOR_STIM_data.StimOnly (red)
    'low Hexanal conc.',          [0 0.5 1 1];          % ODOR_STIM_data.Odor1Conc1 (blue, 0.5 alpha)
    'high Hexanal conc.',         [0 0 1 1];            % ODOR_STIM_data.Odor1Conc2 (blue, 1.0 alpha)
};

% Get the number of fields in ODOR_STIM_data (excluding StimLatencies)
field_names = fieldnames(ODOR_STIM_data);
num_fields = 4; 

% Store scatter handles for lege(num_fields, 1)nd
legend_handles = gobjects; % Preallocate for legend handles

% Plot each field against StimLatencies with larger dots
for i = 1:num_fields
    % Get the field name
    field_name = field_names{i};
    
    % Determine the corresponding label and color
    switch field_name
        case 'Empty'
            label = custom_labels{1, 1};
            color = custom_labels{1, 2}; % Black
        case 'StimOnly'
            label = custom_labels{2, 1};
            color = custom_labels{2, 2}; % Red
        case 'Odor1Conc1'
            label = custom_labels{3, 1};
            color = custom_labels{3, 2}; % Low BENZ (blue, 0.5 alpha)
        case 'Odor1Conc2'
            label = custom_labels{4, 1};
            color = custom_labels{4, 2}; % High BENZ (blue, 1.0 alpha)
        otherwise
            continue; % Skip any fields not explicitly listed
    end
    
    %FILL FUNCTION
    % Calculate the SEM bounds
    mean_values = ODOR_STIM_data.(field_name);
    sem_values = ODOR_STIM_data_sem.(field_name);
    upper_bound = mean_values + sem_values;
    lower_bound = mean_values - sem_values;
    
    % Plot SEM shaded region using fill
    fill([stim_time, fliplr(stim_time)], [upper_bound, fliplr(lower_bound)], ...
         color(1:3), 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 20% transparency for SEM
    
    % Plot lines with alpha transparency
    plot(stim_time, mean_values, 'Color', [color(1:3), color(4)], ...
         'LineWidth', 1, 'LineStyle', '-');

    % Plot markers with alpha transparency
    scatter_handle = scatter(stim_time, ODOR_STIM_data.(field_name), 25, ...
            'MarkerEdgeColor', color(1:3), ...
            'MarkerFaceColor', color(1:3), ...
            'MarkerFaceAlpha', color(4)); % Set marker face alpha

    % Store the scatter handle for the legend
    legend_handles(i) = scatter_handle;
    scatter_handle.DisplayName = label; % Assign label for legend
end

% Dashed vertical line at x = 0
line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

ylim(stim_time_avg_y_lim); % Adjusted to [-0.5, 0.5] since the data is normalized

xticks(stim_time);

% Add labels and title

if strcmp(layer_acr, 'MC')
    xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
end
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex', 'FontSize', axis_size); % Updated y-label

% Add legend without an outline, using only the scatter handles
lgd = legend(legend_handles, 'Location', 'northeast', 'FontSize', legend_size);
lgd.EdgeColor = 'none'; % Remove the outline from the legend

% Add a main title for the entire figure
if strcmp(layer_acr, 'MC')
    sgtitle('Daughter Mitral Cells', 'FontSize', axis_size);
elseif strcmp(layer_acr, 'Glom')
    sgtitle('Selected Glomerulus', 'FontSize', axis_size);
end

% Save the figure
saveas(gcf, fullfile(plot_dir, sprintf('%s_ODOR-ONLY_STIM-ONLY_stim-time_avg_plot.png', layer_acr)));

% Hold off to finish plotting
hold off;



%% ODOR-STIM STIM-TIME AVG.

% no stim or odor (black)
% stim (red)
% low BENZ (blue) (0.5 alpha)
% high BENZ (blue) (1.0 alpha)
% stim + low BENZ (purple) (0.5 alpha)
% stim + high BENZ (purple) (1.0 alpha)

figure_handle = figure('Name', ['ODOR+STIM Plot']);

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 500]);

% Hold on to plot multiple lines
hold on;

% Define custom legend labels and colors
custom_labels = {
    'no stim nor odor',      [0 0 0 1];           % ODOR_STIM_data.Empty (black)
    'stim',                 [1 0 0 1];            % ODOR_STIM_data.StimOnly (red)
    'low Hexanal conc.',          [0 0.5 1 1];          % ODOR_STIM_data.Odor1Conc1 (blue, 0.5 alpha)
    'high Hexanal conc.',         [0 0 1 1];            % ODOR_STIM_data.Odor1Conc2 (blue, 1.0 alpha)
    'stim + low Hexanal conc.',   [0.7, 0.5, 0.8 1];      % ODOR_STIM_data.OdorStim1Conc1 (purple, 0.5 alpha)
    'stim + high Hexanal conc.',  [0.5 0 0.5 1];        % ODOR_STIM_data.OdorStim1Conc2 (purple, 1.0 alpha)
};

% Get the number of fields in ODOR_STIM_data (excluding StimLatencies)
field_names = fieldnames(ODOR_STIM_data);
num_fields = length(field_names); 

% Store scatter handles for legend
legend_handles = gobjects(num_fields, 1); % Preallocate for legend handles

% Plot each field against StimLatencies with larger dots
for i = 1:num_fields
    % Get the field name
    field_name = field_names{i};
    
    % Determine the corresponding label and color
    switch field_name
        case 'Empty'
            label = custom_labels{1, 1};
            color = custom_labels{1, 2}; % Black
        case 'StimOnly'
            label = custom_labels{2, 1};
            color = custom_labels{2, 2}; % Red
        case 'Odor1Conc1'
            label = custom_labels{3, 1};
            color = custom_labels{3, 2}; % Low BENZ (blue, 0.5 alpha)
        case 'Odor1Conc2'
            label = custom_labels{4, 1};
            color = custom_labels{4, 2}; % High BENZ (blue, 1.0 alpha)
        case 'OdorStim1Conc1'
            label = custom_labels{5, 1};
            color = custom_labels{5, 2}; % Stim + low BENZ (purple, 0.5 alpha)
        case 'OdorStim1Conc2'
            label = custom_labels{6, 1};
            color = custom_labels{6, 2}; % Stim + high BENZ (purple, 1.0 alpha)
        otherwise
            continue; % Skip any fields not explicitly listed
    end
    
   %FILL FUNCTION
    % Calculate the SEM bounds
    mean_values = ODOR_STIM_data.(field_name);
    sem_values = ODOR_STIM_data_sem.(field_name);
    upper_bound = mean_values + sem_values;
    lower_bound = mean_values - sem_values;
    
    % Plot SEM shaded region using fill
    fill([stim_time, fliplr(stim_time)], [upper_bound, fliplr(lower_bound)], ...
         color(1:3), 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % 20% transparency for SEM
    
    % Plot lines with alpha transparency
    plot(stim_time, mean_values, 'Color', [color(1:3), color(4)], ...
         'LineWidth', 1, 'LineStyle', '-');

    % Plot markers with alpha transparency
    scatter_handle = scatter(stim_time, ODOR_STIM_data.(field_name), 25, ...
            'MarkerEdgeColor', color(1:3), ...
            'MarkerFaceColor', color(1:3), ...
            'MarkerFaceAlpha', color(4)); % Set marker face alpha

    % Store the scatter handle for the legend
    legend_handles(i) = scatter_handle;
    scatter_handle.DisplayName = label; % Assign label for legend
end

% Dashed vertical line at x = 0
line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

ylim(stim_time_avg_y_lim); % Adjusted to [-0.5, 0.5] since the data is normalized

xticks(stim_time);

% Add labels and title

if strcmp(layer_acr, 'MC')
    xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
end
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex', 'FontSize', axis_size); % Updated y-label

% Add legend without an outline, using only the scatter handles
lgd = legend(legend_handles, 'Location', 'northeast', 'FontSize', legend_size);
lgd.EdgeColor = 'none'; % Remove the outline from the legend

% Add a main title for the entire figure
if strcmp(layer_acr, 'MC')
    sgtitle('Daughter Mitral Cells', 'FontSize', axis_size);
elseif strcmp(layer_acr, 'Glom')
    sgtitle('Selected Glomerulus', 'FontSize', axis_size);
end

% Save the figure
saveas(gcf, fullfile(plot_dir, sprintf('%s_ODOR-STIM_stim-time_avg_plot.png', layer_acr)));

% Hold off to finish plotting
hold off;

%% _________________________ COMBINED PLOTS _______________________________

%% COMBINED PLOTS - STIM-ONLY - AVERAGES ROI's


% Define the filenames for the two PNG files
file1 = sprintf('Glom_STIM-ONLY_avg_plot.png');
file2 = sprintf('MC_STIM-ONLY_avg_plot.png');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Create a blank space for text (adjust the height as needed)
text_height = 100; % Height for the text area
blank_space = uint8(255 * ones(text_height, max(size(img1, 2), size(img2, 2)), 3)); % White space
combined_img = [img1; img2; blank_space]; % Place the blank space at the bottom

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
% title('Single Glomerular Stimulation', 'FontSize', title_size);

% Calculate the vertical position for text placement
baseY = size(combined_img, 1) - text_height / 2;

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), baseY + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), baseY + 100, 'normalized with respect to OB layer & ROIs', ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');


% Save the figure
print(gcf, fullfile(plot_dir, sprintf('stimodor_STIM-ONLY_dF.png')), '-dpng', '-r300');


%% COMBINED PLOTS - ODOR-ONLY - dF traces

% Define the filenames for the two PNG files
file1 = sprintf('Glom_ODOR-ONLY_plot.png');
file2 = sprintf('MC_ODOR-ONLY_plot.png');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Create a blank space for text (adjust the height as needed)
text_height = 100; % Height for the text area
blank_space = uint8(255 * ones(text_height, max(size(img1, 2), size(img2, 2)), 3)); % White space
combined_img = [img1; img2; blank_space]; % Place the blank space at the bottom

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
% title(sprintf('Passive Odor Presentation'), 'FontSize', title_size);

% Calculate the vertical position for text placement
baseY = size(combined_img, 1) - text_height / 2;

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), baseY + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), baseY + 100, 'normalized with respect to OB layer & ROIs', ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');


% Save the figure
print(gcf, fullfile(plot_dir, sprintf('stimodor_ODOR-ONLY_dF_%s.png', chosen_odor)), '-dpng', '-r300');


%% COMBINED - ODOR+STIM FLUORESCENCE TRACE PLOTS

% Define the filenames for the two PNG files
file1 = sprintf('Glom_ODOR-STIM_dF_plot.png');
file2 = sprintf('MC_ODOR-STIM_dF_plot.png');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Define the vertical gap (space) between the images
gap = 0; % Adjust this value as needed for the desired spacing

% Create a blank space of the same width as the images and height of the gap
[height1, width1, ~] = size(img1);
[height2, width2, ~] = size(img2);
blank_space = uint8(255 * ones(gap, max(width1, width2), 3)); % White space

% Combine the images with the blank space
combined_img = [img1; blank_space; img2];

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
% title(sprintf('Odor & Single Glom Stimulation'), 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 100, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add vertical label for the first row ("Glomerular Layer")
text(-80, height1 / 2, 'Selected Glom', 'FontSize', legend_size, 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% Add vertical label for the second row ("Mitral Cell Layer")
text(-80, height1 + height2 / 2, 'Avg. Daughter MCs', 'FontSize', legend_size, 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');


% Save the figure
print(gcf, fullfile(plot_dir, sprintf('stimodor_STIM-ODOR_dF_%s.png', chosen_odor)), '-dpng', '-r300');


%% COMBINED PLOTS - STIM-ONLY STIM-TIME AVG.


% Define the filenames for the two PNG files
file1 = sprintf('Glom_STIM-ONLY_stim-time_avg_plot.png');
file2 = sprintf('MC_STIM-ONLY_stim-time_avg_plot.png');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Create a blank space for text (adjust the height as needed)
text_height = 100; % Height for the text area
blank_space = uint8(255 * ones(text_height, max(size(img1, 2), size(img2, 2)), 3)); % White space
combined_img = [img1; img2; blank_space]; % Place the blank space at the bottom

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
% title(sprintf('Passive Stim Presentation'), 'FontSize', title_size);

% Calculate the vertical position for text placement
baseY = size(combined_img, 1) - text_height / 2;

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), baseY + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), baseY + 100, 'normalized with respect to OB layer & ROIs', ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('stimodor_STIM-ONLY_stim-time_avg.png')), '-dpng', '-r300');

%% COMBINED PLOTS - STIM-ONLY & ODOR-ONLY STIM-TIME AVG.

% Define the filenames for the two PNG files
file1 = sprintf('Glom_ODOR-ONLY_STIM-ONLY_stim-time_avg_plot.png');
file2 = sprintf('MC_ODOR-ONLY_STIM-ONLY_stim-time_avg_plot.png');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Create a blank space for text (adjust the height as needed)
text_height = 100; % Height for the text area
blank_space = uint8(255 * ones(text_height, max(size(img1, 2), size(img2, 2)), 3)); % White space
combined_img = [img1; img2; blank_space]; % Place the blank space at the bottom

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
% title(sprintf('Passive Stim & Odor Presentation | %s', chosen_odor), 'FontSize', title_size);

% Calculate the vertical position for text placement
baseY = size(combined_img, 1) - text_height / 2;

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), baseY + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), baseY + 100, 'normalized with respect to OB layer & ROIs', ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('stimodor_ODOR-ONLY_STIM-ONLY_stim-time_avg_%s.png', chosen_odor)), '-dpng', '-r300');


%% COMBINED PLOTS - ODOR-STIM STIM-TIME AVG.


% Define the filenames for the two PNG files
file1 = sprintf('Glom_ODOR-STIM_stim-time_avg_plot.png');
file2 = sprintf('MC_ODOR-STIM_stim-time_avg_plot.png');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Create a blank space for text (adjust the height as needed)
text_height = 100; % Height for the text area
blank_space = uint8(255 * ones(text_height, max(size(img1, 2), size(img2, 2)), 3)); % White space
combined_img = [img1; img2; blank_space]; % Place the blank space at the bottom

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
% title(sprintf('Passive Stim & Odor Presentation | %s', chosen_odor), 'FontSize', title_size);

% Calculate the vertical position for text placement
baseY = size(combined_img, 1) - text_height / 2;

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), baseY + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), baseY + 100, 'normalized with respect to OB layer & ROIs', ...
    'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('stimodor_STIM-ODOR_stim-time_avg_%s.png', chosen_odor)), '-dpng', '-r300');

end
