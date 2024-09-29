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

%% DEFINE FILE PATHS

% FilePaths = {'mouse0773/240806_Glom_StimOdorMapping/aligned/240806_0773_Glom_StimOdorMapping_Z0_1X_S_v73.mat'};
% ROIindex = [10];

% FilePaths = {'mouse0773/240806_MC_StimOdorMapping/aligned/240806_MC_StimOdorMapping_Z186_3X_S_v73.mat'};
% ROIindex = [3,7,9];

% FilePaths = {'mouse0781/240806_Glom_StimOdorMapping/aligned/240806_0781_Glom_StimOdorMapping_Z0_1X_S_v73.mat'};
% ROIindex = [15]; 

FilePaths = {'mouse0781/240806_MC_StimOdorMapping/aligned/240806_0781_MC_StimOdorMapping_Z105_3X_S_v73.mat'};
ROIindex = [7,9,18,25];

% Extract the parent directory for mouse0773
commonDirectory = fileparts(fileparts(fileparts(FilePaths{1})));

% Create the common "plots" directory
plotsDirectory = fullfile(commonDirectory, 'plots');
if ~exist(plotsDirectory, 'dir')
    mkdir(plotsDirectory);
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

%% DEFINE KEY VARIABLES

% Load the data from the first file in FilePaths
data = load(FilePaths{1}, 'Session');

n_ROI = size(data.Session.OdorResponse{1,1}, 2);
%n_odors = 
%n_stim = 

% TIME VARIABLES
pre_inh = data.Session.Infos.pre_inh;
post_inh = data.Session.Infos.post_inh;
n_frames = pre_inh + post_inh;
n_fps = data.Session.Infos.fps;
tt = linspace(-(pre_inh/n_fps), (post_inh/n_fps), n_frames);

%% NORMALIZE ODOR RESPONSE

ODOR_RESPONSE = data.Session.OdorResponse;

% average accross trials (3rd dimention in double)
ODOR_RESPONSE_avg = cell(size(ODOR_RESPONSE));
ODOR_RESPONSE_avg = cellfun(@(x) mean(x, 3), ODOR_RESPONSE, 'UniformOutput', false);

ODOR_RESPONSE_gmax = -Inf;

% Loop through the averaged cells to find the global maximum
for i = 1:numel(ODOR_RESPONSE_avg)
    currentMax = max(ODOR_RESPONSE_avg{i}(:)); % Find the max in the current matrix
    if currentMax > ODOR_RESPONSE_gmax
        ODOR_RESPONSE_gmax = currentMax;
    end
end

% Normalize ODOR_RESPONSE_avg
ODOR_RESPONSE_norm = cell(size(ODOR_RESPONSE_avg));
for i = 1:numel(ODOR_RESPONSE_avg)
    ODOR_RESPONSE_norm{i} = ODOR_RESPONSE_avg{i} / ODOR_RESPONSE_gmax;
end


%% Replace gating frame fluorescence value with average of before and after frame

% Define the range of row indexes
row_range = pre_inh:(pre_inh + 8); % This will be rows 60 to 68

% Initialize a cell array to store the indices of the minimum values
min_indices = cell(1, numel(ODOR_RESPONSE_norm));

% Loop through each cell in ODOR_RESPONSE_norm
for idx = 1:numel(ODOR_RESPONSE_norm)
    % Extract the current 180x27 double
    current_data = ODOR_RESPONSE_norm{idx};
    
    % Initialize an array to store the indices of the minimum values for the current cell
    min_indices_current = zeros(1, size(current_data, 2)); % Preallocate for each column
    
    % Loop through each column of the current data
    for col_idx = 1:size(current_data, 2)
        % Find the index of the minimum value in the specified row range
        [~, min_row_index] = min(current_data(row_range, col_idx));
        
        % Store the index relative to the entire data (not just the row range)
        min_indices_current(col_idx) = min_row_index + (pre_inh-1); 
    end
    
    % Store the result for the current cell
    min_indices{idx} = min_indices_current;
end

updated_ODOR_RESPONSE_norm = cell(size(ODOR_RESPONSE_norm));

% Loop through each cell in ODOR_RESPONSE_norm
for idx = 1:numel(ODOR_RESPONSE_norm)
    % Extract the current 180x27 double
    current_data = ODOR_RESPONSE_norm{idx};
    
    % Get the indices of the minimum values for the current cell
    min_indices_current = min_indices{idx};
    
    % Loop through each column of the current data
    for col_idx = 1:size(current_data, 2)
        % Get the index of the minimum value for the current column
        min_index = min_indices_current(col_idx);
        
        % Ensure the index is valid for averaging (not at the boundaries)
        if min_index > 1 && min_index < size(current_data, 1) % Check bounds
            % Get the values before and after the minimum index
            before_value = current_data(min_index - 1, col_idx);
            after_value = current_data(min_index + 1, col_idx);
            
            % Calculate the average of the before and after values
            average_value = (before_value + after_value) / 2;
            
            % Replace the value at min_index with the average
            current_data(min_index, col_idx) = average_value;
        end
    end
    
    % Store the updated data in the new cell array
    updated_ODOR_RESPONSE_norm{idx} = current_data;
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
        EMPTY = updated_ODOR_RESPONSE_norm{i};
        
    % Check for STIM_ONLY condition
    elseif strcmp(cond_odors{i}, 'empty') && ~isempty(latency_idx)
        STIM_ONLY{latency_idx} = updated_ODOR_RESPONSE_norm{i};
    
    % Check for ODOR_ONLY condition
    elseif ~strcmp(cond_odors{i}, 'empty') && isempty(cond_latencies{i})
        odor_idx = find(strcmp(unique_odors, cond_odors{i}));
        if ~isempty(odor_idx) && concentration_idx > 0
            ODOR_ONLY{odor_idx, concentration_idx} = updated_ODOR_RESPONSE_norm{i};
        end
    
    % Check for ODOR_STIM condition
    elseif ~strcmp(cond_odors{i}, 'empty') && ~isempty(latency_idx)
        odor_idx = find(strcmp(unique_odors, cond_odors{i}));
        if ~isempty(odor_idx) && concentration_idx > 0
            ODOR_STIM{odor_idx, concentration_idx}{latency_idx} = updated_ODOR_RESPONSE_norm{i};
        end
    end
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
        h = line(tt, currentLatency, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', 'b');
        h.Color(4) = alpha_value;  % Set the alpha value

        hold on;

        xlim([-0.5, 1]);
        ylim([-0.5, 1]); % Adjusted to [0, 1] since the data is normalized

        % Dashed vertical line at x = 0
        line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');

        line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

        % Add labels, title, etc. as needed
        xlabel('Time [s]');
        ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}','Interpreter','latex','FontName','Arial','FontSize',10);

    end
    
    % Plot the no_stim data for the current ROI
    no_stim_ROI = EMPTY(:, ROI_idx); 
    
    % Plot the no_stim data
    plot(tt, no_stim_ROI, 'LineWidth', 1.5, 'DisplayName', 'no stim', 'Color', 'k');
        
    % Add legend
    legend_handle = legend('Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none');
    
    % Set the title for the entire figure
    sgtitle(['Single Glomerular Stimulation | ', layer_acr, ' ', num2str(ROI_idx)]);
    
    saveas(gcf, fullfile(plotsDirectory, sprintf('STIM-ONLY_%s_%d_plot.png', layer_acr, ROI_idx)));

end


%% ODOR-ONLY PLOTS
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
    ylim([-0.5, 1]);

    % Dashed vertical line at x = 0
    line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
    % Horizontal line at y = 0
    line(xlim, [0, 0], 'Color', [0, 0, 0, 0.2], 'LineStyle', '-', 'HandleVisibility', 'off');

    % Add labels, title, etc. as needed
    xlabel('Time [s]');
    ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 10);

    % Add legend with the defined entries
    legend([hLow, hHigh, hNoOdor], {'low conc.', 'high conc.', 'no-odor'}, 'Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none');
    
    hold off;

    % Set the title for the entire figure
    sgtitle(sprintf('Odor Presentation | %s, %s', layer_name, chosen_odor));

    % Save the figure
    saveas(gcf, fullfile(plotsDirectory, sprintf('ODOR-ONLY_%s_%s_plot.png', layer_acr, chosen_odor)));

end


%% ODOR+STIM PLOTS
% (add low and high concentrations)

% Define colors with alpha values
high_concentration_color = [0, 0, 1, 1.0];  % Blue color for higher concentration (RGBA)
low_concentration_color = [0, 0.5, 1, 1.0];   % Same blue for lower concentration (RGBA)

% Create a new figure for all ROIs and concentrations
figure_handle = figure('Name', 'All ODOR+STIM Plots');
% Set the figure dimensions: [left, bottom, width, height]
if strcmp(layer_acr, 'MC')
    set(figure_handle, 'Position', [10, 100, 1000, 1000]);
elseif strcmp(layer_acr, 'Glom')
    set(figure_handle, 'Position', [10, 100, 1000, 300]);
end

% Create a tiled layout with n_ROIs rows and 2 columns for the concentrations
tiledlayout(numel(ROIindex), n_conc, 'Padding', 'compact');

for idx_odor = 1:n_odors
    for ROI_idx = ROIindex % Iterate through each ROI
                  
        for idx_conc = 1:n_conc
            % Get the current odor stimulation data
            current_ODOR_STIM = ODOR_STIM{idx_odor, idx_conc};

            % Create the next tile for the current ROI and concentration
            nexttile;

            hold on;

            for idx_latency = 1:n_latencies
                % Get the column for the current latency and ROI
                ODOR_column = current_ODOR_STIM{idx_latency}(:, ROI_idx);

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
                h = line(tt, ODOR_column, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', color);
                h.Color(4) = alpha_value;  % Set the alpha value
            end

            % Plot the no_stim data for the current ROI
            no_stim_ROI = EMPTY(:, ROI_idx); 

            % Plot the no_stim data with fixed color
            plot(tt, no_stim_ROI, 'LineWidth', 1.5, 'DisplayName', 'no stim', 'Color', 'k');

            % Set limits for x and y axes
            xlim([-0.5, 1]);
            ylim([-0.5, 1]); % Adjusted to [0, 1] since the data is normalized

            % Dashed vertical line at x = 0
            line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
            % Horizontal line at y = 0
            line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

            % Add labels
            xlabel('Time [s]');
            ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}', 'Interpreter', 'latex', 'FontName', 'Arial', 'FontSize', 10);

            % Set the title for each subplot
            title(sprintf('%s %d - %s conc. %s%%', layer_acr, ROI_idx, unique_odors{idx_odor}, unique_concentration{idx_conc}));
            
            % Add a legend only for plots in the first row
            if ROI_idx == ROIindex(1)
                lgd = legend('show', 'Location', 'northwest'); % Show the legend for the first row only, positioned southeast
                lgd.EdgeColor = 'none'; % Remove the outline from the legend
            end
        end
    end
end

% Add a main title for the entire figure
sgtitle(sprintf('Single Glom Stimulation with Odor Background | %s', layer_name));

% Save the figure
saveas(gcf, fullfile(plotsDirectory, sprintf('ODOR-STIM_All_ROIs_%s_plot.png', layer_acr)));

%% COMBINED PLOTS - ORGANIZING DATA
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

%% STIM-ONLY [STIM_ONLY_avg]

for idx_latencies = 1:n_latencies
    selected_columns = STIM_ONLY{1,idx_latencies}(:, ROIindex);
    STIM_ONLY_roi{1,idx_latencies} = mean(selected_columns, 2); 

    for idx = 1:numel(stim_frame)
    % Get the start and end indices for the current range
    start_idx = stim_frame(idx); % Convert to 1-based index
    end_idx = stim_frame_end(idx);
    
    % Average the values in EMPTY_roi from start_idx to end_idx
    STIM_ONLY_avg(idx) = mean(STIM_ONLY_roi{1,idx_latencies}(start_idx:end_idx));
    end
end

%% ODOR-ONLY [ODOR_ONLY_avg]

for idx_odors = 1:n_odors
    for idx_conc = 1:n_conc
         selected_columns = ODOR_ONLY{idx_odors, idx_conc}(:, ROIindex);
         ODOR_ONLY_roi{idx_odors, idx_conc} = mean(selected_columns, 2); 

        % Loop through each pair of indices to compute the average
        for idx = 1:numel(stim_frame)
            % Get the start and end indices for the current range
            start_idx = stim_frame(idx); % Convert to 1-based index
            end_idx = stim_frame_end(idx);
            
            % Average the values in EMPTY_roi from start_idx to end_idx
            ODOR_ONLY_avg{idx_odors, idx_conc}(idx) = mean(ODOR_ONLY_roi{idx_odors, idx_conc}(start_idx:end_idx));
        end
    end
end

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

%% COMBINED PLOTS

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
    'no stim or odor',      [0 0 0 1];           % ODOR_STIM_data.Empty (black)
    'stim',                 [1 0 0 1];           % ODOR_STIM_data.StimOnly (red)
    'low Hexanal',          [0 0 1 0.4];         % ODOR_STIM_data.Odor1Conc1 (blue, 0.5 alpha)
    'high Hexanal',         [0 0 1 1];           % ODOR_STIM_data.Odor1Conc2 (blue, 1.0 alpha)
    'stim + low Hexanal',   [0.5 0 0.5 0.4];     % ODOR_STIM_data.OdorStim1Conc1 (purple, 0.5 alpha)
    'stim + high Hexanal',  [0.5 0 0.5 1];       % ODOR_STIM_data.OdorStim1Conc2 (purple, 1.0 alpha)
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
    
    % Plot lines with alpha transparency
    plot(stim_time, ODOR_STIM_data.(field_name), 'Color', [color(1:3), color(4)], ...
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

ylim([-0.25, 0.25]); % Adjusted to [-0.5, 0.5] since the data is normalized

xticks(stim_time);

% Add labels and title
xlabel('Stim Onset Latency [ms]');
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex'); % Updated y-label

% Add legend without an outline, using only the scatter handles
lgd = legend(legend_handles, 'Location', 'best');
lgd.EdgeColor = 'none'; % Remove the outline from the legend

% Add a main title for the entire figure
sgtitle(sprintf('Single Glom Stimulation with Odor Background | %s', layer_name));

% Save the figure
saveas(gcf, fullfile(plotsDirectory, sprintf('ODOR-STIM_avgd_ROIs_%s_plot.png', layer_acr)));

% Hold off to finish plotting
hold off;