clearvars
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
title_size = 14;

%% CHOSEN ODOR TO DISPLAY

chosen_odor = {'Hexanal'};
% chosen_odor = {'EthylTiglate','EthylButyrate','Hexanal','BenzAldehyde','MethylValerate'};

% Extended Odor Panel
% chosen_odor = {'2MBA', '5M2H','EthylButyrate','Hexanal','BenzAldehyde','Heptanal','MethylValerate','PropionicAcid'};
% chosen_odor = {'2MBA','Hexanal','BenzAldehyde','Heptanal'}; % chosen from 240731 session

% chosen_odor = {'Heptanal','Hexanal'}; 

%% DEFINE FILE PATHS

FilePaths = {'mouse0773/240703_Glom_100flow/aligned/240703_Glom_OdorMapping_Z0_1X_100flow_S_v73.mat',
             'mouse0773/240705_Glom_20flow/aligned/240705_Glom_OdorMapping_Z0_1X_20flow_S_v73.mat'};
ROIindex_single = [10];
% ROIindex = [2,5,4,7,10]; % BenzAldehyde
% cross_y_lim = [-0.3, 0.3] % BenzAldehyde
% avg_y_lim = [-0.3, 0.3] % BenzAldehyde
% ROIindex = [10,13,16]; % EthylButyrate
% cross_y_lim = [-0.2, 1] % Ethyl Butyrate
% avg_y_lim = [-0.2, 1] % Ethyl Butyrate
ROIindex = [10,11,12,14,15,17,18]; % Hexanal
cross_y_lim = [-0.2, 1] % Hexanal
avg_y_lim = [-0.3, 0.2] % Hexanal
% ROIindex = [7,9,10,11,12,13,16]; % MethylValerate
% cross_y_lim = [-0.2, 1] % MethylValerate
% avg_y_lim = [-0.3, 0.3] % MethylValerate

% ROIindex = [2,7,9,10,18]; % EthylTiglate - not great
% 
% FilePaths = {'mouse0773/240703_MC_100flow/aligned/240703_MC_OdorMapping_Z196_3X_100flow_S_v73.mat',
%               'mouse0773/240705_MC_20flow/aligned/240703_MC_OdorMapping_Z196_3X_20flow_S_v73.mat'};
% ROIindex = [4,7,9];
% ROIindex_single = [4,7,9];
% cross_y_lim = [-0.2, 1] % Ethyl Butyrate
% avg_y_lim = [-0.2, 1] % Ethyl Butyrate
% cross_y_lim = [-0.3, 0.3] % Hexanal
% avg_y_lim = [-0.3, 0.3] % Hexanal
% cross_y_lim = [-0.3, 0.3] % BenzAldehyde
% avg_y_lim = [-0.3, 0.3] % BenzAldehyde
% cross_y_lim = [-0.3, 0.3] % MethylValerate
% avg_y_lim = [-0.3, 0.3] % MethylValerate

% FilePaths = {'mouse0781/240705_Glom_100flow/aligned/240705_Glom_OdorMapping_Z0_1X_100flow_S_v73.mat',
%             'mouse0781/240705_Glom_20flow/aligned/240705_Glom_OdorMapping_Z0_1X_20flow_S_v73.mat'};
% ROIindex_single = [13];
% ROIindex = [12,13]; % BenzAldehyde
% ROIindex = [8,9,12,13]; % EthylButyrate
% ROIindex = [11,12,13,15,16,22,23]; % Hexanal
% ROIindex = [11,12,13,16]; % MethylValerate
% ROIindex = [12,13,24]; % EthylTiglate
% 
% FilePaths = {'mouse0781/240731_Glom_100flow/aligned/240731_Glom_OdorMapping_Z0_1X_100flow_S_v73.mat',
%              'mouse0781/240731_Glom_20flow/aligned/240731_0781_Glom_OdorMapping_Z0_1X_20flow_S_v73.mat'};
% ROIindex_single = [15];
% ROIindex = [7,9,8,10,11,15]; % 2MBA
% cross_y_lim = [-0.2, 1] % 2MBA
% avg_y_lim = [-0.2, 0.6] % 2MBA
% ROIindex = [7,11,15,18,20]; % Hexanal
% cross_y_lim = [-0.2, 1] % Hexanal
% avg_y_lim = [-0.2, 0.3] % Hexanal
% ROIindex = [7,10,11,15]; % BenzAldehyde
% cross_y_lim = [-0.2, 0.5] % BenzAldehyde
% avg_y_lim = [-0.2, 0.3] % BenzAldehyde
% ROIindex = [7,10,11,15,21,28]; % Heptanal
% cross_y_lim = [-0.5, 1] % Heptanal
% avg_y_lim = [-0.2, 0.4] % Heptanal
% 
% FilePaths = {'mouse0781/240801_MC_100flow/aligned/240801_0781_MC_OdorMapping_Z94_3X_100flow_S_v73.mat',
%              'mouse0781/240801_MC_20flow/aligned/240801_0781_MC_OdorMapping_Z94_3X_20flow_S_v73.mat'};
% ROIindex_single = [8,10,18,26];
% ROIindex = [8,10,18,26];
% cross_y_lim = [-0.2, 1] % 2MBA
% avg_y_lim = [-0.2, 0.6] % 2MBA
% cross_y_lim = [-0.3, 0.3] % Hexanal
% avg_y_lim = [-0.2, 0.3] % Hexanal
% cross_y_lim = [-0.3, 0.5] % Heptanal
% % avg_y_lim = [-0.3, 0.3] % Heptanal
% cross_y_lim = [-0.2, 0.5] % BenzAldehyde
% avg_y_lim = [-0.2, 0.3] % BenzAldehyde

% Extract the parent directory for mouse0773
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
    elseif contains(FilePaths{FilePathIndex}, 'MC')
        layer_name = 'Mitral Cell Layer';
    end

    if contains(FilePaths{FilePathIndex}, 'Glom')
        layer_acr = 'Glom';
    elseif contains(FilePaths{FilePathIndex}, 'MC')
        layer_acr = 'MC';
    end

end

%% PRINCIPLE SCALAR VALUES

n_FilePaths = numel(FilePaths);

for i = 1:n_FilePaths
    data = load(FilePaths{i}, 'Session');
end

n_odors = size(data.Session.OdorResponse, 2);
n_ROI = size(data.Session.OdorResponse{1,1}, 2);
n_trials = size(data.Session.OdorResponse{1,1}, 3);

% TIME VARIABLES
pre_inh = data.Session.Infos.pre_inh;
post_inh = data.Session.Infos.post_inh;
n_frames = pre_inh + post_inh;
n_fps = data.Session.Infos.fps;
tt = linspace(-(pre_inh/n_fps), (post_inh/n_fps), n_frames);

%% ORGANIZE FLUORESENCE RESPONSES
ODOR_RESPONSE_MEAN = cell(n_FilePaths,1);
ODOR_RESPONSE_UNC = cell(n_FilePaths,1);

odor = cell(n_odors, 1);
concentrations = [];
concentration = zeros(1, n_odors);

for idxFilePath = 1:n_FilePaths  % Loop through one concentration at a time
    % Load the session data
    load(FilePaths{idxFilePath});
   
    % Initialize concentrations as a double matrix on the first iteration
    if isempty(concentrations)
        concentrations = zeros(n_FilePaths, n_odors);
    end

    % Create cell arrays for mean and uncertainty with dimensions 1xn_odors
    OdorResponse_mean = cell(1, n_odors);
    OdorResponse_unc = cell(1, n_odors);
   
    for i = 1:n_odors
        matrix = Session.OdorResponse{i};
        % 'matrix' has 3 dimensions. Columns are different ROIs, rows are time,
        % the different layers are iterations (10)

        % Calculate the mean across the third dimension
        meanMatrix = mean(matrix, 3);

        % Calculate the uncertainty (standard error of the mean) across the third dimension
        stdMatrix = std(matrix, 0, 3);
        semMatrix = stdMatrix / sqrt(size(matrix, 3));

        % Store the results in the respective cell arrays
        OdorResponse_mean{i} = meanMatrix;
        OdorResponse_unc{i} = semMatrix;

    % RECONFIGURING UNIQUECONDS INTO ODOR AND CONCENTRATION CELLS

        str = Session.UniqueConds{i};

        % Split the string based on colon
        parts = strsplit(str, ':');

        % Extract odor (characters before the colon)
        odor{i} = parts{1};

        % Extract concentration (characters after the colon)
        concs = strsplit(parts{2}, '/');

        % Convert concentration values to numbers
        numbers = str2double(concs);

        % Remove NaN values
        non_nan_values = nonzeros(numbers);

        % Keep only one value
        if ~isempty(non_nan_values)
           non_nan_values = non_nan_values(1);
        end

        % Divide by 10
        non_nan_values = non_nan_values / 10;

        % Store the concentration value in the concentration variable
        concentration(i) = non_nan_values;
    end

    concentrations(idxFilePath, :) = concentration;
   
    %% EXTRACT ODOR NAMES
    % Initialize cell array to store substrings
    odorconcat = cell(n_odors, 1);

    % Loop through each element in 'odor'
    for i = 1:numel(odor)
        % Split the current string at each '/'
        substrings = strsplit(odor{i}, '/');

        % Store the substrings in 'odorconcat'
        odorconcat{i} = substrings;
    end

    % Convert 'odorconcat' to a cell array
    odorconcat = cellfun(@(x) x.', odorconcat, 'UniformOutput', false);

    % Convert the cell array to a matrix
    odorconcat = vertcat(odorconcat{:});

    % Define the values to be excluded
    excludeValues = {'--', 'xx', 'yy', 'zz'};

    % Create logical indices for entries to keep
    entriesToKeep = ~any(contains(odorconcat, excludeValues), 2);

    % Extract the entries to keep
    odor = odorconcat(entriesToKeep, :);
   
    % Initialize an array to store the indices
    chosen_odor_idx = zeros(length(chosen_odor), 1);

    % Loop through each chosen_odor and find the corresponding index in odor
    for i = 1:length(chosen_odor)
        % Find the index of the current chosen_odor in the odor array
        index = find(strcmp(odor, chosen_odor{i}));

        % Store the index in the chosen_odor_indices array
        if ~isempty(index)
            chosen_odor_idx(i) = index;
        else
            chosen_odor_idx(i) = NaN; % If the chosen_odor is not found, store NaN
        end
    end
   
    ODOR_RESPONSE_MEAN{idxFilePath} = OdorResponse_mean;
    ODOR_RESPONSE_UNC{idxFilePath} = OdorResponse_unc;
    
%     ODOR_RESPONSE_MEAN contains the averaged fluorescence traces for all odors between the n_trials
% 
%     ODOR_RESPONSE_MEAN is a 2x1 cell array
%     2 represents the number of concentrations being considered
%     Each cell contains another 1x6 cell
%     6 represents the number of odors being considered (n_odors)
%     Each of the 1x6 cells contains a 180X29 double
%     180 corresponds to the number of time points being considered (n_frames)
%     29 corresponds to the number of ROIs being considered (n_ROI)
   
end % end looping through both filepaths (idxFilePath)

%% NORMALIZE FLUORESENCE RESPONSES (ODOR_RESPONSE_MEAN)

global_max = -inf;
ODOR_RESPONSE_MEAN_norm = cell(size(ODOR_RESPONSE_MEAN));

% Find maximum value in ODOR_RESPONSE_MEAN
for row = 1:size(ODOR_RESPONSE_MEAN, 1)
    for col = 1:size(ODOR_RESPONSE_MEAN, 2)
        % Access the 1x6 cell within each element
        for inner_idx = 1:length(ODOR_RESPONSE_MEAN{row, col})
           
            current_matrix = ODOR_RESPONSE_MEAN{row, col}{inner_idx};
           
            current_max = max(current_matrix(:)); % Use (:) to convert matrix to a single vector
           
            % Update the global_max if current_max is larger
            if current_max > global_max
                global_max = current_max;
            end
        end
       
% ODOR_RESPONSE_MEAN_norm: Divide ODOR_RESPONSE_MEAN by its global maximum
        % Initialize the inner cell array to store the normalized matrices
        ODOR_RESPONSE_MEAN_norm{row, col} = cell(size(ODOR_RESPONSE_MEAN{row, col}));
       
        % Access the 1x6 cell within each element
        for inner_idx = 1:length(ODOR_RESPONSE_MEAN{row, col})
            % Access the 180x32 double matrix
            current_matrix = ODOR_RESPONSE_MEAN{row, col}{inner_idx};
           
            % Normalize the matrix using the global maximum value
            normalized_matrix = current_matrix / global_max;
           
            % Store the normalized matrix in the new cell array
            ODOR_RESPONSE_MEAN_norm{row, col}{inner_idx} = normalized_matrix;
        end
    end
end

% ODOR_RESPONSE_MEAN_norm divides measurements from all concentrations by the global maximum in ODOR_RESPONSE_MEAN
% 
% ODOR_RESPONSE_MEAN_norm is a 2x1 cell array
% 2 represents the number of concentrations being considered
% Each cell contains another 1x6 cell
% 6 represents the number of odors being considered (n_odors)
% Each of the 1x6 cells contains a 180X29 double
% 180 corresponds to the number of time points being considered (n_frames)
% 29 corresponds to the number of ROIs being considered (n_ROI)

%% NO-ODOR CODE

ODOR_RESPONSE_MEAN_norm_empty = cell(2,1);

for cell_array_index = 1:n_FilePaths
    ODOR_RESPONSE_MEAN_norm_empty{cell_array_index,1} = ODOR_RESPONSE_MEAN_norm{cell_array_index,1}{1,6};
end

ODOR_RESPONSE_MEAN_norm_empty_avg = (ODOR_RESPONSE_MEAN_norm_empty{1} + ODOR_RESPONSE_MEAN_norm_empty{2}) / 2;

ODOR_RESPONSE_MEAN_norm_empty_avg_rois = ODOR_RESPONSE_MEAN_norm_empty_avg(:, ROIindex_single);


%% PLOT CHECK
% 
% plot_check_high = ODOR_RESPONSE_MEAN_norm{2,1}{1,chosen_odor_idx}(:,ROIindex_single);
% plot_check_low = ODOR_RESPONSE_MEAN_norm{1,1}{1,chosen_odor_idx}(:,ROIindex_single);
% 
% % Create figure and plot
% figure;
% plot(tt, plot_check_high, 'LineWidth', 1.5);
% hold on;
% plot(tt, plot_check_low, 'LineWidth', 1.5);
% hold on;
% plot(tt, ODOR_RESPONSE_MEAN_norm_empty_avg_rois, 'LineWidth', 1.5);
% 
% % Plot formatting
% xlim([-0.5, 1.5]);
% ylim(cross_y_lim); % Make sure cross_y_lim is defined beforehand
% line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '-', 'LineWidth', 1, 'HandleVisibility', 'off');
% line(xlim, [0, 0], 'Color', [0, 0, 0, 0.2], 'LineStyle', '-', 'HandleVisibility', 'off');
% 
% % Add labels
% xlabel('Time [s]', 'FontSize', 12);
% ylabel('\Delta F/F_0', 'FontSize', 12);
% 
% % Optional: add a title
% title('Odor Response', 'FontSize', 14);


%% LOOP THROUGH CHOSEN ODORS FOR PLOTS

for k = 1:numel(chosen_odor_idx)
    idxChosenOdor = chosen_odor_idx(k);
    chosen_odor_conc = concentrations(:, idxChosenOdor);

    %% Time Cross-Section Plot
    [sorted_conc, sort_idx] = sort(chosen_odor_conc, 'descend');
    figure_handle = figure('Name', ['Chosen Odor ', num2str(idxChosenOdor)], 'Position', [10, 50, 600, 400]);
    hold on;

    % Initialize flags for legend labels
    added_single_high = false;
    added_single_low = false;
    added_neighbor_high = false;
    added_neighbor_low = false;

    handles = [];
    labels = [];

    for sorted_idx = 1:length(sort_idx)
        cell_array_index = sort_idx(sorted_idx);
        OdorResponse_mean = ODOR_RESPONSE_MEAN_norm{cell_array_index, 1}{1, idxChosenOdor};

        for j = 1:numel(ROIindex)
            roi_index = ROIindex(j);
            currentColumn = OdorResponse_mean(:, roi_index);

            % Determine legend label and color
            is_single = ismember(roi_index, ROIindex_single);
            is_low_conc = sorted_idx == 2;

            if is_single
                if strcmp(layer_acr, 'MC')
                    if is_low_conc
                        legend_label = 'daughter MCs - low conc.';
                    else
                        legend_label = 'daughter MCs - high conc.';
                    end
                else
                    if is_low_conc
                        legend_label = 'single Glom - low conc.';
                    else
                        legend_label = 'single Glom - high conc.';
                    end
                end
                color = [0, 0, 1, 0.4 + 0.6 * ~is_low_conc];
                linestyle = '-';

                % Plot single ROIs
                if is_low_conc && ~added_single_low
                    h_single_low = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', color, 'LineStyle', linestyle);
                    added_single_low = true;
                    handles(end+1) = h_single_low;
                    labels{end+1} = legend_label;
                elseif ~is_low_conc && ~added_single_high
                    h_single_high = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', color, 'LineStyle', linestyle);
                    added_single_high = true;
                    handles(end+1) = h_single_high;
                    labels{end+1} = legend_label;
                else
                    plot(tt, currentColumn, 'LineWidth', 1.5, 'Color', color, 'LineStyle', linestyle, 'HandleVisibility', 'off');
                end
            else
                if strcmp(layer_acr, 'Glom')
                    if is_low_conc
                        legend_label = 'neighbor Gloms - low conc.';
                    else
                        legend_label = 'neighbor Gloms - high conc.';
                    end
                else
                    legend_label = [layer_acr, ' ', num2str(roi_index)];
                end
                color = [0.5, 0.5, 0.5, 0.3 + 0.7 * ~is_low_conc];

                % Plot neighboring ROIs
                if is_low_conc && ~added_neighbor_low
                    h_neighbor_low = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', color);
                    added_neighbor_low = true;
                    handles(end+1) = h_neighbor_low;
                    labels{end+1} = legend_label;
                elseif ~is_low_conc && ~added_neighbor_high
                    h_neighbor_high = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', color);
                    added_neighbor_high = true;
                    handles(end+1) = h_neighbor_high;
                    labels{end+1} = legend_label;
                else
                    plot(tt, currentColumn, 'LineWidth', 1.5, 'Color', color, 'HandleVisibility', 'off');
                end
            end
        end
    end

    % Plot no-odor line
    h_no_odor = plot(tt, ODOR_RESPONSE_MEAN_norm_empty_avg_rois, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'no-odor');
    handles(end+1) = h_no_odor;
    labels{end+1} = 'no-odor';

    % Plot formatting
    xlim([-0.5, 1.5]);
    ylim(cross_y_lim);
    line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '-', 'LineWidth', 1, 'HandleVisibility', 'off');
    line(xlim, [0, 0], 'Color', [0, 0, 0, 0.2], 'LineStyle', '-', 'HandleVisibility', 'off');

    ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}','Interpreter','latex','FontName','Arial','FontSize', axis_size);
    if strcmp(layer_acr, 'MC')
        xlabel('Time With Respect to Inhalation-Onset [s]', 'FontSize', axis_size);
        sgtitle('Daughter Mitral Cells', 'FontSize', axis_size);
    elseif strcmp(layer_acr, 'Glom')
        sgtitle('Selected Glomerulus', 'FontSize', axis_size);
    end

    % Add legend
    legend(handles, labels, 'Location', 'northeast', 'Orientation', 'vertical', 'EdgeColor', 'none');

    % Save figure
    print(gcf, fullfile(plot_dir, sprintf('%s_dF_%s_plot.png', layer_acr, chosen_odor{k})), '-dpng', '-r300');
    hold off;


%% FLUORESCENCE AMPLITUDE PLOTS FOR EACH ROI

    % Initialize the amplitude matrix
    ampl_chosen_odor = zeros(n_FilePaths, n_ROI);

    for FilePathIndex = 1:n_FilePaths

        % Organize data
        chosen_odor_response = ODOR_RESPONSE_MEAN_norm{FilePathIndex,1}{1,idxChosenOdor};
       
        % Compute the average of rows 60 to 75 (0.5 sec post inhilation)
        avg_response = mean(chosen_odor_response(pre_inh : pre_inh+15, :), 1);
        plot_empty = mean(ODOR_RESPONSE_MEAN_norm_empty_avg(pre_inh : pre_inh+15, :), 1);

        ampl_chosen_odor(FilePathIndex, :) = avg_response;  

    end

    %% CREATE MASKS MATRIX

    Masks = cell(n_FilePaths, n_ROI);

    for i = 1:n_FilePaths
        data = load(FilePaths{i}, 'Session');

        % Store the 512x512 matrices into the corresponding row of Masks
        for j = 1:n_ROI
             Masks{i, j} = reshape(data.Session.CellMask(:, j), [512, 512]);
        end

    end

    %% Attribute Fluorescence value to Masks matrix

    % Loop through the values in ampl_chosen_odor and replace the 1s in Masks
    for i = 1:n_FilePaths
        for j = 1:n_ROI
            if ~isempty(Masks{i, j})

                % Replace 1s in the matrix with the amplification value
                mask_matrix = zeros(size(Masks{i, j}));
                mask_matrix(Masks{i, j} ~= 0) = ampl_chosen_odor(i, j);

                % Store the updated matrix back in Masks
                Masks{i, j} = mask_matrix;
            end
        end
    end

    Masks_empty = Masks(1, :);

    for j = 1:n_ROI
        if ~isempty(Masks_empty{1, j})

            % Replace 1s in the matrix with the amplification value
            mask_matrix = zeros(size(Masks_empty{1, j}));
            mask_matrix(Masks{1, j} ~= 0) = plot_empty(j);

            % Store the updated matrix back in Masks
            Masks_empty{1, j} = mask_matrix;
        end
    end

    %% Add all Masks Together

    % Initialize a cell array to store the summed results
    Masks_sum = cell(n_FilePaths, 1);
    Masks_empty_sum = cell(1, 1);

    % Loop through each row of Masks
    for i = 1:n_FilePaths
        % Initialize a matrix with zeros of the same size as the first mask in the row
        sum_matrix = zeros(size(Masks{i, 1}));

        % Add each mask in the current row to the summed_matrix
        for j = 1:n_ROI
            sum_matrix = sum_matrix + Masks{i, j};
        end

        % Store the summed matrix in the summed_masks cell array
        Masks_sum{i} = sum_matrix;
    end

    sum_matrix_empty = zeros(size(Masks{i, 1}));

    % Add each mask in the current row to the summed_matrix
    for j = 1:n_ROI
        sum_matrix_empty = sum_matrix_empty + Masks_empty{1, j};
    end

    % Store the summed matrix in the summed_masks cell array
    Masks_empty_sum{1} = sum_matrix_empty;

   
    %% Normalize for Plotting for the specific odor in question

    % Find the overall maximum value from Masks_sum and Masks_empty_sum

    % Compute the maximum value across all Masks_sum matrices
    max_values_sum = cellfun(@(x) max(x(:)), Masks_sum);

    % Compute the maximum value from the Masks_empty_sum matrix
    max_value_empty = max(Masks_empty_sum{1}(:));

    % Combine max_values_sum and max_value_empty into a single vector
    all_max_values = [max_values_sum(:); max_value_empty];

    % Find the overall maximum value
    max_value_overall = max(all_max_values);


    %% Create colormap

    % Define the range of values
    valueRange = -1:0.01:1;

    % Define the start, middle, and end colors
    startColor = [0 0 1]; % Blue for minimum value
    midColor = [1 1 1];   % White for zero value
    endColor = [1 0 0];   % Red for maximum value

    % Generate the color scales
    colorScale = generateColorScale(startColor, midColor, endColor, valueRange);

    %% ROI FIGURE

    fig = figure;
    set(fig, 'Position', [20, 50, 900, 300]);

    % Create a tiled layout with two columns (one for each subplot)
    plt_layout = tiledlayout(1, size(Masks_sum,1) + 1, 'TileSpacing', 'none', 'Padding', 'none'); % +1 for the no-odor subplot
    plt_layout.TileSpacing = 'compact';
    plt_layout.Padding = 'compact';

    % Plot the no-stimulus condition
    data_matrix_empty = Masks_empty_sum{1,1}; % Normalize with the overall max for this odor
    color_indices_no_stim = arrayfun(@(x) findClosestIndex(x, valueRange), data_matrix_empty);
    color_image_no_stim = ind2rgb(color_indices_no_stim, colorScale);
   
    % Plot Masks_empty_sum in the first tile
    nexttile;
    imshow(color_image_no_stim);
    ylabel(sprintf('%s', layer_name),'FontSize', axis_size);
    if strcmp(layer_acr, 'MC')
        xlabel(sprintf('no-odor'),'FontSize', axis_size);
    end
    

    [~, sorted_indices] = sort(chosen_odor_conc);    

    for idx = 1: length(sorted_indices)
        i = sorted_indices(idx); % Get the actual index in the original order
   
        data_matrix = Masks_sum{i}; % / max_value_overall; % Normalize with the overall max for this odor
   
        % Convert the values to indices in the colorScale
        color_indices = arrayfun(@(x) findClosestIndex(x, valueRange), data_matrix);
   
        % Create an RGB image based on the colorScale
        color_image = ind2rgb(color_indices, colorScale);
   
        % Create subplots for each path
        nexttile;
        imshow(color_image);
   
        hold on;
        % Loop through each ROIindex and plot the contours
        for j = 1:length(ROIindex_single)
            logical_mask = Masks{i, ROIindex_single(j)} ~= 0;
            contour(logical_mask, [1, 1], 'LineColor', 'k', 'LineWidth', 1.00);
        end
   
        chosen_odor_conc_curr = chosen_odor_conc(i);
        
        if strcmp(layer_acr, 'MC')
        % Add x-axis label to each subplot
            xlabel(sprintf('%s Concentration %.0f%%', chosen_odor{k}, chosen_odor_conc_curr),'FontSize', axis_size);
        end

    end

    % Adjust colormap to use the generated color scale
    colormap(gca, colorScale);
    cb = colorbar('eastoutside');
    caxis([min(valueRange) max(valueRange)]);
    cb.Ticks = [-1, 0, 1];

    chosen_odor_conc_curr = chosen_odor_conc(i);

    % Save the figure
    print(gcf, fullfile(plot_dir, sprintf('%s_%s_mask_plot.png', layer_acr, chosen_odor{k})), '-dpng', '-r300');

end



%% Find the average fluorescence 15 frames (0.5sec) post stim presentation

stim_time = [10, 30, 60, 120, 180, 240];
n_latencies = length(stim_time);

% stim_time to frame_time
stim_frame = (stim_time * 30 / 1000);
stim_frame = pre_inh + round(stim_frame);
stim_frame_end = stim_frame + 15;
stim_frame_end = round(stim_frame_end);

n_conc = size(ODOR_RESPONSE_MEAN_norm, 1); 

% Initialize ODOR_ONLY_mean
ODOR_ONLY_mean = cell(n_conc, 1);

% Loop through concentrations
for idx_conc = 1:n_conc
    % Initialize the double for this concentration
    ODOR_ONLY_mean{idx_conc, 1} = zeros(length(ROIindex_single), length(stim_time));
    
    % Loop through ROI indices
    for idx_ROI = 1:length(ROIindex_single)
        % Extract the column corresponding to the current ROI index
        selected_column = ODOR_RESPONSE_MEAN_norm{idx_conc, 1}{1, chosen_odor_idx}(:, ROIindex_single(idx_ROI));
        
        % Loop through stim_time
        for idx_stim_time = 1:length(stim_time)
            row_start = stim_frame(idx_stim_time);
            row_end = stim_frame_end(idx_stim_time);

            % Calculate the mean value in the selected column over the specified frames
            mean_value = mean(selected_column(row_start:row_end));
            
            % Store the mean value in the appropriate location in ODOR_ONLY_ROI_response
            ODOR_ONLY_mean{idx_conc, 1}(idx_ROI, idx_stim_time) = mean_value;
        
        end
    end
end

ODOR_ONLY_mean_high = ODOR_ONLY_mean{1,:}(:,:)
ODOR_ONLY_mean_high_avg = mean(ODOR_ONLY_mean_high, 1);  % average for all ROI's

ODOR_ONLY_mean_low = ODOR_ONLY_mean{2,:}(:,:)
ODOR_ONLY_mean_low_avg = mean(ODOR_ONLY_mean_low, 1); % average for all ROI's

% Initialize arrays to hold standard deviations and SEMs
ODOR_ONLY_std = cell(n_conc, 1);
ODOR_ONLY_sem = cell(n_conc, 1);

for idx_conc = 1:n_conc
    % Preallocate space for standard deviations
    ODOR_ONLY_std{idx_conc, 1} = zeros(size(ODOR_ONLY_mean{idx_conc, 1}));
    
    % Loop through ROIs and latencies
    for i = 1:length(ROIindex_single)
        for j = 1:n_latencies
            % Get the specific ROI index
            r = ROIindex_single(i);
            
            % Calculate standard deviation over the frame range
            frame_std = std(ODOR_RESPONSE_MEAN_norm{idx_conc, 1}{1, chosen_odor_idx}(stim_frame(j):stim_frame_end(j), r), 0, 1);
            
            % Store STD in the preallocated matrix (use parentheses for matrix operations)
            ODOR_ONLY_std{idx_conc, 1}(i,j) = frame_std;
        end
    end

    % Calculate SEM (using n_trials as the divisor)
    ODOR_ONLY_sem{idx_conc, 1} = ODOR_ONLY_std{idx_conc, 1} / sqrt(n_trials);
end

ODOR_ONLY_sem_high = ODOR_ONLY_sem{1,:}(:,:);
ODOR_ONLY_sem_high_avg = mean(ODOR_ONLY_sem_high, 1); % average for all ROI's

ODOR_ONLY_sem_low = ODOR_ONLY_sem{2,:}(:,:);
ODOR_ONLY_sem_low_avg = mean(ODOR_ONLY_sem_low, 1); % average for all ROI's

% Average accross the neccessary ROIindexes and the time-frames in question
% Compute this average on ODOR_RESPONSE_MEAN_norm_empty_avg
% ODOR_RESPONSE_MEAN_norm_empty_avg : empty trials averaged accross both concentrations
% Initialize ODOR_ONLY_mean_empty
ODOR_ONLY_empty_mean = zeros(1, n_latencies);  % n_latencies = length(stim_time)

% Loop through each latency (stim_time) and calculate the average between stim_frame and stim_frame_end
for i = 1:n_latencies
    % Extract the data between stim_frame(i) and stim_frame_end(i) for the selected ROIs (columns from ROIindex_single)
    selected_data = ODOR_RESPONSE_MEAN_norm_empty_avg(stim_frame(i):stim_frame_end(i), ROIindex_single);
    
    % Calculate the mean of the selected rows and columns, then store in ODOR_ONLY_mean_empty
    ODOR_ONLY_empty_mean(i) = mean(selected_data(:));  % Averages all values within the range across the selected ROIs
end

ODOR_ONLY_empty_std = zeros(length(ROIindex_single), n_latencies);  % STD matrix
ODOR_ONLY_empty_sem = zeros(length(ROIindex_single), n_latencies);  % SEM matrix
    
% Loop through ROIs and latencies
for i = 1:length(ROIindex_single)
    for j = 1:n_latencies
        % Get the specific ROI index
        r = ROIindex_single(i);

        % Calculate standard deviation over the frame range
        frame_std = std(ODOR_RESPONSE_MEAN_norm_empty_avg(stim_frame(j):stim_frame_end(j), r), 0, 1);

        % Store STD in the matrix
        ODOR_ONLY_empty_std(i, j) = frame_std;
    end
end

% Calculate SEM (using n_trials as the divisor)
ODOR_ONLY_empty_sem = ODOR_ONLY_empty_std / sqrt(n_trials);  % SEM matrix

% Step to average the standard deviation and SEM across their rows/
% ROIs bc we only want to plot a proxy for the 'no stim nor odor' curve

% Average standard deviations across the rows (i.e., across ROIs) for each latency
ODOR_ONLY_empty_std_avg = mean(ODOR_ONLY_empty_std, 1);  % Averages along rows (dimension 1)

% Average SEM across the rows (i.e., across ROIs) for each latency
ODOR_ONLY_empty_sem_avg = mean(ODOR_ONLY_empty_sem, 1);  % Averages along rows (dimension 1)


%% PLOT STIM-TIME AVERAGE PLOT - AVERAGE ACCROSS ROI's

% no stim or odor (black)
% stim (red)
% low BENZ (blue) (0.5 alpha)
% high BENZ (blue) (1.0 alpha)
% stim + low BENZ (purple) (0.5 alpha)
% stim + high BENZ (purple) (1.0 alpha)

% Initialize the figure
figure_handle = figure('Name', 'ODOR AVG Plot - average accross ROIs');

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 400]);

% Hold on to plot multiple lines
hold on;

% LEGEND CUSTOM LABEL
if strcmp(layer_acr, 'MC')
    custom_labels = {
        sprintf('high %s, daughter MCs', chosen_odor{k}), [0 0 1];  % Blue (1.0 alpha)
        sprintf('low %s, daughter MCs', chosen_odor{k}), [0.5 0.5 1];  % Blue (0.5 alpha)
        'no stim nor odor', [0 0 0];                                  % Black for no stim or odor
    };
elseif strcmp(layer_acr, 'Glom')
    custom_labels = {     
        sprintf('high %s, single Glom', chosen_odor{k}), [0 0 1];  % Blue (1.0 alpha)
        sprintf('low %s, single Glom', chosen_odor{k}), [0.5 0.5 1];   % Blue (0.5 alpha)
        'no stim nor odor', [0 0 0];                                   % Black for no stim or odor
    };
end

% NO STIM NOR ODOR: plot average data and store the handle
h1 = plot(stim_time, ODOR_ONLY_empty_mean, 'Color', custom_labels{3, 2}, ...
          'LineWidth', 1, 'LineStyle', '-');
      
% Plot shaded area for SEM (uncertainty) for 'no stim nor odor'
fill([stim_time, fliplr(stim_time)], ...
     [ODOR_ONLY_empty_mean - ODOR_ONLY_empty_sem_avg, fliplr(ODOR_ONLY_empty_mean + ODOR_ONLY_empty_sem_avg)], ...
     custom_labels{3, 2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');   

% ODOR - HIGH CONC.: plot average and SEM area
h2 = plot(stim_time, ODOR_ONLY_mean_high_avg, 'Color', custom_labels{1, 2}, ...
          'LineWidth', 1, 'LineStyle', '-');
fill([stim_time, fliplr(stim_time)], ...
     [ODOR_ONLY_mean_high_avg - ODOR_ONLY_sem_high_avg, fliplr(ODOR_ONLY_mean_high_avg + ODOR_ONLY_sem_high_avg)], ...
     custom_labels{1, 2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% ODOR - LOW CONC.: plot average and SEM area
h3 = plot(stim_time, ODOR_ONLY_mean_low_avg, 'Color', custom_labels{2, 2}, ...
          'LineWidth', 1, 'LineStyle', '-');
fill([stim_time, fliplr(stim_time)], ...
     [ODOR_ONLY_mean_low_avg - ODOR_ONLY_sem_low_avg, fliplr(ODOR_ONLY_mean_low_avg + ODOR_ONLY_sem_low_avg)], ...
     custom_labels{2, 2}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

% Adjust y-axis limits as needed
ylim(avg_y_lim);

% Set x-axis ticks
xticks(stim_time);

% Add labels and title
if strcmp(layer_acr, 'MC')
    xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
end
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex', 'FontSize', axis_size);

if strcmp(layer_acr, 'MC')
    sgtitle('Daughter Mitral Cells', 'FontSize', axis_size);
elseif strcmp(layer_acr, 'Glom')
    sgtitle('Selected Glomerulus', 'FontSize', axis_size);
end

% Create legend with all relevant handles
legend([h2, h3, h1], custom_labels{:, 1}, 'Location', 'northeast', 'EdgeColor', 'none');

% Save the figure in high resolution (e.g., 300 DPI)
print(gcf, fullfile(plot_dir, sprintf('%s_%s_ODOR_plot_avg.png', layer_acr, chosen_odor{k})), '-dpng', '-r300');

% Finish plotting
hold off;

%% PLOT STIM-TIME AVERAGE PLOT 

% no stim or odor (black)
% stim (red)
% low BENZ (blue) (0.5 alpha)
% high BENZ (blue) (1.0 alpha)
% stim + low BENZ (purple) (0.5 alpha)
% stim + high BENZ (purple) (1.0 alpha)

% Initialize the figure
figure_handle = figure('Name', 'ODOR AVG Plot');

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 400]);

% Hold on to plot multiple lines
hold on;

% LEGEND CUSTOM LABEL
if strcmp(layer_acr, 'MC')
    custom_labels = {
        sprintf('high %s, daughter MCs', chosen_odor{k}), [0 0 1 1];  % Blue (1.0 alpha)
        sprintf('low %s, daughter MCs', chosen_odor{k}), [0.5 0.5 1 0.5];  % Blue (0.5 alpha)
        'no stim nor odor', [0 0 0];                                  % Black for no stim or odor
    };
elseif strcmp(layer_acr, 'Glom')
    custom_labels = {     
        sprintf('high %s, single Glom', chosen_odor{k}), [0 0 1 1];  % Blue (1.0 alpha)
        sprintf('low %s, single Glom', chosen_odor{k}), [0.5 0.5 1 0.5];   % Blue (0.5 alpha)
        'no stim nor odor', [0 0 0];                                   % Black for no stim or odor
    };
end

% ODOR - HIGH CONC.: plot each row of STIM_ONLY_numeric and fill the SEM area
for i = 1:size(ODOR_ONLY_mean_high, 1)
    % Plot each ROI line and store the handle
    h1(i) = plot(stim_time, ODOR_ONLY_mean_high(i, :), 'Color', custom_labels{1, 2}(1:3), ...
                 'LineWidth', 1, 'LineStyle', '-');
    
    % Plot shaded area for SEM (uncertainty) for the current ROI
    fill([stim_time, fliplr(stim_time)], ...
         [ODOR_ONLY_mean_high(i, :) - ODOR_ONLY_sem_high(i, :), fliplr(ODOR_ONLY_mean_high(i, :) + ODOR_ONLY_sem_high(i, :))], ...
         custom_labels{1, 2}(1:3), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% ODOR - LOW CONC.: plot each row of STIM_ONLY_numeric and fill the SEM area
for i = 1:size(ODOR_ONLY_mean_low, 1)
    % Plot each ROI line and store the handle
    h2(i) = plot(stim_time, ODOR_ONLY_mean_low(i, :), 'Color', custom_labels{2, 2}(1:3), ...
                 'LineWidth', 1, 'LineStyle', '-');
    
    % Plot shaded area for SEM (uncertainty) for the current ROI
    fill([stim_time, fliplr(stim_time)], ...
         [ODOR_ONLY_mean_low(i, :) - ODOR_ONLY_sem_low(i, :), fliplr(ODOR_ONLY_mean_low(i, :) + ODOR_ONLY_sem_low(i, :))], ...
         custom_labels{2, 2}(1:3), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

% NO STIM NOR ODOR: plot average data and store the handle
h3 = plot(stim_time, ODOR_ONLY_empty_mean, 'Color', custom_labels{3, 2}(1:3), ...
          'LineWidth', 1, 'LineStyle', '-');
      
% Plot shaded area for SEM (uncertainty) for the current ROI
fill([stim_time, fliplr(stim_time)], ...
     [ODOR_ONLY_empty_mean - ODOR_ONLY_empty_sem_avg, fliplr(ODOR_ONLY_empty_mean + ODOR_ONLY_empty_sem_avg)], ...
     custom_labels{3, 2}(1:3), 'FaceAlpha', 0.1, 'EdgeColor', 'none');   


% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

% Adjust y-axis limits as needed
ylim(avg_y_lim);

% Set x-axis ticks
xticks(stim_time);

% Add labels and title
xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex', 'FontSize', axis_size);

if strcmp(layer_acr, 'MC')
    sgtitle('Daughter Mitral Cells', 'FontSize', axis_size);
elseif strcmp(layer_acr, 'Glom')
    sgtitle('Selected Glomerulus', 'FontSize', axis_size);
end

% Create legend with all relevant handles
legend([h1(1), h2(1), h3], custom_labels{:, 1}, 'Location', 'northeast', 'EdgeColor', 'none');

% Save the figure in high resolution (e.g., 300 DPI)
print(gcf, fullfile(plot_dir, sprintf('%s_%s_ODOR_plot.png', layer_acr, chosen_odor{k})), '-dpng', '-r300');

% Finish plotting
hold off;


%% COMBINED PLOTS - dF

% Extract the odor from the cell array
odor = chosen_odor{1};

% Define the filenames for the two PNG files
file1 = sprintf('%s_dF_%s_plot.png', 'Glom', odor);
file2 = sprintf('%s_dF_%s_plot.png', 'MC', odor);

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
title(sprintf('Passive Odor Presentation | %s', chosen_odor{1}), 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 100, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('odor_ODOR-ONLY_dF_%s.png', chosen_odor{k})), '-dpng', '-r300');


%% COMBINED PLOTS - ROI Masks

% Extract the odor from the cell array
odor = chosen_odor{1};

% Define the filenames for the two PNG files
file1 = sprintf('%s_%s_mask_plot.png', 'Glom', odor);
file2 = sprintf('%s_%s_mask_plot.png', 'MC', odor);

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

title(sprintf('Passive Odor Presentation |  %s', chosen_odor{k}'), 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 100, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 150, sprintf('post-inhilation onset (0.5 sec) normalized fluorescence avg.'), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure in high resolution (e.g., 300 DPI)
print(gcf, fullfile(plot_dir, sprintf('odor_ODOR-ONLY_ROIs_%s.png', chosen_odor{k})), '-dpng', '-r300');


%% COMBINED STIM-TIME AVERAGE PLOT

% Define the filenames for the two PNG files
file1 = sprintf('Glom_%s_ODOR_plot_avg.png', chosen_odor{k});
file2 = sprintf('MC_%s_ODOR_plot_avg.png', chosen_odor{k});

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Create a blank space for text (adjust the height as needed)
% Create a blank space for text (adjust the height as needed)
text_height = 100; % Height for the text area
blank_space = uint8(255 * ones(text_height, max(size(img1, 2), size(img2, 2)), 3)); % White space

% Combine the images with the blank space
combined_img = [img1; img2; blank_space]; % Place the blank space at the bottom

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
title(sprintf('Passive Odor Presentation | %s', chosen_odor{1}), 'FontSize', title_size);

% Calculate the vertical position for text placement
baseY = size(combined_img, 1) - text_height / 2;

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), baseY + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), baseY + 100, 'normalized with respect to OB layer', ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');


% Save the figure
print(gcf, fullfile(plot_dir, sprintf('odor_ODOR-ONLY_stim-time_avg_%s.png', chosen_odor{k})), '-dpng', '-r300');

%% COLORSCALE FUNCTION

function colorScale = generateColorScale(startColor, midColor, endColor, valueRange)
    % Generate a color scale between startColor, midColor, and endColor for the specified valueRange
    midIndex = find(valueRange == 0, 1);
    numSteps = length(valueRange);
    colorScale = zeros(numSteps, 3);
    for i = 1:3
        colorScale(1:midIndex, i) = linspace(startColor(i), midColor(i), midIndex);
        colorScale(midIndex:end, i) = linspace(midColor(i), endColor(i), numSteps - midIndex + 1);
    end
end

%% FIND CLOSEST INDEX FUNCTION

function idx = findClosestIndex(value, range)
    [~, idx] = min(abs(range - value));
end