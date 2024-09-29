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

axis_size = 12;
title_size = 14;

%% CHOSEN ODOR TO DISPLAY

chosen_odor = {'Heptanal'};
% chosen_odor = {'EthylTiglate','EthylButyrate','Hexanal','BenzAldehyde','MethylValerate'};

% Extended Odor Panel
% chosen_odor = {'2MBA', '5M2H','EthylButyrate','Hexanal','BenzAldehyde','Heptanal','MethylValerate','PropionicAcid'};
% chosen_odor = {'2MBA','Hexanal','BenzAldehyde','Heptanal'}; % chosen from 240731 session

% chosen_odor = {'Heptanal','Hexanal'}; 

%% DEFINE FILE PATHS

% FilePaths = {'mouse0773/240703_Glom_100flow/aligned/240703_Glom_OdorMapping_Z0_1X_100flow_S_v73.mat',
%              'mouse0773/240705_Glom_20flow/aligned/240705_Glom_OdorMapping_Z0_1X_20flow_S_v73.mat'};
% ROIindex_single = [10];
% ROIindex = [2,5,4,7,10]; % BenzAldehyde
% % ROIindex = [10,13,16]; % EthylButyrate
% ROIindex = [10,11,12,14,15,17,18]; % Hexanal
% ROIindex = [7,9,10,11,12,13,16]; % MethylValerate
% % ROIindex = [2,7,9,10,18]; % EthylTiglate - not great

% FilePaths = {'mouse0773/240703_MC_100flow/aligned/240703_MC_OdorMapping_Z196_3X_100flow_S_v73.mat',
%               'mouse0773/240705_MC_20flow/aligned/240703_MC_OdorMapping_Z196_3X_20flow_S_v73.mat'};
% ROIindex = [4,7,9];
% ROIindex_single = [4,7,9];

% FilePaths = {'mouse0781/240705_Glom_100flow/aligned/240705_Glom_OdorMapping_Z0_1X_100flow_S_v73.mat',
%             'mouse0781/240705_Glom_20flow/aligned/240705_Glom_OdorMapping_Z0_1X_20flow_S_v73.mat'};
% ROIindex_single = [13];
% ROIindex = [12,13]; % BenzAldehyde
% ROIindex = [8,9,12,13]; % EthylButyrate
% ROIindex = [11,12,13,15,16,22,23]; % Hexanal
% ROIindex = [11,12,13,16]; % MethylValerate
% ROIindex = [12,13,24]; % EthylTiglate

% FilePaths = {'mouse0781/240731_Glom_100flow/aligned/240731_Glom_OdorMapping_Z0_1X_100flow_S_v73.mat',
%              'mouse0781/240731_Glom_20flow/aligned/240731_0781_Glom_OdorMapping_Z0_1X_20flow_S_v73.mat'};
% ROIindex_single = [15];
% ROIindex = [7,9,8,10,11,15]; % 2MBA
% ROIindex = [7,11,15,18,20]; % Hexanal
% ROIindex = [7,10,11,15]; % BenzAldehyde
% ROIindex = [7,10,11,15,21,28]; % Heptanal
% glom 12 is highly sensitive to the sniff cycle and could be removed from the cross-trace plots

FilePaths = {'mouse0781/240801_MC_100flow/aligned/240801_0781_MC_OdorMapping_Z94_3X_100flow_S_v73.mat',
             'mouse0781/240801_MC_20flow/aligned/240801_0781_MC_OdorMapping_Z94_3X_20flow_S_v73.mat'};
ROIindex_single = [8,10,18,26];
ROIindex = [8,10,18,26];


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



%% PRINCIPLE VARIABLES

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

% CHECK NORMALIZATION MATRIX
% target_value = 1;  % Replace with the value you want to search for
%
% % Initialize a cell array to store the locations where the value is found
% value_locations = {};
%
% % Loop through the cell array to find the target value
% for row = 1:size(ODOR_RESPONSE_MEAN_norm, 1)
%     for col = 1:size(ODOR_RESPONSE_MEAN_norm, 2)
%         % Access the inner cell array
%         inner_cell_array = ODOR_RESPONSE_MEAN_norm{row, col};
%
%         % Loop through each matrix within the inner cell array
%         for inner_idx = 1:length(inner_cell_array)
%             % Access the matrix
%             current_matrix = inner_cell_array{inner_idx};
%
%             % Find the indices where the target value appears
%             [row_idx, col_idx] = find(current_matrix == target_value);
%
%             % If the value is found, store the location details
%             if ~isempty(row_idx)
%                 % Store the details in a structured format
%                 value_locations{end+1} = struct( ...
%                     'OuterRow', row, ...
%                     'OuterCol', col, ...
%                     'InnerIndex', inner_idx, ...
%                     'RowIndices', row_idx, ...
%                     'ColIndices', col_idx ...
%                 );
%             end
%         end
%     end
% end
%
% % Display the results
% if isempty(value_locations)
%     disp('The target value was not found in ODOR_RESPONSE_MEAN_norm.');
% else
%     disp('The target value was found at the following locations:');
%     for k = 1:length(value_locations)
%         loc = value_locations{k};
%         fprintf('Outer cell [%d, %d], Inner cell %d:\n', loc.OuterRow, loc.OuterCol, loc.InnerIndex);
%         for idx = 1:length(loc.RowIndices)
%             fprintf('  Row: %d, Column: %d\n', loc.RowIndices(idx), loc.ColIndices(idx));
%         end
%     end
% end


%% NO-ODOR CODE

ODOR_RESPONSE_MEAN_norm_empty = cell(2,1);

for cell_array_index = 1:n_FilePaths
    ODOR_RESPONSE_MEAN_norm_empty{cell_array_index,1} = ODOR_RESPONSE_MEAN_norm{cell_array_index,1}{1,6};
end

plot_no_odor = cell(2, 1); % Initialize plot_no_odor as a 2x1 cell array

for cell_array_index = 1:size(ODOR_RESPONSE_MEAN_norm_empty, 1)
    % Extract the relevant data matrix
    data_matrix = ODOR_RESPONSE_MEAN_norm_empty{cell_array_index};

    % Extract columns corresponding to ROIindex
    selected_columns = data_matrix(:, 1:n_ROI);

    % Average the selected columns along the second dimension (columns)
    averaged_columns = mean(selected_columns, 2); % Calculate mean along columns

    % Store the result in plot_no_odor
    plot_no_odor{cell_array_index} = averaged_columns;
end

% plot_no_odor averaged
plot_no_odor_avg = mean(cell2mat(plot_no_odor'), 2);  % Averages along each corresponding element of the two doubles

%% LOOP THROUGH CHOSEN ODORS FOR PLOTS

% Loop through the chosen odor indices
for k = 1:numel(chosen_odor_idx)
    idxChosenOdor = chosen_odor_idx(k);
    chosen_odor_conc = concentrations(:, idxChosenOdor);

%% TiME CROSS-SECTION PLOT

% Define a colormap to ensure distinct colors for each ROI
colors = lines(numel(ROIindex));

% Sort the concentrations and get the sorting indices in descending order
[sorted_conc, sort_idx] = sort(chosen_odor_conc, 'descend');

% Create a new figure for each chosen odor index
figure_handle = figure('Name', ['Chosen Odor ', num2str(idxChosenOdor)]);
% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 50, 600, 400]);

% Initialize flags to track if specific legend labels have been added
added_single_high = false;
added_single_low = false;
added_neighbor_high = false;
added_neighbor_low = false;

% Initialize empty arrays for plot handles and labels
handles = [];
labels = [];

for sorted_idx = 1:length(sort_idx)
    % Get the index from the sorted indices
    cell_array_index = sort_idx(sorted_idx);
    chosen_odor_conc_curr = sorted_conc(sorted_idx);
    OdorResponse_mean = ODOR_RESPONSE_MEAN_norm{cell_array_index, 1}{1, idxChosenOdor};
    
    % Hold on to plot on the same figure
    hold on;

    % Loop through each ROI index and plot the corresponding columns
    for j = 1:numel(ROIindex)
        roi_index = ROIindex(j);
        currentColumn = OdorResponse_mean(:, roi_index);
        
        % Check if the ROI is in ROIindex_single
        if ismember(roi_index, ROIindex_single)
            color = 'r'; % Red for matching ROIs
            if sorted_idx == 2  % Low concentration
                if strcmp(layer_acr, 'MC')
                    legend_label = 'daughter MC - low conc.';
                else
                    legend_label = 'single Glom - low conc.';
                end
                if ~added_single_low
                    h_single_low = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', [0, 0, 1, 0.4], 'LineStyle', '-');
                    added_single_low = true;
                    handles(end+1) = h_single_low; %#ok<*AGROW>
                    labels{end+1} = legend_label;
                else
                    plot(tt, currentColumn, 'LineWidth', 1.5, 'Color', [0, 0, 1, 0.4], 'LineStyle', '-', 'HandleVisibility', 'off');
                end
            else  % High concentration
                if strcmp(layer_acr, 'MC')
                    legend_label = 'daughter MC - high conc.';
                else
                    legend_label = 'single Glom - high conc.';
                end
                if ~added_single_high
                    h_single_high = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color',  [0, 0, 1, 1]);
                    added_single_high = true;
                    handles(end+1) = h_single_high;
                    labels{end+1} = legend_label;
                else
                    plot(tt, currentColumn, 'LineWidth', 1.5, 'Color',  [0, 0, 1, 1], 'HandleVisibility', 'off');
                end
            end
        else
            color = [0.5, 0.5, 0.5]; % Grey for non-matching ROIs
            if sorted_idx == 2  % Low concentration
                if strcmp(layer_acr, 'Glom')
                    legend_label = 'neighbor Gloms - low conc.';
                else
                    legend_label = [layer_acr, ' ', num2str(roi_index)];
                end
                if ~added_neighbor_low
                    h_neighbor_low = plot(tt, currentColumn, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', [0.5, 0.5, 0.5, 0.3]);
                    added_neighbor_low = true;
                    handles(end+1) = h_neighbor_low;
                    labels{end+1} = legend_label;
                else
                    plot(tt, currentColumn, 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5, 0.3], 'HandleVisibility', 'off');
                end
            else  % High concentration
                if strcmp(layer_acr, 'Glom')
                    legend_label = 'neighbor Gloms - high conc.';
                else
                    legend_label = [layer_acr, ' ', num2str(roi_index)];
                end
                if ~added_neighbor_high
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
end

% Plot a single 'no-odor' line and make sure it is black
h_no_odor = plot(tt, plot_no_odor_avg, 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'no-odor');
handles(end+1) = h_no_odor;
labels{end+1} = 'no-odor';

% Set plot limits and labels
xlim([-0.5, 1.5]);
ylim([-0.3, 0.5]);

% Dashed vertical line at x = 0
line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '-', 'LineWidth', 1, 'HandleVisibility', 'off');
% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.2], 'LineStyle', '-', 'HandleVisibility', 'off');

% Add labels, title, etc. as needed
ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}','Interpreter','latex','FontName','Arial','FontSize', axis_size);

% Conditionally add x-axis label based on layer_acr
if strcmp(layer_acr, 'MC')
    xlabel('Time With Respect to Inhilation-Onset [s]','FontSize', axis_size);
end

% Add legend in the desired order
legend(handles, labels, 'Location', 'northwest', 'Orientation', 'vertical', 'EdgeColor', 'none');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('%s_dF_%s_plot.png', layer_acr, chosen_odor{k})), '-dpng', '-r300');

hold off;


%% FLUORESCENCE AMPLITUDE PLOTS FOR EACH ROI

ODOR_RESPONSE_MEAN_norm_empty_avg = (ODOR_RESPONSE_MEAN_norm_empty{1} + ODOR_RESPONSE_MEAN_norm_empty{2}) / 2;

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
text(size(combined_img, 2), size(combined_img, 1) + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('%s_combined_dF_plot.png', chosen_odor{k})), '-dpng', '-r300');

%% COMBINED PLOTS - Mask

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
text(size(combined_img, 2), size(combined_img, 1) + 50, sprintf('Post-Inhilation Onset (0.5 sec) Normalized Fluorescence Avg.'), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure in high resolution (e.g., 300 DPI)
print(gcf, fullfile(plot_dir, sprintf('%s_combined_mask_plot.png', chosen_odor{k})), '-dpng', '-r300');



%% Find the average fluorescence 15 frames (0.5sec) post stim presentation

stim_time = [10, 30, 60, 120, 180, 240];

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

% Initialize ODOR_ONLY_mean_empty
ODOR_ONLY_mean_empty = zeros(1, length(stim_time));

% STEP 1: Average the two doubles in ODOR_RESPONSE_MEAN_norm_empty
averaged_double = mean(cat(3, ODOR_RESPONSE_MEAN_norm_empty{1}, ODOR_RESPONSE_MEAN_norm_empty{2}), 3);

% STEP 2: Extract and average the columns corresponding to ROIindex_single
selected_columns_avg = mean(averaged_double(:, ROIindex_single), 2);

% STEP 3: Loop through stim_time to calculate and store the mean values
stim_frame = pre_inh + round(stim_time * 30 / 1000);
stim_frame_end = round(stim_frame + 15);

for idx_stim_time = 1:length(stim_time)
    row_start = stim_frame(idx_stim_time);
    row_end = stim_frame_end(idx_stim_time);
    
    % Calculate the mean value in the averaged column over the specified frames
    mean_value = mean(selected_columns_avg(row_start:row_end));
    
    % Store the mean value in the appropriate location in ODOR_ONLY_mean_empty
    ODOR_ONLY_mean_empty(1, idx_stim_time) = mean_value;
end

%% PLOT STIM-TIME AVERAGE PLOT

% no stim or odor (black)
% stim (red)
% low BENZ (blue) (0.5 alpha)
% high BENZ (blue) (1.0 alpha)
% stim + low BENZ (purple) (0.5 alpha)
% stim + high BENZ (purple) (1.0 alpha)

% Initialize the figure
figure_handle = figure('Name', 'ODOR+STIM Plot');

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 400]);

% Hold on to plot multiple lines
hold on;

% LEGEND CUSTOM LABEL
if strcmp(layer_acr, 'MC')
    custom_labels = {
        'no stim or odor',                   [0 0 0 1];           % Black for no stim or odor
        sprintf('low %s, daughter MC', chosen_odor{k}),           [0 0 1 0.4];         % (blue, 0.4 alpha)
        sprintf('high %s, daughter MC', chosen_odor{k}),         [0 0 1 1];           % (blue, 1.0 alpha)
    };
elseif strcmp(layer_acr, 'Glom')
    custom_labels = {
        'no stim or odor',                   [0 0 0 1];           % Black for no stim or odor
        sprintf('low %s, single Glom', chosen_odor{k}),          [0 0 1 0.4];         % (blue, 0.4 alpha)
        sprintf('high %s, single Glom', chosen_odor{k}),         [0 0 1 1];           % (blue, 1.0 alpha)
    };
end

% Plot the 'no stim or odor' data from ODOR_ONLY_mean_empty
plot(stim_time, ODOR_ONLY_mean_empty, 'Color', custom_labels{1, 2}, ...
     'LineWidth', 1, 'LineStyle', '-');

% Plot the data with different colors and transparencies
for idx_conc = 1:n_conc
    for i = 1:length(ROIindex_single)
        data_to_plot = ODOR_ONLY_mean{idx_conc, 1}(i, :);
        
        % Determine the color based on the concentration
        if idx_conc == 1
            color = custom_labels{2, 2}; % low BENZ (blue, 0.4 alpha)
        elseif idx_conc == 2
            color = custom_labels{3, 2}; % high BENZ (blue, 1.0 alpha)
        end
        
        % Plot the data
        plot(stim_time, data_to_plot, 'Color', color, ...
             'LineWidth', 1, 'LineStyle', '-');
    end
end

% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

% Adjust y-axis limits as needed
ylim([-0.3, 0.3]);

% Set x-axis ticks
xticks(stim_time);

% Add labels and title
if strcmp(layer_acr, 'MC')
    xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
end
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex');

% Create legend
legend(custom_labels{:, 1}, 'Location', 'northeast', 'EdgeColor', 'none');

% Save the figure in high resolution (e.g., 300 DPI)
print(gcf, fullfile(plot_dir, sprintf('%s_%s_ODOR-STIM_plot.png', layer_acr, chosen_odor{k})), '-dpng', '-r300');

% Finish plotting
hold off;

%% COMBINED STIM-TIME AVERAGE PLOT

% Define the filenames for the two PNG files
file1 = sprintf('Glom_%s_ODOR-STIM_plot.png', chosen_odor{k});
file2 = sprintf('MC_%s_ODOR-STIM_plot.png', chosen_odor{k});

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
title(sprintf('Passive Odor Presentation | %s', chosen_odor{k}), 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure in high resolution (e.g., 300 DPI)
print(gcf, fullfile(plot_dir, sprintf('%s_ODOR-STIM_plot.png', chosen_odor{k})), '-dpng', '-r300');


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