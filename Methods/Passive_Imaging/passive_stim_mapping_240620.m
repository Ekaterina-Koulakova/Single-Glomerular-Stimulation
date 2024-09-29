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

axis_size = 12;
title_size = 14;

%% DEFINE FILE PATHS

% FilePaths = {'mouse0773/240702_Glom_Stim/aligned/240702_0773_StimMapping_Glom_Z0_1X_S_v73.mat'};
% ROIindex = [16];

% FilePaths = {'mouse0773/240702_MC_Stim/aligned/240702_0773_StimMapping_MC_Z196_3X_S_v73.mat'};
% ROIindex = [3,7,10]; % MC exc and inh highly dependent on latency

% FilePaths = {'mouse0781/240528_Glom_Stim/aligned/240528_Glom_StimMapping_S_v73.mat'};
% ROIindex = [12];

FilePaths = {'mouse0781/240528_MC_Stim/aligned/240528_MC_StimMapping_S_v73.mat'};
ROIindex = [5,10,22,29];


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

%% DEFINE KEY VARIABLES

% Load the data from the first file in FilePaths
data = load(FilePaths{1}, 'Session');

% TIME VARIABLES
pre_inh = data.Session.Infos.pre_inh;
post_inh = data.Session.Infos.post_inh;
n_frames = pre_inh + post_inh;
n_fps = data.Session.Infos.fps;
tt = linspace(-(pre_inh/n_fps), (post_inh/n_fps), n_frames);

n_ROI = size(data.Session.OdorResponse{1,1}, 2);
n_trials = size(data.Session.OdorResponse{1,1}, 3);


%% ORGANIZING DATA

% Initialize the concentrations array as empty
stimpresent = []; 

% Initialize cell arrays for trial_stim_type and trial_latencies
numConds = length(data.Session.UniqueConds);
trial_stim_type = cell(numConds, 1);
trial_latencies = cell(numConds, 1);
    
 % Extract the required parts from UniqueConds
for i = 1:numConds
    cond = data.Session.UniqueConds{i};
    % Extract stim type
    stim_type_match = regexp(cond, 'id:(.*?)(/L:|$)', 'tokens');
    if ~isempty(stim_type_match)
        trial_stim_type{i} = stim_type_match{1}{1};
    else
        trial_stim_type{i} = ''; % If no match, leave empty
    end
    
    % Extract latency
    latency_match = regexp(cond, '/L:(\d+)', 'tokens');
    if ~isempty(latency_match)
        trial_latencies{i} = latency_match{1}{1};
    else
        trial_latencies{i} = ''; % If no match, leave empty
    end
end

% ORDER DATA BY LATENCIES

trial_lat_ROI_1 = {};
trial_lat_ROI_2 = {};

% Loop through each element in trial_stim_type
for i = 1:length(trial_stim_type)
    if strcmp(trial_stim_type{i}, 'ROI_1.tif')
        trial_lat_ROI_1{1, end+1} = i;            % Store the index
        trial_lat_ROI_1{2, end} = trial_latencies{i};  % Store the value
    end
    
    if strcmp(trial_stim_type{i}, 'ROI_2.tif')
        trial_lat_ROI_2{1, end+1} = i;            % Store the index
        trial_lat_ROI_2{2, end} = trial_latencies{i};  % Store the value
    end
end

% Sort trial_lat_ROI_1 by the values in the 2nd row
if ~isempty(trial_lat_ROI_1)
    [~, sort_idx_ROI_1] = sort(cellfun(@str2double, trial_lat_ROI_1(2, :)));
    trial_lat_ROI_1 = trial_lat_ROI_1(:, sort_idx_ROI_1);
end

% Sort trial_lat_ROI_2 by the values in the 2nd row
if ~isempty(trial_lat_ROI_2)
    [~, sort_idx_ROI_2] = sort(cellfun(@str2double, trial_lat_ROI_2(2, :)));
    trial_lat_ROI_2 = trial_lat_ROI_2(:, sort_idx_ROI_2);
end


latencies = unique(trial_latencies(~cellfun('isempty', trial_latencies)));
n_latencies = size(latencies,1);
stimtype = unique(trial_stim_type);
n_stimtype = size(stimtype,1)-1;

count_ROI_1 = 1;
count_ROI_2 = 1;

for i = 1:length(data.Session.UniqueConds)
    cond = data.Session.UniqueConds{i};
    if contains(cond, 'ROI_1')
        ROI_1_idx{count_ROI_1} = i;
        count_ROI_1 = count_ROI_1 + 1;
    elseif contains(cond, 'ROI_2')
        ROI_2_idx{count_ROI_2} = i;
        count_ROI_2 = count_ROI_2 + 1;
    end
end

% These arrays contain data.Sessions.OdorResponse divided by the
% ROI that was used for stimulation
StimResponse_ROI_1 = cell(1,size(ROI_1_idx,2));
StimResponse_ROI_2 = cell(1,size(ROI_2_idx,2));

for i = 1:length(ROI_1_idx)
    idx = ROI_1_idx{i};
    if ~isempty(idx)
        StimResponse_ROI_1{i} = data.Session.OdorResponse{idx};
    end
end

for i = 1:size(ROI_2_idx,2);
    idx = ROI_2_idx{i};
    StimResponse_ROI_2{i} = data.Session.OdorResponse{idx};
end

StimResponse_MEAN_ROI_1 = cell(1,size(ROI_1_idx,2));
StimResponse_UNC_ROI_1 = cell(1,size(ROI_1_idx,2));
StimResponse_MEAN_ROI_2 = cell(1,size(ROI_2_idx,2));
StimResponse_UNC_ROI_2 = cell(1,size(ROI_1_idx,2));

% Calculate the mean across the 3rd dimension
for i = 1:length(StimResponse_ROI_1)
    if ~isempty(StimResponse_ROI_1{i})
        StimResponse_MEAN_ROI_1{i} = mean(StimResponse_ROI_1{i}, 3);
        StimResponse_UNC_ROI_1{i} = std(StimResponse_ROI_1{i}, 0, 3) / sqrt(size(StimResponse_ROI_1{i}, 3));
    end
end

for i = 1:length(StimResponse_ROI_2)
    if ~isempty(StimResponse_ROI_2{i})
        StimResponse_MEAN_ROI_2{i} = mean(StimResponse_ROI_2{i}, 3);
        StimResponse_UNC_ROI_2{i} = std(StimResponse_ROI_2{i}, 0, 3) / sqrt(size(StimResponse_ROI_2{i}, 3));
    end
end

StimResponse_MEAN_ROI_1 = StimResponse_MEAN_ROI_1(sort_idx_ROI_1);
StimResponse_MEAN_ROI_2 = StimResponse_MEAN_ROI_2(sort_idx_ROI_2);

STIM_RESPONSE_MEAN = cell(n_stimtype, length(ROI_1_idx));
STIM_RESPONSE_UNC = cell(n_stimtype, length(ROI_1_idx));

STIM_RESPONSE_MEAN(1, :) = StimResponse_MEAN_ROI_1;
STIM_RESPONSE_UNC(1, :) = StimResponse_UNC_ROI_1;
STIM_RESPONSE_MEAN(2, :) = StimResponse_MEAN_ROI_2;
STIM_RESPONSE_UNC(2, :) = StimResponse_UNC_ROI_2;

STIM_latencies(1,:) = trial_lat_ROI_1(2,:);
STIM_latencies(2,:) = trial_lat_ROI_2(2,:);

%% NORMALIZE STIM_RESPONSE_MEAN

global_max = -inf;

% Find the global maximum value in the entire cell array
for row = 1:n_stimtype
    for col = 1:n_latencies
        current_max = max(max(STIM_RESPONSE_MEAN{row, col}));
        if current_max > global_max
            global_max = current_max;
        end
    end
end

% Normalize the entire cell array using the global maximum value
for row = 1:n_stimtype
    for col = 1:n_latencies
        STIM_RESPONSE_MEAN_norm{row, col} = STIM_RESPONSE_MEAN{row, col} / global_max;
    end
end

%% Creating no_stim variable

% Initialize the matrix to store the sum of all matrices
accumulated_sum = zeros(size(STIM_RESPONSE_MEAN_norm{2, 1}));

% Loop through each cell and sum the matrices
for col = 1:n_latencies
    accumulated_sum = accumulated_sum + STIM_RESPONSE_MEAN_norm{2, col};
end

% Calculate the average matrix
no_stim = accumulated_sum / n_latencies;


%% TiME CROSS-SECTION PLOT - STIM

% Create a new figure
figure_handle = figure('Name', 'ROI Subplots');

if strcmp(layer_acr, 'MC')
% Set the figure dimensions: [left, bottom, width, height]
    set(figure_handle, 'Position', [10, 100, 400 * length(ROIindex), 400]);
elseif strcmp(layer_acr, 'Glom')
    set(figure_handle, 'Position', [10, 100, 1000, 400]);
end

% Create a tiled layout with one row and num_ROIs columns
tiledlayout(1, length(ROIindex), 'Padding', 'compact', 'TileSpacing', 'compact');

% Loop through each ROI
for i = 1:length(ROIindex);
    ROI_idx = ROIindex(i);
    
    % Create a new subplot for each ROI
    nexttile;
    
    hold on;

    % Loop through each latency and plot the current ROI column of each matrix in OdorResponse_mean
    for latency_idx = 1:n_latencies
        currentMean = STIM_RESPONSE_MEAN_norm{1, latency_idx}(:, ROI_idx);           
        legend_label = sprintf('%s ms', STIM_latencies{n_stimtype, latency_idx});
        
        % Set the alpha value based on the latency index
        alpha_value = 1 - (latency_idx - 1) * 0.15;
        
        % Create a plot with transparency
        h = line(tt, currentMean, 'LineWidth', 1.5, 'DisplayName', legend_label, 'Color', 'r');
        h.Color(4) = alpha_value;  % Set the alpha value

        hold on;

        xlim([-0.5, 1]);
        ylim([-0.25, 1.2]); % Adjust as needed

        % Dashed vertical line at x = 0
        line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '-', 'LineWidth', 1, 'HandleVisibility', 'off');
        line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

        % Add labels, title, etc. as needed
        xlabel('Time Post Inhilation Onset [s]');

        if i == 1
            ylabel('Normalized \textbf{\textsf{$\Delta F/F_0$}}','Interpreter','latex','FontName','Arial','FontSize', axis_size);
        end
    end
    
    % Plot the no_stim data for the current ROI
    no_stim_ROI = no_stim(:, ROI_idx); 

    % Plot the no_stim data as a dashed line
    plot(tt, no_stim_ROI, 'LineWidth', 1.5, 'DisplayName', 'no stim', 'Color', 'k', 'LineStyle', '-');

    % Add legend for the first subplot only
    if strcmp(layer_acr, 'Glom')
        legend('Location', 'southeastoutside', 'Orientation', 'vertical', 'EdgeColor', 'none');
    end

end

% Set the main title for the entire figure
if strcmp(layer_acr, 'MC')
    sgtitle('Daughter Mitral Cells');
elseif strcmp(layer_acr, 'Glom')
    sgtitle('Target Glomerulus');
end

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('%s_dF_plot.png', layer_acr)), '-dpng', '-r300');

%% FLUORESCENCE AMPLITUDE PLOTS FOR EACH ROI [0.5 sec average post-stim presentation]

% ORGANIZE DATA

%Convert Stim Latencies [ms] to [frame count] 
stim_latency_frame = cell(size(STIM_latencies));

for i = 1:size(STIM_latencies, 1)
    for j = 1:size(STIM_latencies, 2)
        X = str2double(STIM_latencies{i, j});
        stim_latency_frame{i, j} = 30 + ((30 * X) / 1000);
    end
end

% Round all frame values to the nearest whole number
stim_latency_frame_aprx = cellfun(@round, stim_latency_frame, 'UniformOutput', false);


%% Find the average fluorescence 15 frames (0.5sec) post stim presentation 
POST_STIM_MEAN = cell(size(STIM_RESPONSE_MEAN));

stim_time = [10, 30, 60, 120, 180, 240];

% stim_time to frame_time
stim_frame = (stim_time * 30 / 1000);
stim_frame = pre_inh + round(stim_frame);
stim_frame_end = stim_frame + 10;
stim_frame_end = round(stim_frame_end);

for i = 1:size(STIM_RESPONSE_MEAN_norm, 1)
    for j = 1:n_latencies
        % Get the current STIM_RESPONSE_MEAN matrix
        response_mean = STIM_RESPONSE_MEAN_norm{i, j};
        
        % Get the corresponding latency frame from stim_latencies_frame_aprx
        latency_frame = stim_latency_frame_aprx{i, j};
        
        % Initialize an array to store the averages
        avg_values = zeros(1, size(response_mean, 2));
        
        % Loop through each column and calculate the average of the 15 rows
        % after the latency frame index
        for col = 1:size(response_mean, 2)
            
            % Calculate the mean of the 15 rows (or fewer if at the matrix edge)
            avg_values(col) = mean(response_mean(stim_frame(j):stim_frame_end(j), col));
        end
        
        % Store the result in POST_STIM_MEAN
        POST_STIM_MEAN{i, j} = avg_values;
    end
end  

%% NORMALIZE POST_STIM_MEAN

global_max = -inf;

% Find the global maximum value in the entire cell array
for row = 1:n_stimtype
    for col = 1:n_latencies
        current_max = max(POST_STIM_MEAN{row, col});
        if current_max > global_max
            global_max = current_max;
        end
    end
end

% Normalize the entire cell array using the global maximum value
for row = 1:n_stimtype
    for col = 1:n_latencies
        POST_STIM_MEAN_norm{row, col} = POST_STIM_MEAN{row, col} / global_max;
    end
end

%% CREATE MASKS MATRIX

Masks = cell(n_stimtype, n_ROI);

for i = 1:n_stimtype
    % Store the 512x512 matrices into the corresponding row of Masks
    for j = 1:n_ROI
         Masks{i, j} = reshape(data.Session.CellMask(:, j), [512, 512]);
    end
   
end
   
% imshow(Masks{2,6},[])

% reshape Masks
Masks_F = cell(size(Masks, 1), size(Masks, 2), n_latencies);

for i = 1:size(Masks, 1)
    for j = 1:size(Masks, 2)
        % Get the original matrix from the Masks array
        original_matrix = Masks{i, j};
        
        % Store the original matrix in each cell of the third dimension
        for k = 1:n_latencies
            Masks_F{i, j, k} = original_matrix;
        end
    end
end

% imshow(Masks_F{2,6,1}{1}, [])


%% Attribute Fluorecence value to TIFF files

% Loop through each column of POST_STIM_MEAN
for lat = 1: n_latencies
    % Get the 2x1 column(i) attributed to POST_STIM_MEAN
    post_stim_mean_col = POST_STIM_MEAN_norm(:, lat);
    
    % Loop through the matrix found at Masks_F ( : , : , i)
    for stmtyp = 1:size(Masks_F, 1)
        for roi = 1:size(Masks_F, 2)
            % Get the current mask matrix
            current_mask = Masks_F{stmtyp, roi, lat};
            
            % Get the value to assign from POST_STIM_MEAN
            value_to_assign = post_stim_mean_col{stmtyp}(roi);
            
            % Loop through the mask matrix and update non-zero values
            for s = 1:size(current_mask, 1)
                for r = 1:size(current_mask, 2)
                    if current_mask(s, r) ~= 0
                        current_mask(s, r) = value_to_assign;
                    end
                end
            end
            
            % Store the updated mask back into Masks_F
            Masks_F{stmtyp, roi, lat} = current_mask;
        end
    end
end

% imshow(Masks_F{2,20,4}, [])

%  unique(Masks_F{2,20,2}) = 0, 0.0089
% POST_STIM_MEAN{2,2}(20) = 0.0089
% unique(Masks_F{1,12,3}) = 0, 0.0384
% POST_STIM_MEAN{1,3}(12) = 0.0384
% This code works properly

%% Add all Masks Together

sum_Masks_F = cell(size(Masks_F, 1), size(Masks_F, 3));

% Loop through each element in the third dimension
for i = 1:size(Masks_F, 1)
    for k = 1:size(Masks_F, 3)
        % Initialize the sum matrix with zeros
        sum_matrix = zeros(size(Masks_F{i, 1, k}));
        
        % Sum across the second dimension
        for j = 1:size(Masks_F, 2)
            sum_matrix = sum_matrix + Masks_F{i, j, k};
        end
        
        % Store the result in the sum_Masks_F array
        sum_Masks_F{i, k} = sum_matrix;
    end
end

%imshow(Masks_F{2,20,4}, [])

% unique(sum_Masks_F{1,3}) = -0.0067
% size(unique(sum_Masks_F{1,3})) = 26 1
% matches POST_STIM_MEAN{1,3}

%% Create no-stim subplot

accumulated_sum = zeros(512, 512);

% Loop through the 2nd row of the cell array and accumulate the sum
for col = 1:6
    accumulated_sum = accumulated_sum + sum_Masks_F{2, col};
end

% Calculate the average
Mask_F_no_stim = accumulated_sum / 6;


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
set(fig, 'Position', [20, 50, 1600, 300]);

% Create a tiled layout with two columns (one for each subplot)
plt_layout = tiledlayout(1, size(sum_Masks_F,2) + 1, 'TileSpacing', 'none', 'Padding', 'none'); % +1 for the no-stim subplot
plt_layout.TileSpacing = 'compact';
plt_layout.Padding = 'compact';

% Convert the values to indices in the colorScale
color_indices_no_stim = arrayfun(@(x) findClosestIndex(x, valueRange), Mask_F_no_stim);

% Create an RGB image based on the colorScale
color_image_no_stim = ind2rgb(color_indices_no_stim, colorScale);

% Plot Mask_F_no_stim in the first tile
nexttile;
imshow(color_image_no_stim);
ylabel(sprintf('%s', layer_name),'FontSize', axis_size);
if strcmp(layer_acr, 'MC')
    xlabel(sprintf('no-stim'),'FontSize', axis_size);
end

for k = 1:size(sum_Masks_F, 2)
    % Map the values directly to the range of valueRange
    data_matrix = sum_Masks_F{1, k};

    % Convert the values to indices in the colorScale
    color_indices = arrayfun(@(x) findClosestIndex(x, valueRange), data_matrix);

    % Create an RGB image based on the colorScale
    color_image = ind2rgb(color_indices, colorScale);

    % Create subplots for each path
    nexttile;
    imshow(color_image);

    hold on;
    for j = 1:length(ROIindex)
        logical_mask = Masks_F{1, ROIindex(j), k} ~= 0;
        contour(logical_mask, [1, 1], 'LineColor', 'k', 'LineWidth', 1.00);
    end

    % Ensure latency is a numeric value
    latency = STIM_latencies{1, k};
    if iscell(latency)
        latency = cell2mat(latency);
    end

    if strcmp(layer_acr, 'MC')
    % Add x-axis label to each subplot
        xlabel(sprintf('%s ms Latency', latency),'FontSize', axis_size);
    end
    
end


% Adjust colormap to use the generated color scale
colormap(gca, colorScale);
cb = colorbar('eastoutside');
caxis([min(valueRange) max(valueRange)]);
cb.Ticks = [-1, 0, 1];

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('%s_mask_plot.png', layer_acr)), '-dpng', '-r300');

%% STIM-TIME AVERAGE 
% plot POST_STIM_MEAN_norm
% reconfigure no-stim variable 
stim_time = [10, 30, 60, 120, 180, 240];

% STIM_ONLY dimention = (n_ROI , n_latencies)
STIM_ONLY = cell(length(ROIindex), n_latencies);

for i = 1:length(ROIindex)
    for j = 1:n_latencies
        r = ROIindex(i);  % Get the specific ROI index
        STIM_ONLY{i,j} = POST_STIM_MEAN{1,j}(r);  % Assign the value
    end
end

% Convert cell array STIM_ONLY_no_stim to a numeric array
STIM_ONLY_numeric = cellfun(@(x) x, STIM_ONLY);

% Calculate the mean across the second dimension
STIM_ONLY_mean = mean(STIM_ONLY_numeric, 1);


STIM_ONLY_no_stim = cell(length(ROIindex), n_latencies);

for i = 1:length(ROIindex)
    for j = 1:n_latencies
        r = ROIindex(i);  % Get the specific ROI index
        STIM_ONLY_no_stim{i,j} = POST_STIM_MEAN{2,j}(r);  % Assign the value
    end
end


% Convert cell array STIM_ONLY_no_stim to a numeric array
STIM_ONLY_no_stim_numeric = cellfun(@(x) x, STIM_ONLY_no_stim);

% Calculate the mean across the second dimension
STIM_ONLY_no_stim_mean = mean(STIM_ONLY_no_stim_numeric, 1);

%% PLOT STIM-TIME AVERAGE PLOT

% no stim or odor (black)
% stim (red)
% low BENZ (blue) (0.5 alpha)
% high BENZ (blue) (1.0 alpha)
% stim + low BENZ (purple) (0.5 alpha)
% stim + high BENZ (purple) (1.0 alpha)

figure_handle = figure('Name', 'ODOR+STIM Plot');

% Set the figure dimensions: [left, bottom, width, height]
set(figure_handle, 'Position', [10, 100, 600, 400]);

% Hold on to plot multiple lines
hold on;

% LEGEND CUSTOM LABEL
if strcmp(layer_acr, 'MC')
    custom_labels = {
        'no stim or odor', [0 0 0 1];    % Black for no stim or odor
        'stim, daughter MC', [1 0 0 1];  % Red for stim
    };
elseif strcmp(layer_acr, 'Glom')
    custom_labels = {
        'no stim or odor', [0 0 0 1];    % Black for no stim or odor
        'stim, single Glom', [1 0 0 1];  % Red for stim
    };
end

% Plot the 'no stim or odor' data
plot(stim_time, STIM_ONLY_no_stim_mean, 'Color', custom_labels{1, 2}, ...
     'LineWidth', 1, 'LineStyle', '-');

% Plot each row of STIM_ONLY with the 'stim' label
for i = 1:size(STIM_ONLY, 1)
    % Convert cell content to a numeric array before plotting
    plot(stim_time, cell2mat(STIM_ONLY(i, :)), 'Color', custom_labels{2, 2}, ...
         'LineWidth', 1, 'LineStyle', '-');
end

% Dashed vertical line at x = 0
line([0 0], ylim, 'Color', [0.8, 0.8, 0.8], 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Horizontal line at y = 0
line(xlim, [0, 0], 'Color', [0, 0, 0, 0.1], 'LineStyle', '-', 'HandleVisibility', 'off');

% Adjust y-axis limits as needed
ylim([-0.2, 1.0]);

% Set x-axis ticks
xticks(stim_time);

% Add labels and title
if strcmp(layer_acr, 'MC')
    xlabel('Stim Onset Latency [ms]', 'FontSize', axis_size);
end
ylabel('Normalized $$\Delta F/F_o$$ (0.5s Averaged)', 'Interpreter', 'latex');

% Create legend
legend(custom_labels{:, 1}, 'Location', 'best', 'EdgeColor', 'none');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('%s_ODOR-STIM_plot.png', layer_acr)), '-dpng', '-r300');

% Finish plotting
hold off;

%% COMBINED STIM-TIME AVERAGE PLOT

% Define the filenames for the two PNG files
file1 = sprintf('Glom_ODOR-STIM_plot.png');
file2 = sprintf('MC_ODOR-STIM_plot.png');

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
blank_space = uint8(0 * ones(gap, max(width1, width2), 3)); % White space

% Combine the images with the blank space
combined_img = [img1; blank_space; img2];

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
title('Single Glomerular Stimulation', 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), size(combined_img, 1) + 30, 'normalized with respect to OB layer', ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('ODOR-STIM_plot.png')), '-dpng', '-r300');

%% COMBINED PLOTS - Fluorescence Cross-Traces

% Define the filenames for the two PNG files
file1 = sprintf('%s_dF_plot.png', 'Glom');
file2 = sprintf('%s_dF_plot.png', 'MC');

% Create full paths to the PNG files
path1 = fullfile(plot_dir, file1);
path2 = fullfile(plot_dir, file2);

% Read the PNG files
img1 = imread(path1);
img2 = imread(path2);

% Determine the dimensions of the two images
[height1, width1, ~] = size(img1);
[height2, width2, ~] = size(img2);

% Determine the width of the combined image
max_width = max(width1, width2);

% Pad the images to the same width
if width1 < max_width
    % Calculate the padding for left and right to center img1
    pad_left = floor((max_width - width1) / 2);
    pad_right = max_width - width1 - pad_left;
    img1 = padarray(img1, [0 pad_left], 255, 'pre');  % Pad img1 on the left
    img1 = padarray(img1, [0 pad_right], 255, 'post'); % Pad img1 on the right
end
if width2 < max_width
    img2 = padarray(img2, [0 max_width - width2], 255, 'post'); % Pad img2 to the right
end

% Define the vertical gap (space) between the images
gap = 50; % Adjust this value as needed for the desired spacing

% Create a blank space of the same width as the images and height of the gap
blank_space = uint8(255 * ones(gap, max_width, 3)); % White space

% Combine the images with the blank space
combined_img = [img1; blank_space; img2];

% Create a new figure
figure;

% Display the combined image
imshow(combined_img);

% Set the title if needed
title('Single Glomerular Stimulation Fluorescence Cross-Traces', 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), size(combined_img, 1) + 30, 'normalized with respect to OB layer', ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('dF_plot.png')), '-dpng', '-r300');


%% COMBINED PLOTS - ROI Masks

% Define the filenames for the two PNG files
file1 = sprintf('Glom_mask_plot.png');
file2 = sprintf('MC_mask_plot.png');

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

title(sprintf('Passive Single Glom Stim Presentation | Normalized 0.5 sec Fluorescence Avg. Post-Stim Delivery'), 'FontSize', title_size);

% Add label at the bottom right for the number of trials
text(size(combined_img, 2), size(combined_img, 1) + 50, sprintf('number of trials: %d', n_trials), ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Add an additional line of text just above the "number of trials" label
text(size(combined_img, 2), size(combined_img, 1) + 30, 'normalized with respect to OB layer', ...
    'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Save the figure
print(gcf, fullfile(plot_dir, sprintf('combined_mask_plot.png')), '-dpng', '-r300');


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