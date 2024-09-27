import os
import h5py
import re
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt


string_patterns = ["*mouse0691*.h5", "*mouse0070*.h5", "*mouse0071*.h5", "*mouse0721*.h5", "*mouse0722*.h5", "*mouse0773*.h5", "*mouse0781*.h5"]

avg_SUCCESS_plot = []
mouse = []
n_stim_sess = []
n_blank_sess = []
frac_SUCCESS_stim = []
frac_SUCCESS_blank = []

for string_pattern in string_patterns:

    directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim/Training'
    home_directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim'
    
    file_list = os.listdir(directory)
    
    # Extract 'mouse0064' from the 'string_pattern' using regular expressions
    match = re.search(r'mouse(\d+)', string_pattern)
    mouse_number = match.group(1)
    mouse.append(mouse_number)
    
    file_paths = glob.glob(os.path.join(directory, string_pattern))
    file_paths = [path.replace("\\", "/") for path in file_paths]
    
    
    # =============================================================================
    # SORT FILES BY SESSION NUMBER (indicated in title of file)
    # =============================================================================
    
    def sort_key(file_path):
       # Extract the 4 digits preceding "V_"
       start_index = file_path.find("sess") + 4
       end_index = start_index + 2
       number_str = file_path[start_index:end_index]
       
       # Convert the number string to an integer
       try:
           number = int(number_str)
       except ValueError:
           number = 0
       
       return number
    
    sorted_files_sess = sorted(file_paths, key=sort_key)
    
    # Create a list of dataframes using the various file_paths.
    data = []
    
    for file_path in sorted_files_sess:
        with h5py.File(file_path, 'r') as h5_file:
            dataset = h5_file['Trials']
            data_frame = pd.DataFrame(dataset[:]) 
            data.append(data_frame)
              
    # =============================================================================
    # # iterating the columns
    # for col in data[1].columns:
    #     print(col) 
    # =============================================================================
    
    #SELECT THE RANGE OF TRIALS THAT YOU WOULD LIKE TO CONSIDER FOR EACH DATAFRAME IN 'data'
    # Define the desired slicing ranges for each dataframe
    slicing_ranges = []
        
    if mouse_number == '0691':    
        slicing_ranges = [(20,250),  #sess01 23.09.07
                          (20,250),  #sess02 23.09.12
                          (100,250), #sess03 23.09.13
                          (40,110,130,260),  #sess04 23.09.14
                          (0,265),   #sess05 23.09.15
                          (60,230)]  #sess06 23.09.19
        
    if mouse_number == '0070':    
        slicing_ranges = [(50,300),  #sess01 23.10.16
                          (50,300),  #sess02 23.10.17
                          (100,250), #sess03 23.10.18
                          (50,225),  #sess04 23.10.19
                          (180,400), #sess05 23.10.20
                          (180,390), #sess06 23.10.23
                          (50,360),  #sess07 23.10.24
                          (20,150)]  #sess08 23.11.02
        
    if mouse_number == '0071':    
        slicing_ranges = [(50,300),  #sess01 23.10.16
                          (50,300),  #sess02 23.10.17
                          (50,300),  #sess03 23.10.18
                          (150,300), #sess04 23.10.19
                          (50,400),  #sess05 23.10.20
                          (30,320)]  #sess06 23.10.23
        
    if mouse_number == '0721':    
        slicing_ranges = [(50,400),  #sess01 24.02.12
                          (50,300),  #sess02 24.02.13
                          (50,350),  #sess03 24.02.14
                          (50,370),  #sess04 24.02.15
                          (40,110),  #sess05 24.02.21
                          (50,250)]  #sess06 24.02.22                     
    
    if mouse_number == '0722':    
        slicing_ranges = [(50,175),  #sess01 24.02.12
                          (50,220),  #sess02 24.02.13
                          (40,200),  #sess03 24.02.14
                          (50,370),  #sess04 24.02.15
                          (100,300), #sess05 24.02.26
                          (150,400), #sess06 24.02.27
                          (25,180)]  #sess07 24.02.28
    
    if mouse_number == '0773':    
        slicing_ranges = [(50,300),  #sess01 24.08.13
                          (50,500),  #sess02 24.08.14
                          (50,600),  #sess03 24.08.15
                          (50,500),  #sess04 24.08.16
                          (50,400),  #sess05 24.08.19
                          (50,145)]  #sess06 24.08.20

        
    if mouse_number == '0781':    
        slicing_ranges = [(20,120),  #sess01 24.08.14
                          (20,120),  #sess02 24.08.15
                          (50,450),  #sess03 24.08.16
                          (50,600),  #sess04 24.08.19
                          (50,550),  #sess05 24.08.20
                          (50,225),  #sess06 24.08.21
                          (25,225),  #sess07 24.08.22
# =============================================================================
#                           (25,200),  #sess08 24.08.23
#                           (50,400),  #sess09 24.08.26
#                           (50,500),  #sess10 24.08.27
# =============================================================================
                          (50,400)]  #sess11 24.08.28
                              
    
    # Apply the slicing to each dataframe in the 'data' list
    for i, ranges in enumerate(slicing_ranges):
        if len(ranges) == 2:
            start_idx, end_idx = ranges
            data[i] = data[i].loc[start_idx:end_idx]
        elif len(ranges) == 4:
            start_idx1, end_idx1, start_idx2, end_idx2 = ranges
            data[i] = pd.concat([data[i].loc[start_idx1:end_idx1], data[i].loc[start_idx2:end_idx2]]) 
    
    # =============================================================================
    # REMOVE TRIALS WHERE PAVLOV = 1
    # =============================================================================

    for df in data:
        df.drop(df[df['pavlov'] == 1].index, inplace=True) # inplace=True : ensures that changes are made directly in the original dataframes. 
    
    # =============================================================================
    # NAMING VARIABLES FROM DATA 
    # =============================================================================
    
    # trialtype is a list of arrays, where each array represents the "trialtype" column from each dataframe in the "data" list
    trialtype = [df['trialtype'].values for df in data]
    # iterates over each dataframe in the "data" list. For each dataframe, it 
    # accesses the 'trialtype' column using df['trialtype'] and retrieves its values 
    # using the .values attribute.
    
    result = [df['result'].values for df in data]
    mask_id = [df['mask_id'].values for df in data]
    
    # Create empty arrays to store the sums
    n_stim = []
    n_blank = []
    
    for a in mask_id:
    # Count occurrences of 99 (stimulus) and 0 (blank)
        n_stim_count = np.sum(a == 0)
        n_blank_count = np.sum(a == 99)

    # Append counts to the respective arrays
        n_stim.append(n_stim_count)
        n_blank.append(n_blank_count)

    # Convert the lists to NumPy arrays if needed
    n_stim = np.array(n_stim)
    n_blank = np.array(n_blank)
    length_sess = n_stim + n_blank
    # Append n_stim and n_blank arrays for the current string_pattern
    n_stim_sess.append(n_stim)
    n_blank_sess.append(n_blank)
   
    # Ensure that tot_trials_sess will store the element-wise sum of n_stim_sess and n_blank_sess
    tot_trials_sess = []
    
    # Loop through the lists of arrays (n_stim_sess and n_blank_sess)
    for stim_array, blank_array in zip(n_stim_sess, n_blank_sess):
        # Ensure both arrays have the same shape before adding
        if stim_array.shape == blank_array.shape:
            combined_array = stim_array + blank_array  # Element-wise sum
            tot_trials_sess.append(combined_array)     # Append the result
    
    # number of sessions conducted
    n_sessions = len(data)
    
    # =============================================================================
    # PRINT ALL INDECES PRESENT IN mask_id
    # =============================================================================
    
    # 1: Convert numpy arrays to pandas DataFrames if needed
    mask_id_df = [pd.DataFrame(df) if isinstance(df, np.ndarray) else df for df in mask_id]
    
    # 2: Concatenate all dataframes into a single dataframe or series
    concatenated_data = pd.concat(mask_id_df, ignore_index=True)
    
    # 3: Remove duplicate values and sort in numerical order
    ROI_values = sorted(concatenated_data[0].unique())
    
    n_ROIs = len(ROI_values)
    
    print(ROI_values)
    
    # =============================================================================
    # FIND OUT THE SUCCESS RATES of EACH ROI
    # =============================================================================

    # =============================================================================
    # SUM of how many times each mask_id is present in each session.
    # =============================================================================
        
    number_counts = {number: [] for number in ROI_values}
    
    # Iterate through each DataFrame in 'mask_id'
    for df in mask_id_df:
        if 'Numbers' in df.columns:
            column_name = 'Numbers'
        else:
            column_name = df.columns[0]  # Get the column name dynamically
    
        for number in ROI_values:
            count = df[column_name].eq(number).sum()
            number_counts[number].append(count)
    
    # This df contains columns with ROI indeces and rows that store the number of times each ROI was presented for each voltage
    ROI_sum_df = pd.DataFrame(number_counts, index=range(1, len(mask_id)+1))
    
    # =============================================================================
    # TABULATE SUCCESS RATES OF VARIOUS ROIS
    # =============================================================================
    
    # Create a list of SUC_mask_id_## variables by looping through the ROIs listed in 'mask_id_values'
    for value in ROI_values:
        var_name = f"SUC_ROI_{value:02}"
        locals()[var_name] = [[] for _ in range(n_sessions)]
    
    # Create a dictionary to store the lists of SUC_ROI values
    SUC_ROI_dict = {value: [[] for _ in range(n_sessions)] for value in ROI_values}
    
    
    # BLANK TRIALS SUCCESS RATE
        
    for i in range(n_sessions):
        for x in range(len(trialtype[i])):         
                
                if mask_id[i][x] == 99 and result[i][x] == 1:
                    SUC_ROI_dict[99][i].append(1)
                    
                if mask_id[i][x] == 99 and result[i][x] != 1:
                    SUC_ROI_dict[99][i].append(0)
                
                if mask_id[i][x] != 99:
                    SUC_ROI_dict[99][i].append(0)
                    
    
    # STIM TRIALS SUCCESS RATE
    
    ROI_values_stim = [value for value in ROI_values if value != 99]
    keys_no_99 = [key for key in SUC_ROI_dict.keys() if key != 99]
    
    
    for i in range(n_sessions):
        for x in range(len(trialtype[i])):
            for v in ROI_values_stim:
                
                if result[i][x] != 2:
                    SUC_ROI_dict[v][i].append(0)
    
                if mask_id[i][x] == v and result[i][x] == 2:
                    SUC_ROI_dict[v][i].append(1)
                            
                if mask_id[i][x] != v and result[i][x] == 2:
                    SUC_ROI_dict[v][i].append(0)
                    
                                    
    # CODE TO SUM SUC_ROI_dict
    
    sum_SUC_ROI_dict = {}
    
    for key, binary_lists in SUC_ROI_dict.items():
        sum_lists = [sum(binary_list) for binary_list in binary_lists]
        sum_SUC_ROI_dict[key] = sum_lists
    
    
    # Convert the dictionary into a DataFrame
    sum_SUC_ROI_df = pd.DataFrame(sum_SUC_ROI_dict, index=range(1, len(mask_id)+1))
                
    # Fraction of successful trials for each ROI
    frac_SUCCESS_ROI = sum_SUC_ROI_df/ROI_sum_df
    
    # Array of fraction of successful trials for Blanks (mask_id = 99)
    frac_SUCCESS_Blanks = frac_SUCCESS_ROI[99].to_numpy()
    
    # Array of fraction of successful trials for Blanks (mask_id = 0)
    frac_SUCCESS_Stim = frac_SUCCESS_ROI[0].to_numpy()
    # This line of code replaces the NaN values in frac_SUCCESS_Stim with zeros
    # NaN values are there bc for two sessions for one of the mice, there were no Stim trials without the use of pavlov.
    frac_SUCCESS_Stim = np.nan_to_num(frac_SUCCESS_Stim, nan=0)
    
    # =============================================================================
    # FIND SUCCESS RATE OF ALL ROIs COMBINED
    # =============================================================================
    
    avg_SUCCESS_frac = (frac_SUCCESS_Stim + frac_SUCCESS_Blanks) / 2.0
    
    # This line of code is meant to append the success rate arrays of various other mice.
    avg_SUCCESS_plot.append(avg_SUCCESS_frac)
    
    
# =============================================================================
# CONFIDENCE INTERVALS
# =============================================================================
# 1. Calculate 95% confidence interval
# https://www.dummies.com/article/academics-the-arts/science/biology/the-confidence-interval-around-a-proportion-149351/


# =============================================================================
# To calculate the Standard Error (SE) for each value in avg_SUCCESS_plot (which 
# corresponds to p in the equation) using the values in tot_trials_sess (which 
# corresponds to N), we can loop through each element of both avg_SUCCESS_plot 
# and tot_trials_sess. The formula for SE is:
# SE = sqrt(p(1-p)/N)
# This equation should be applied to each element in avg_SUCCESS_plot (p) and 
# tot_trials_sess (N) and the output should be stored in SE_plot
# Multiply SE by a preselected k-value to form CI (which will now represent 
# the confidence interval.
# =============================================================================


# Define the constant k for calculating the confidence interval
k = 1.96
# k is 1.96 for normal-based 95 percent confidence limits.

# Initialize CI to store the calculated confidence intervals
CI = []

# Loop over the avg_SUCCESS_plot (p) and tot_trials_sess (N) lists
for avg_success_array, tot_trials_array in zip(avg_SUCCESS_plot, tot_trials_sess):
    # Create an array to store the confidence intervals for the current session
    ci_array = []
    
    # Loop over each value of p and N in the current arrays
    for p, N in zip(avg_success_array, tot_trials_array):
        # Calculate the standard error using the given formula
        if N > 0:  # Ensure N is positive to avoid division by zero
            SE = np.sqrt(p * (1 - p) / N)
        else:
            SE = 0  # Handle cases where N is zero
        
        # Multiply by the scaling factor k = 1.96 to get the confidence interval
        CI_value = SE * k
        
        # Append the CI value to the array
        ci_array.append(CI_value)
    
    # Append the CI array for this session to CI
    CI.append(ci_array)
        

    
# =============================================================================
#     Create the confidence_intervals list of arrays, where each element is an 
#     array of two values (upper and lower bounds) for each corresponding value 
#     in avg_SUCCESS_plot and CI
# =============================================================================
    
    confidence_intervals = []

    # Loop through each session in avg_SUCCESS_plot and CI
    for avg_success_array, ci_array in zip(avg_SUCCESS_plot, CI):
        # Create an array to store confidence intervals for the current session
        session_confidence_intervals = []
        
        # Loop through each p-value in avg_SUCCESS_array and its corresponding CI value
        for p, ci in zip(avg_success_array, ci_array):
            # Calculate the upper and lower bounds of the confidence interval
            upper_bound = p + ci
            lower_bound = p - ci
            
            # Store the upper and lower bounds as an array [upper_bound, lower_bound]
            session_confidence_intervals.append([upper_bound, lower_bound])
        
        # Append the confidence intervals for this session to the main list
        confidence_intervals.append(session_confidence_intervals)


# =============================================================================
# PLOT SHOWING THE FRACTION OF SUCCESSFUL TRIALS 
# =============================================================================
# Create a single figure with one subplot
fig, ax = plt.subplots(1, 1, figsize=(9, 5))

# Get the maximum length of the sublists in avg_SUCCESS_plot to set x-axis limits
max_x_length = max(len(sublist) for sublist in avg_SUCCESS_plot)

# Loop through avg_SUCCESS_plot and plot each list for each mouse
for i, sublist in enumerate(avg_SUCCESS_plot):
    x = np.arange(len(sublist)) + 1  # x-axis values (session numbers)
    
    # Plot the data points and capture the color of the line
    line, = ax.plot(x, sublist, 'o-', label=f'Mouse {mouse[i]}', alpha=0.7)
    line_color = line.get_color()  # Get the color of the current line

    # Extract the confidence intervals for the current session
    confidence_interval_array = confidence_intervals[i]
    
    # Separate the upper and lower bounds from confidence_interval_array
    lower_bounds = [ci[1] for ci in confidence_interval_array]
    upper_bounds = [ci[0] for ci in confidence_interval_array]

    # Plot the shaded error region for the confidence interval
    ax.fill_between(x, lower_bounds, upper_bounds, color=line_color, alpha=0.3)

# Add dashed lines at y = 0.0, 0.5, and 1.0 for reference
ax.axhline(y=0.0, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8)

# Set labels, title, and axis limits
ax.set_title(f'Training Single Glomerular Stimulation Detection', fontsize=16)
ax.set_ylabel('Fraction of Successful Trials', fontsize=14)
ax.set_xlabel('Session Number', fontsize=14)
ax.set_ylim(0.4, 1.1)  # Adjust y-axis limits for clarity
ax.set_xlim(0.5, max_x_length + 0.5)  # Set x-axis limits
ax.set_xticks(np.arange(1, max_x_length + 1))  # Ensure x-ticks match session numbers
ax.legend(loc='lower right', fontsize='large')  # Add a legend for the mice

# Display the plot
plt.show()

# =============================================================================
# SAVE PRODUCED PLOTS
# =============================================================================

# Save the figure as a .png file in the specified directory
save_directory = os.path.join(home_directory, 'plots_behavior')
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# Define the full path for the plot
save_path = os.path.join(save_directory, f'SGS_Training_combinedmice.png')

# Save the plot
fig.savefig(save_path, dpi=300)