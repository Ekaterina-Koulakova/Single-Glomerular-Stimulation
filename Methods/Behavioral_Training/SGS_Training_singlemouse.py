for var in list(locals().keys()):
    del locals()[var]


import h5py
import re
import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

# =============================================================================
# SELECT MOUSE TO ANALYZE by uncommenting its respective string pattern
# =============================================================================

string_pattern = "*mouse0691*.h5"
# string_pattern = "*mouse0070*.h5"
# string_pattern = "*mouse0071*.h5"
# string_pattern = "*mouse0773*.h5"
# string_pattern = "*mouse0781*.h5"


directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim/Training'
home_directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim'

file_list = os.listdir(directory)

# Extract 'mouse0064' from the 'string_pattern' using regular expressions
match = re.search(r'mouse(\d+)', string_pattern)
mouse_number = match.group(1)


file_paths = glob.glob(os.path.join(directory, string_pattern))
file_paths = [path.replace("\\", "/") for path in file_paths]

print(file_paths)

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
print(sorted_files_sess)

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
    slicing_ranges = [(20,250),  #sess01 230907
                      (20,250),  #sess02 230912
                      (100,250), #sess03 230913
                      (40,110,130,260),  #sess04 230914
                      (0,265),   #sess05 230915
                      (60,230)]  #sess06 230919
    
if mouse_number == '0070':    
    slicing_ranges = [(50,300),  #sess01 231016
                      (50,300),  #sess02 231017
                      (100,250), #sess03 231018
                      (50,225),  #sess04 231019
                      (180,400), #sess05 231020
                      (180,390), #sess06 231023
                      (50,360)]  #sess07 231024
    
if mouse_number == '0071':    
    slicing_ranges = [(50,300),  #sess01 231016
                      (50,300),  #sess02 231017
                      (50,300),  #sess03 231018
                      (150,300), #sess04 231019
                      (50,400),  #sess05 231020
                      (30,320)]  #sess06 231023
    
    
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
# NAMING VARIABLES FROM DATA 
# =============================================================================

# trialtype is a list of arrays, where each array represents the "trialtype" column from each dataframe in the "data" list
trialtype = [df['trialtype'].values for df in data]
# iterates over each dataframe in the "data" list. For each dataframe, it 
# accesses the 'trialtype' column using df['trialtype'] and retrieves its values 
# using the .values attribute.

result = [df['result'].values for df in data]
mask_id = [df['mask_id'].values for df in data]
pavlov = [df['pavlov'].values for df in data]
# number of sessions conducted
n_sessions = len(data)

n_stim = np.array([len(arr) - np.sum(arr) for arr in trialtype])
n_blank = np.array([np.sum(arr) for arr in trialtype])
n_trialtype = pd.DataFrame({'n_stim': n_stim, 'n_blank': n_blank})

tot_trials_sess = np.array([len(trials) for trials in trialtype], dtype=float)


# =============================================================================
# CREATE AN ARRAY OF PERCENT PAVLOV FOR EACH SESSION
# =============================================================================
# sum of pavlov/sum of match trials
# nonmatch trials were included as a check (variables created should store 0's)

pavlov_percent = []

for i in range(n_sessions):
    p_pavlov = (sum(pavlov[i])/len(pavlov[i]))*100
    pavlov_percent.append(p_pavlov) 

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
                

# MATCH TRIALS SUCCESS RATE

ROI_values_match = [value for value in ROI_values if value != 99]
keys_no_99 = [key for key in SUC_ROI_dict.keys() if key != 99]


for i in range(n_sessions):
    for x in range(len(trialtype[i])):
        for v in ROI_values_match:
            
            if result[i][x] != 2:
                SUC_ROI_dict[v][i].append(0)

            if mask_id[i][x] == v and result[i][x] == 2:
                SUC_ROI_dict[v][i].append(1)
                        
            if mask_id[i][x] != v and result[i][x] == 2:
                SUC_ROI_dict[v][i].append(0)
                
# =============================================================================
# # CHECK
# result_list = SUC_ROI_dict[0][3]
# length_of_list = len(result_list)
# print(length_of_list)
# length_of_trialtype_0 = len(trialtype[3])
# print(length_of_trialtype_0)
# =============================================================================
                                
# CODE TO SUM SUC_ROI_dict

sum_SUC_ROI_dict = {}

for key, binary_lists in SUC_ROI_dict.items():
    sum_lists = [sum(binary_list) for binary_list in binary_lists]
    sum_SUC_ROI_dict[key] = sum_lists


# Convert the dictionary into a DataFrame
sum_SUC_ROI_df = pd.DataFrame(sum_SUC_ROI_dict, index=range(1, len(mask_id)+1))
            
# Fraction of successful trials for each ROI
frac_SUCCESS_ROI = sum_SUC_ROI_df/ROI_sum_df

# Fraction of successful trials for blanks
frac_SUCCESS_Blanks = frac_SUCCESS_ROI[99].to_numpy()

# Fraction of successful trials for stim sets
frac_SUCCESS_ROI_stim = frac_SUCCESS_ROI.drop(columns=99)

# =============================================================================
# FIND SUCCESS RATE OF ALL ROIs COMBINED
# =============================================================================
# .drop(columns=[99]) to exclude Blank trials

sum_match_stim = ROI_sum_df.drop(columns=[99]).sum(axis=1).to_numpy()
# axis 0 represents rows and axis 1 represents columns. Here, we are summing over the columns.

sum_match_stim_SUC = sum_SUC_ROI_df.drop(columns=[99]).sum(axis=1).to_numpy()

SUCCESS_match_stim = sum_match_stim_SUC/sum_match_stim

avg_SUCCESS_plot = (SUCCESS_match_stim + frac_SUCCESS_Blanks) / 2.0


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

# Loop over the avg_SUCCESS_plot (p) and tot_trials_sess (N) arrays
for p, N in zip(avg_SUCCESS_plot, tot_trials_sess):
    # Calculate the standard error using the given formula
    if N > 0:  # Ensure N is positive to avoid division by zero
        SE = np.sqrt(p * (1 - p) / N)
    else:
        SE = 0  # Handle cases where N is zero
    
    # Multiply by the scaling factor k = 1.96 to get the confidence interval
    CI_value = SE * k
    
    # Append the CI value to the CI list
    CI.append(CI_value)

# Convert CI to a NumPy array for further processing
CI = np.array(CI)


    
# =============================================================================
#     Create the confidence_intervals list of arrays, where each element is an 
#     array of two values (upper and lower bounds) for each corresponding value 
#     in avg_SUCCESS_plot and CI
# =============================================================================

# Create the confidence_intervals list of arrays, where each element is an 
# array of two values (upper and lower bounds) for each corresponding value 
# in avg_SUCCESS_plot and CI
confidence_intervals = []

# Loop through each p-value in avg_SUCCESS_plot and its corresponding CI value
for p, ci in zip(avg_SUCCESS_plot, CI):
    # Calculate the upper and lower bounds of the confidence interval
    upper_bound = p + ci
    lower_bound = p - ci
    
    # Store the upper and lower bounds as an array [upper_bound, lower_bound]
    confidence_intervals.append([upper_bound, lower_bound])

# Convert confidence_intervals to a NumPy array
confidence_intervals = np.array(confidence_intervals) 


# =============================================================================
# PLOT SHOWING THE FRACTION OF SUCCESSFUL TRIALS 
# =============================================================================
# Create a single figure with two vertically stacked subplots

x = np.arange(n_sessions) + 1

fig, ax1 = plt.subplots(figsize=(9, 5))

# Plot Average Success Rate on the left y-axis
line1, = ax1.plot(x, avg_SUCCESS_plot, 'o-', label='Avg. Success Rate', color='purple')
ax1.axhline(y=0.0, color='gray', linestyle='--', linewidth=0.8)
ax1.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)
ax1.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8)

# Extract the confidence intervals for the average success rate
lower_bounds = confidence_intervals[:, 1]  # Lower bounds
upper_bounds = confidence_intervals[:, 0]  # Upper bounds

# Plot the shaded error region for the confidence interval
ax1.fill_between(x, lower_bounds, upper_bounds, color='red', alpha=0.3)

# Left y-axis (Fraction of Successful Match Trials)
ax1.set_ylabel('Fraction of Successful Trials', fontsize=14)
ax1.set_ylim(-0.1, 1.1)  # Setting limits from slightly below 0 to slightly above 1
ax1.set_xlim(0.5, n_sessions + 0.5)
ax1.set_xticks(x)

# Create a second y-axis on the right for Percent Pavlov
ax2 = ax1.twinx()
ax2.plot(x, pavlov_percent, 'g*', label='Percent Pavlov', markersize=10)

# Right y-axis (Percentage of Pavlov Trials)
ax2.set_ylabel('Percentage of Pavlov Trials', fontsize=14)
ax2.set_ylim(-10, 110)  # Pavlov percentage range from 0 to 100
ax2.set_yticks(np.arange(0, 101, 20))  # Ticks at 0, 20, 40, 60, 80, 100

# Title and legend
ax1.set_title(f'Training Single Glomerular Stimulation | Mouse {mouse_number}', fontsize=16)

# Add legends for both lines
fig.legend(loc='lower left', bbox_to_anchor=(0.12, 0.12), fontsize='large')

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
save_path = os.path.join(save_directory, f'SGS_Training_mouse{mouse_number}.png')

# Save the plot
fig.savefig(save_path, dpi=300)