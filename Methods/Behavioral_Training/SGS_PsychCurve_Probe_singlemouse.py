import h5py
import re
import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit # for simoid function
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker


# =============================================================================
# SELECT MOUSE TO ANALYZE by uncommenting its respective string pattern
# =============================================================================
# string_pattern = "*mouse0691*.h5"
string_pattern = "*mouse0070*.h5"

directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim/PsychCurve/probe'
home_directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim'

file_list = os.listdir(directory)

#SELECT THE RANGE OF TRIALS THAT YOU WOULD LIKE TO CONSIDER FOR EACH DATAFRAME IN 'data'
# Define the desired slicing ranges for each dataframe

# Extract 'mouse0064' from the 'string_pattern' using regular expressions
match = re.search(r'mouse(\d+)', string_pattern)
mouse_number = match.group(1)


file_paths = glob.glob(os.path.join(directory, string_pattern))
file_paths = [path.replace("\\", "/") for path in file_paths]

print(file_paths)


# Create a list of dataframes using the various file_paths.
data = []

for fp in file_paths:
    with h5py.File(fp, 'r') as h5_file:
        dataset = h5_file['Trials']
        data_frame = pd.DataFrame(dataset[:]) 
        data.append(data_frame)


# =============================================================================
# ORGANIZE DATA
# =============================================================================

# 1. Slice data
# 2. Remove trials with pavlov = 1
# 3. Concatenate all trials
# 4. Organize data into list of dataframes ordered by the parcentage of maximum power

# =============================================================================
# 1. SLICE DATA ACCORDING TO THE TRIAL RANGES DEFINED
# =============================================================================

slicing_ranges = []

if mouse_number == '0691':
    slicing_ranges = [(0,370),  #sess01    231209
                      (0,350),  #sess02    231210
                      (0,350),  #sess03    231211
                      (0,90)]   #sess04    231212

if mouse_number == '0070':
    slicing_ranges = [(0,355),  #sess01    231209
                      (0,350),  #sess02    231210
                      (0,250),  #sess03    231211
                      (0,350)]  #sess04    231212
    

# Apply the slicing to each dataframe in the 'data' list
for i, ranges in enumerate(slicing_ranges):
    if len(ranges) == 2:
        start_idx, end_idx = ranges
        data[i] = data[i].loc[start_idx:end_idx]
    elif len(ranges) == 4:
        start_idx1, end_idx1, start_idx2, end_idx2 = ranges
        data[i] = pd.concat([data[i].loc[start_idx1:end_idx1], data[i].loc[start_idx2:end_idx2]]) 
        
# =============================================================================
# 2. REMOVE TRIALS WHERE PAVLOV = 1
# =============================================================================

for df in data:
    df.drop(df[df['pavlov'] == 1].index, inplace=True) # inplace=True : ensures that changes are made directly in the original dataframes. 
    df.drop(df[df['trialtype'] == 1].index, inplace=True) # inplace=True : ensures that changes are made directly in the original dataframes. 

# =============================================================================
# 3. CONCATENATE ALL TRIALS IN ALL DF's LISTED IN 'data'
# =============================================================================

all_data = pd.concat(data, ignore_index=True)

# =============================================================================
# 4. Organize data into list of dataframes ordered by the parcentage of maximum power
# =============================================================================

# Sort the all_data dataframe by 'amplitude_2' column in increasing order
all_data_sort = all_data.sort_values(by='amplitude_2')

# Initialize an empty list to store the dataframes
AOM_data = []

# Iterate through unique values of 'amplitude_2'
unique_amplitudes = all_data_sort['amplitude_2'].unique()

for amplitude_value in unique_amplitudes:
    # Slice the dataframe based on the current amplitude_value
    sliced_df = all_data_sort[all_data_sort['amplitude_2'] == amplitude_value]
    # Append the sliced dataframe to the AOM_data list
    AOM_data.append(sliced_df)
    
    
# =============================================================================
# DEFINING VARIABLES FROM DATA 
# =============================================================================

# trialtype is a list of arrays, where each array represents the "trialtype" column from each dataframe in the "data" list
trialtype = [df['trialtype'].values for df in AOM_data]
# iterates over each dataframe in the "data" list. For each dataframe, it 
# accesses the 'trialtype' column using df['trialtype'] and retrieves its values 
# using the .values attribute.

result = [df['result'].values for df in AOM_data]
pavlov = [df['pavlov'].values for df in AOM_data]
n_powers = len(AOM_data) # number of AOM settings presented


n_stim = []
for i in range(len(AOM_data)):
    n_stim_level = len(AOM_data[i])
    n_stim.append(n_stim_level)
    
stim_present = pd.DataFrame({'AOM' : unique_amplitudes, 'times_presented': n_stim})


# =============================================================================
# TABULATE SUCCESS RATES OF VARIOUS POWER LEVELS
# =============================================================================

correct_trials = []

for i in range(n_powers):
    correct_trials.append([])  # Initialize a list for this power level
    for x in range(n_stim[i]):
        if result[i][x] == 2:
            correct_trials[i].append(1)
        if result[i][x] != 2:
            correct_trials[i].append(0)


sum_correct_trials = []  # Create a list to store the sums of each inner list

for i in range(n_powers):
    sublist_sum = sum(correct_trials[i])  
    sum_correct_trials.append(sublist_sum)   # Append the sum to the 'sums' list
    
# FRACTION OF CORRECT TRIALS

frac_correct = []  # Create a list to store the division results

for x, y in zip(sum_correct_trials, n_stim):
    if y != 0:
        frac_correct.append(x / y)
    else:
        # Handle division by zero or any other appropriate action
        frac_correct.append(float('inf'))  # Using infinity as an example

# stim_present = pd.DataFrame({'AOM' : unique_amplitudes, 'n_presented': n_stim, 'n_correct': sum_correct_trials, 'fraction_correct': frac_correct})


# =============================================================================
# POWER CALLIBRATION FOR AOM AMPLITUDES
# =============================================================================

# =============================================================================
# # Calibration percentages for each amplitude
# calibration_dict = {
#     3424: 100,
#     3179: 95,
#     2845: 90,
#     2642: 85,
#     2485: 80,
#     2349: 75,
#     2224: 70,
#     2100: 65,
#     1986: 60,
#     1874: 55,
#     1764: 50,
#     1655: 45,
#     1543: 40,
#     1428: 35,
#     1312: 30,
#     1189: 25,
#     1054: 20,
#     905: 15,
#     731: 10,
#     501: 5,
# }
# 
# =============================================================================


cal_power = []

for amp_array in unique_amplitudes:
    amp_percent = []
    for digit_char in str(amp_array):
        digit = int(digit_char)
        if amp_array == 3424:
            amp_percent.append(round((100/100)*58, 2))
        elif amp_array == 3179:
            amp_percent.append(round((95/100)*58, 2))
        elif amp_array == 2845:
            amp_percent.append(round((90/100)*58, 2))
        elif amp_array == 2642:
            amp_percent.append(round((85/100)*58, 2))
        elif amp_array == 2485:
            amp_percent.append(round((80/100)*58, 2))
        elif amp_array == 2349:
            amp_percent.append(round((75/100)*58, 2))
        elif amp_array == 2224:
            amp_percent.append(round((70/100)*58, 2))
        elif amp_array == 2100:
            amp_percent.append(round((65/100)*58, 2))
        elif amp_array == 1986:
            amp_percent.append(round((60/100)*58, 2))
        elif amp_array == 1874:
            amp_percent.append(round((55/100)*58, 2))
        elif amp_array == 1764:
            amp_percent.append(round((50/100)*58, 2))
        elif amp_array == 1655:
            amp_percent.append(round((45/100)*58, 2))
        elif amp_array == 1543:
            amp_percent.append(round((40/100)*58, 2))
        elif amp_array == 1428:
            amp_percent.append(round((35/100)*58, 2))
        elif amp_array == 1312:
            amp_percent.append(round((30/100)*58, 2))
        elif amp_array == 1189:
            amp_percent.append(round((25/100)*58, 2))
        elif amp_array == 1054:
            amp_percent.append(round((20/100)*58, 2))
        elif amp_array == 905:
            amp_percent.append(round((15/100)*58, 2))
        elif amp_array == 731:
            amp_percent.append(round((10/100)*58, 2))
        elif amp_array == 501:
            amp_percent.append(round((5/100)*58, 2))
        else:
            amp_percent.append(round((0/100)*58, 2))

    cal_power.append(amp_percent)


# Flatten the list of arrays
conc_cal_power = np.concatenate(cal_power)

# Get the unique values
power_levels = np.unique(conc_cal_power)


# =============================================================================
# RESULTS
# =============================================================================


SUMMARY = pd.DataFrame({
    'power_levels': power_levels,
    'n_stim': n_stim,
    'n_correct': sum_correct_trials,
    'frac_correct': frac_correct
})

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
for p, N in zip(frac_correct, n_stim):
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
for p, ci in zip(frac_correct, CI):
    # Calculate the upper and lower bounds of the confidence interval
    upper_bound = p + ci
    lower_bound = p - ci
    
    # Store the upper and lower bounds as an array [upper_bound, lower_bound]
    confidence_intervals.append([upper_bound, lower_bound])

# Convert confidence_intervals to a NumPy array
confidence_intervals = np.array(confidence_intervals) 


# =============================================================================
# DEFINE SIGMOID FUNCTION
# =============================================================================

# =============================================================================
# # Define the sigmoid function: 
# # y = L / (1 + np.exp(-k * (x - x0)))
# # L = maximum value or plateau that the function reaches.
# # x0 = midpoint or inflection point along the x-axis.
# # k = slope or steepness of the curve.
# def sigmoid(x, L, x0, k):
#     y = L / (1 + np.exp(-k * (x - x0)))
#     return y
# 
# # Define input data
# x = power_levels
# y_value = frac_correct
# 
# # curve_fit models function to data
# params, covariance = curve_fit(sigmoid, x, y_value, p0=(1,35,1)) 
# L, x0, k = params
# 
# # Create a range of x values for the sigmoid curve
# x_fit = np.linspace(min(x), max(x), 100) # 1000 = numb of points I'd like to specify in the given region
# y_fit = sigmoid(x_fit, *params)
# 
# =============================================================================

# =============================================================================
# PLOT SHOWING THE FRACTION OF SUCCESSFUL TRIALS 
# =============================================================================
# Create a single figure with two vertically stacked subplots
fig, axs = plt.subplots(2, 1, figsize=(9, 9))

x = power_levels

# Plot data in the first subplot (axs[0])
axs[0].plot(x, frac_correct, 'o-', color='red', label='Data')

# Extract the upper and lower bounds from the confidence_intervals
lower_bounds = [ci[1] for ci in confidence_intervals]
upper_bounds = [ci[0] for ci in confidence_intervals]

# Plot the confidence interval as a shaded region
axs[0].fill_between(x, lower_bounds, upper_bounds, color='red', alpha=0.2, label='95% Confidence Interval')

axs[0].axhline(y=0.0, color='gray', linestyle='--', linewidth=0.8)
axs[0].axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)
axs[0].axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8)
axs[0].set_title(f'Psychometric Curve for Single Glomerular Stimulation | Mouse {mouse_number}', fontsize=16)
axs[0].set_ylabel('Fraction of Successful Detection', fontsize=14)
axs[0].set_ylim(-0.1, 1.1)
axs[0].set_xlim(min(x)-5, max(x)+5)
axs[0].set_xticks(range(0, 61, 5))  # Set x-ticks from 0 to 60 in increments of 5
axs[0].set_xticklabels(range(0, 61, 5))  # Label the ticks
axs[0].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

# Plot an 'n_stim' bar plot with adjusted x-coordinates in the second subplot (axs[1])
bar_width = 2
axs[1].bar(x[:-1], n_stim[:-1], width=bar_width, align='center', label='Number of Light Trials', color='red')
axs[1].set_xlabel('Power Level [mW/mm²]', fontsize=14)
axs[1].set_ylabel('Number of Trials', fontsize=14)
axs[1].set_xlim(min(x)-5, max(x)+5)
axs[1].set_xticks(range(0, 61, 5))  # Set x-ticks from 0 to 60 in increments of 5
axs[1].set_xticklabels(range(0, 61, 5))  # Label the ticks
axs[1].yaxis.set_major_locator(MaxNLocator(integer=True))
axs[1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
axs[1].set_ylim(0, 14)
axs[1].set_yticks(range(0, 14, 2)) 

# Adjust layout and display the plot
plt.tight_layout()
plt.show()

# =============================================================================
# SAVE PRODUCED PLOTS
# =============================================================================

# Save the figure as a .png file in the specified directory
save_directory = os.path.join(home_directory, 'plots_behavior')
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# Define the full path for the plot
save_path = os.path.join(save_directory, f'SGS_PsychCurve_Probe_mouse{mouse_number}.png')

# Save the plot
fig.savefig(save_path, dpi=300)