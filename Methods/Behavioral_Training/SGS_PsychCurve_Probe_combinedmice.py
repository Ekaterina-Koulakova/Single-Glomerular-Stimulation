import h5py
import re
import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


string_patterns = ["*mouse0691*.h5","*mouse0070*.h5"]

frac_correct_plot = []
confidence_intervals_all = []
mouse = []

for string_pattern in string_patterns:

    directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim/PsychCurve/probe'
    home_directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim'
    
    file_list = os.listdir(directory)
    
    #SELECT THE RANGE OF TRIALS THAT YOU WOULD LIKE TO CONSIDER FOR EACH DATAFRAME IN 'data'
    # Define the desired slicing ranges for each dataframe
    
    # Extract 'mouse0064' from the 'string_pattern' using regular expressions
    match = re.search(r'mouse(\d+)', string_pattern)
    mouse_number = match.group(1)
    mouse.append(mouse_number)
    
    file_paths = glob.glob(os.path.join(directory, string_pattern))
    file_paths = [path.replace("\\", "/") for path in file_paths]
    
    
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
    
    stim_present = pd.DataFrame({'AOM' : unique_amplitudes, 'n_presented': n_stim, 'n_correct': sum_correct_trials, 'fraction_correct': frac_correct})
    
    frac_correct_plot.append(frac_correct)
    
    # =============================================================================
    # POWER CALLIBRATION FOR AOM AMPLITUDES
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

    k = 1.96  # For 95% confidence
    
    CI = []
    
    for p, N in zip(frac_correct, n_stim):
        SE = np.sqrt(p * (1 - p) / N) if N > 0 else 0
        CI_value = SE * k
        CI.append(CI_value)
    
    # Store confidence intervals for this mouse
    confidence_intervals = [[p + ci, p - ci] for p, ci in zip(frac_correct, CI)]
    confidence_intervals_all.append(np.array(confidence_intervals))
    

# =============================================================================
# PLOT SHOWING THE FRACTION OF SUCCESSFUL TRIALS 
# =============================================================================
# Create a single figure with two vertically stacked subplots

fig, ax = plt.subplots(1, 1, figsize=(9, 5))
x = power_levels

for i, sublist in enumerate(frac_correct_plot):
    # Plot the frac_correct for each mouse
    line, = ax.plot(x, sublist, 'o-', label=f'Mouse {mouse[i]}')
    line_color = line.get_color()
    
    # Plot each mouse's confidence intervals
    confidence_intervals = confidence_intervals_all[i]
    lower_bounds = confidence_intervals[:, 1]
    upper_bounds = confidence_intervals[:, 0]
    
    ax.fill_between(x, lower_bounds, upper_bounds, alpha=0.2, color=line_color)

ax.axhline(y=0.0, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8)

# Set title and axis labels
ax.set_title(f'Psychometric Curve for Single Glomerular Stimulation | Probe Trials', fontsize=16)
ax.set_ylabel('Fraction of Successful Detection', fontsize=14)
ax.set_xlabel('Power Density [mW/mm$^2$]', fontsize=14)
ax.set_ylim(-0.1, 1.1)
ax.set_xlim(min(x) - 5, max(x) + 5)

# Set x-ticks from 0 to 60 in increments of 5
ax.set_xticks(range(0, 61, 5))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

# Add a legend
ax.legend(loc='upper left', fontsize='large')

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
save_path = os.path.join(save_directory, f'SGS_PsychCurve_Probe_combinedmice.png')

# Save the plot
fig.savefig(save_path, dpi=300)