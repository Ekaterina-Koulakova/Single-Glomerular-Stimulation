import h5py
import re
import os
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


string_patterns = ["*mouse0070*.h5"]

frac_correct_plot = []
mouse = []

for string_pattern in string_patterns:

    directory = 'C:/Users/Ekaterina/Dropbox/NYU/SingleGlomStim/Behavior_SingleGlomStim/Psych Curve/AOM_block'
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
    # 3. Organize data by stimulation amplitude blocks
    # 4. Concatenate dataframs with the same amplitude
    
    # =============================================================================
    # 1. SLICE DATA ACCORDING TO THE TRIAL RANGES DEFINED
    # =============================================================================
# DEFINE SESSION TYPE: Random (0) or Graded (1)
  
slicing_ranges = []    

if mouse_number == '0070':
    slicing_ranges = [(20,300),  #sess01    240319 
                      (20,249),  #sess02    240320
                      (20,300),  #sess03    240321
                      (20,248),  #sess04    240322
                      (20,300),  #sess05    240325
                      (20,300),  #sess06    240326
                      (20,248)]  #sess07    240327
                      
    
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
    
    # =============================================================================
    # 3. ORGANIZE DATA TO CREATE TRIAL BLOCKS BASED ON STIM AMPLITUDE USED 
    # =============================================================================
    
    amplitude_2_unique_values = []
    slice_on = []
    
    for df in data:
        # Extract unique values from the 'amplitude_2' column and append to 'amplitude_2_unique_values' list
        unique_values = df['amplitude_2'].unique()
        amplitude_2_unique_values.append(unique_values)
        
        # Initialize a dictionary to track the first occurance of each value.
        first_occurance ={}
        
        # Iterate through each value in unique_value
        for value in unique_values:
            # Find the index of the first occurance of the value in the dataframe
            idx = df[df['amplitude_2'] == value].index[0]
            #Store the index in the first_occurance dictionary.
            first_occurance[value] = idx
            
        # Convert the dictionary to a list and appeand in to 'slice_on'
        slice_on.append(list(first_occurance.values()))
        
    # 'amplitude_2_unique_values' contains arrays of unique values from the 'amplitude_2' column of each dataframe in 'data'
    # 'slice_on' contains arrays of the first occurrence index of each unique value in 'amplitude_2' for each dataframe in 'data'
       
    # Initialize an empty list to store slice_off
    slice_off = slice_on
    
    # Remove the first element from each sublist in slice_off
    slice_off = [sublist[1:] for sublist in slice_off]
    
    # Subtract 1 from each member of slice_off
    slice_off = [[value - 1 for value in sublist] for sublist in slice_off]

    for sublist in slice_off:
        sublist.append(len(data[slice_off.index(sublist)]))
    
    # Subtract 20 from each value in slice_on and slice_off
    adjusted_slice_on = [[val - 20 for val in sublist] for sublist in slice_on]
    adjusted_slice_off = [[val - 20 for val in sublist] for sublist in slice_off]

    # Create block_ranges based on adjusted_slice_on and adjusted_slice_off
    block_ranges = []
    for i in range(len(adjusted_slice_on)):
        ranges = [(adjusted_slice_on[i][j], adjusted_slice_off[i][j]) for j in range(len(adjusted_slice_on[i]))]
        block_ranges.append(ranges)

    data_sliced = []
    for df, ranges in zip(data, block_ranges):
        sliced_dfs = [df.iloc[start:end] for start, end in ranges]
        data_sliced.extend(sliced_dfs)
    
    data_sliced.sort(key=lambda df: df['amplitude_2'].iloc[0])
    
    # =============================================================================
    # 4. CONCATENATE ALL DFs WITH EQUAL STIM 'amplitude_2' VALUES 
    # =============================================================================
    
    data_org = []
    i = 0
    while i < len(data_sliced):
        curr_df = data_sliced[i]
        curr_value = curr_df['amplitude_2'].iloc[0]
        j = i + 1
        while j < len(data_sliced) and data_sliced[j]['amplitude_2'].iloc[0] == curr_value:
            curr_df = pd.concat([curr_df, data_sliced[j]], ignore_index=True)
            j += 1
        data_org.append(curr_df)
        i = j
        
     
    # =============================================================================
    # DEFINING VARIABLES FROM DATA 
    # =============================================================================
    
    # trialtype is a list of arrays, where each array represents the "trialtype" column from each dataframe in the "data" list
    trialtype = [df['trialtype'].values for df in data_org]
    # iterates over each dataframe in the "data" list. For each dataframe, it 
    # accesses the 'trialtype' column using df['trialtype'] and retrieves its values 
    # using the .values attribute.
    
    result = [df['result'].values for df in data_org]
    
    n_stim = []
    
    for arr in trialtype:
        count_ones = sum(1 for val in arr if val == 1)
        n_stim.append(count_ones)
        
    len_data_org = []
    
    for df in data_org:
        length = len(df)
        len_data_org.append(length)
        
        
    unique_amplitudes = []

    # Iterate through each dataframe in 'data_org'
    for df in data_org:
        # Extract unique values from the 'amplitude_2' column
        unique_values = df['amplitude_2'].unique().tolist()
        # Extend the 'unique_amplitudes' list with unique values from this dataframe
        unique_amplitudes.extend(unique_values)

    # Remove duplicates by converting to set and back to list
    unique_amplitudes = list(set(unique_amplitudes))

    # Sort the list in ascending order
    unique_amplitudes.sort()
        
    
    # =============================================================================
    # TABULATE SUCCESS RATES OF VARIOUS POWER LEVELS
    # =============================================================================
    
    correct_trials = []
    
    for i in range(len(data_org)):
        correct_trials.append([])  # Initialize a list for this dataframe
        for x in range(len(trialtype[i])):
            if trialtype[i][x] == 0 and result[i][x] == 2:
                correct_trials[i].append(1)  # Correct trial
            elif trialtype[i][x] == 1 and result[i][x] == 1:
                correct_trials[i].append(1)  # Correct trial
            else:
                correct_trials[i].append(0)  # Incorrect trial

    sum_correct_trials = []  # Create a list to store the sums of each inner list
    
    for i in range(len(data_org)):
        sublist_sum = sum(correct_trials[i])  
        sum_correct_trials.append(sublist_sum)   # Append the sum to the 'sums' list
        
    sum_trials = []
        
    for i in range(len(data_org)):
        sum_trials.append(len(data_org[i]))
        
        
    # FRACTION OF CORRECT TRIALS
    
    frac_correct = []

    # Calculate the fraction of correct trials for each dataframe
    for sum_correct, sum_trials in zip(sum_correct_trials, sum_trials):
        if sum_trials != 0:
            frac_correct.append(sum_correct / sum_trials)
        else:
            # Handle division by zero or any other appropriate action
            frac_correct.append(float('inf')) 
    
    summary = pd.DataFrame({'AOM' : unique_amplitudes, 
                            'fraction_success': frac_correct, 
                            'times_presented': n_stim})
    
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
#     ERROR CALCULATION 
# =============================================================================
    
    error = []
    confidence_interval = []

    # Iterate over the arrays in frac_SUCCESS_stim and n_stim_sess
    for frac_array, n_array in zip(frac_correct, len_data_org):
        # Calculate Error for each corresponding array
        error_array = np.sqrt((frac_array * (1 - frac_array)) / n_array)
        corr = 1.963 * error_array
    
        # Append the calculated array to the Error list
        error.append(error_array)
        confidence_interval.append(corr)
    
    
    # =============================================================================
    # RESULTS
    # =============================================================================
    
    
    SUMMARY = pd.DataFrame({
        'AOM' : unique_amplitudes, 
        'power_levels': power_levels,
        'n_stim': n_stim,
        'frac_correct': frac_correct
    })

# =============================================================================
# PLOT SHOWING THE FRACTION OF SUCCESSFUL TRIALS 
# =============================================================================
# Create a single figure with two vertically stacked subplots

fig, ax = plt.subplots(1, 1, figsize=(9,5))

x = power_levels
y = np.array(frac_correct)
confidence_interval = np.array(confidence_interval)

ax.plot(x, y, 'o-', color='red', label=f'Mouse 0070', alpha=0.7)

ax.fill_between(x, y - confidence_interval, y + confidence_interval, color='red', alpha=0.2)

ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8)

ax.set_title(f'Psychometric Curve for Single Glomerular Stimulation', fontsize=16)
ax.set_ylabel('Fraction of Successful Detection', fontsize=14)
ax.set_xlabel('Power Density [mW/mm$^2$]', fontsize=14)
ax.set_ylim(0.4, 1.1)
ax.set_xlim(min(x)-5, max(x)+5)
ax.set_xticks(range(0, 61, 5))  # Set x-ticks from 0 to 60 in increments of 5
ax.set_xticklabels(range(0, 61, 5))  # Label the ticks
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax.legend(loc='lower right', fontsize='large')

# Adjust layout and display the plot
plt.show()

# =============================================================================
# SAVE PRODUCED PLOTS
# =============================================================================

# Save the figure as a .png file in the specified directory
save_directory = os.path.join(home_directory, 'plots_behavior')
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

# Define the full path for the plot
save_path = os.path.join(save_directory, f'SGS_PsychCurve_combinedmice.png')

# Save the plot
fig.savefig(save_path, dpi=300)