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
# 3. Organize data into list of dataframes ordered by the parcentage of maximum power


# =============================================================================
# 1. SLICE DATA ACCORDING TO THE TRIAL RANGES DEFINED
# =============================================================================

slicing_ranges = []

if mouse_number == '0691':
    slicing_ranges = [(0,370),  #sess01    231209
                      (0,350),  #sess02    231210
                      (0,400),  #sess03    231211
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
# 3. Organize data by sessions into a list of dataframes ordered by the parcentage of maximum power
# =============================================================================

# Sort the all_data dataframe by 'amplitude_2' column in increasing order
AOM_data = []
unique_amplitudes = []

for df in data:
    sorted_df = df.sort_values(by='amplitude_2')
    AOM_data.append(sorted_df)
    unique_amp = sorted_df['amplitude_2'].unique()
    unique_amplitudes.append(unique_amp)
    

# =============================================================================
# SORT AOM DATA into SESSIONS & AOM POWER LEVELS
# =============================================================================

# Forloop that iterates over each dataframe in the AOM_data list.
# For each dataframe, it slices the data based on unique values in the 
# 'amplitude_2' column and appends those slices to the AOM_data_session list.

AOM_data_session = []

# Iterate over each dataframe in AOM_data
for df in AOM_data:
    # Get unique values from the 'amplitude_2' column
    u = df['amplitude_2'].unique()

    # Create a list to store the sliced dataframes
    session_list = []

    # Slice the dataframe based on each unique amplitude value
    for amp in u:
        session = df[df['amplitude_2'] == amp]
        session_list.append(session)

    # Append the list of sliced dataframes to AOM_data_session
    AOM_data_session.append(session_list)
    
    
len_stim = []

for session_list in AOM_data_session:
    inner_len_stim = []
    for session in session_list:
        len_session = len(session)
        inner_len_stim.append(len_session)
    len_stim.append(inner_len_stim)

# =============================================================================
# DEFINING VARIABLES FROM DATA 
# =============================================================================

trialtype = []

for session_list in AOM_data_session:
    session_trialtype = []
    for session in session_list:
        if 'trialtype' in session.columns:
            session_trialtype.append(session['trialtype'].values)
        else:
            session_trialtype.append(None)  # Append None for missing 'trialtype' column
    trialtype.append(session_trialtype)

result = []

for session_list in AOM_data_session:
    session_result = []
    for session in session_list:
        if 'result' in session.columns:
            session_result.append(session['result'].values)
        else:
            session_result.append(None)  # Append None for missing 'trialtype' column
    result.append(session_result)
    
mask_id = []

for session_list in AOM_data_session:
    session_mask_id = []
    for session in session_list:
        if 'mask_id' in session.columns:
            session_mask_id.append(session['mask_id'].values)
        else:
            session_mask_id.append(None)  # Append None for missing 'trialtype' column
    mask_id.append(session_mask_id)

pavlov = []

for session_list in AOM_data_session:
    session_pavlov = []
    for session in session_list:
        if 'pavlov' in session.columns:
            session_pavlov.append(session['pavlov'].values)
        else:
            session_pavlov.append(None)  # Append None for missing 'trialtype' column
    pavlov.append(session_pavlov)

# Number of sessions conducted

n_sessions = len(AOM_data_session)

# Array of how many AOM powers were included in each session.

n_powers = []

for amp in unique_amplitudes:
    n = len(amp)
    n_powers.append(n)

n_powers_array = np.array(n_powers)


# =============================================================================
# TABULATE SUCCESS RATES OF VARIOUS POWER LEVELS
# =============================================================================

correct_trials = []

for i in range(n_sessions):
    session_correct_trials = []  # Initialize a list for this session
    for x in range(len(result[i])):
        power_level_correct_trials = []  # Initialize a list for this power level
        for y in range(len(result[i][x])):
            if result[i][x][y] == 2:
                power_level_correct_trials.append(1)
            else:
                power_level_correct_trials.append(0)
        session_correct_trials.append(power_level_correct_trials)
    correct_trials.append(session_correct_trials)

sum_correct_trials = []

for i in range(n_sessions):
    session_sum_correct_trials = []  # Initialize a list for this session
    for x in range(len(correct_trials[i])):
        power_level_sum_correct_trials = sum(correct_trials[i][x])
        session_sum_correct_trials.append(power_level_sum_correct_trials)
    sum_correct_trials.append(session_sum_correct_trials)
    

total_trials = []

for i in range(n_sessions):
    session_total_trials = []  # Initialize a list for this session
    for x in range(len(correct_trials[i])):
        power_level_total_trials = len(correct_trials[i][x])
        session_total_trials.append(power_level_total_trials)
    total_trials.append(session_total_trials)
  
# FRACTION OF CORRECT TRIALS

frac_correct = []  # Create a list to store the division results


frac_correct = []

for i in range(n_sessions):
    session_frac_correct = []  # Initialize a list for this session
    for x in range(len(correct_trials[i])):
        if total_trials[i][x] != 0:
            power_level_frac_correct = sum(correct_trials[i][x]) / total_trials[i][x]
        else:
            power_level_frac_correct = float('inf')  # Using infinity as an example
        session_frac_correct.append(power_level_frac_correct)
    frac_correct.append(session_frac_correct)

# =============================================================================
# POWER CALLIBRATION FOR AOM AMPLITUDES
# =============================================================================

amplitude_to_power = {
    3424: (round((100/100)*58, 2)),
    3179: (round((95/100)*58, 2)),
    2845: (round((90/100)*58, 2)),
    2642: (round((85/100)*58, 2)),
    2485: (round((80/100)*58, 2)),
    2349: (round((75/100)*58, 2)),
    2224: (round((70/100)*58, 2)),
    2100: (round((65/100)*58, 2)),
    1986: (round((60/100)*58, 2)),
    1874: (round((55/100)*58, 2)),
    1764: (round((50/100)*58, 2)),
    1655: (round((45/100)*58, 2)),
    1543: (round((40/100)*58, 2)),
    1428: (round((35/100)*58, 2)),
    1312: (round((30/100)*58, 2)),
    1189: (round((25/100)*58, 2)),
    1054: (round((20/100)*58, 2)),
    905:  (round((15/100)*58, 2)),
    731:  (round((10/100)*58, 2)),
    501:  (round((5/100)*58, 2)),
}

cal_power = []

# Iterate over each list of amplitudes in unique_amplitudes
for amp_array in unique_amplitudes:
    amp_percent = []
    for amp in amp_array:  # treating amp_array as an iterable of amplitudes
        if amp in amplitude_to_power:
            amp_percent.append(amplitude_to_power[amp])  # Get the calibrated power
        else:
            amp_percent.append(round((0 / 100) * 58, 2))  # If not found, append 0 power
    cal_power.append(amp_percent)  # Append the list of calibrated powers

# Here, cal_power remains a list of lists with variable lengths
# Now if you need unique power levels, we can keep it a list of unique arrays
power_levels = [np.unique(amps) for amps in cal_power]  # Unique for each sublist


# =============================================================================
# PLOT SHOWING THE FRACTION OF SUCCESSFUL TRIALS 
# =============================================================================

# Plot the data in cal_percent and frac_correct
# cal_percent is a list of arrays and frac_correct is a list of lists


# Create a single figure
fig, ax = plt.subplots(figsize=(9, 5))

# Determine the number of lines to plot
num_lines = len(AOM_data)

# Use a colormap that transitions from purple to green
cmap = plt.get_cmap('viridis')

for i in range(num_lines):
    x = power_levels[i]
    y = frac_correct[i]

    # Convert y (a list) to a NumPy array for plotting
    y = np.array(y)

    # Get a color from the viridis colormap, where 0 is purple and 1 is green
    color = cmap(i / num_lines)

    ax.plot(x, y, 'o-', label=f'Session {i+1}', color=color)

# Add dashed lines at y = 0.0, 0.5, and 1.0
ax.axhline(y=0.0, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8)

# Set title and axis labels
ax.set_title(f'Psychometric Curve for Single Glomerular Stimulation | Mouse {mouse_number}', fontsize=16)
ax.set_ylabel('Fraction of Successful Detection', fontsize=14)
ax.set_xlabel('Power Density [mW/mm$^2$]', fontsize=14)
ax.set_ylim(-0.1, 1.1)
ax.set_xlim(-5, 65)

# Set x-ticks from 0 to 60 in increments of 5
ax.set_xticks(range(0, 61, 5))

# Add a legend
ax.legend(loc='upper left', fontsize='large')

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
save_path = os.path.join(save_directory, f'SGS_PsychCurve_Probe_mouse{mouse_number}.png')

# Save the plot
fig.savefig(save_path, dpi=300)
