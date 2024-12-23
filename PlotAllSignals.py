# -*- coding: utf-8 -*-
"""
Created on Fri May  3 18:50:03 2024

@author: mkolisnyk
"""


# calculate statistics for CCA, CA and ISO.
from scipy.io import loadmat,savemat
import neurokit2 as nk
import os
import numpy as np
import mne
import matplotlib.pyplot as plt
from mne_connectivity import spectral_connectivity_epochs
import pandas as pd
import seaborn as sns
from scipy.stats import zscore
from scipy.signal import find_peaks
    
# Define frequency bands
fmin = [1, 4, 8, 12, 30]
fmax = [4, 8, 12, 30, 50]
ecg_idx = [6,7]; nc_idx = [5]; eeg_idx = [0,1,2,3,4]; tcd_idx = [8]; art_idx = [9]
fs = 100
# Window sizes in seconds
window_sizes = [1, 10, 30, 60]

crit_times = [(4000,8000),(13000,17000),(4000,6000)]


 # Set the global font size, figure size, and line width for plots
plt.rcParams.update({
    'font.size': 8,  # Adjust based on your journal's requirements
    #'figure.figsize': (12.68, 6.34),  # Typical figure size, but adjust as needed
    'figure.dpi': 600,
    'lines.linewidth': 0.5,
    'axes.labelsize': 8,  # Adjust label size as necessary
    'font.family': 'Arial',
    'axes.titlesize': 8,
    'xtick.labelsize':8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.autolayout': True,  # Adjust layout automatically
    'axes.spines.top': False,  # Remove top and right spines
    'axes.spines.right': False,
    'axes.grid': True,  # Enable grid (optional, based on your preference)
    'grid.color': 'gray',  # Grid color
    'grid.linestyle': 'None',  # Grid linestyle
    'grid.linewidth': 0.5,
    'grid.alpha': 0.5  # Grid line transparency
    
})


def compute_zero_durations(signal, sfreq, window_sizes):
    """
    Computes the duration of zero values within specified window sizes for a given signal.

    :param signal: numpy array, the signal data (e.g., TCD signal)
    :param sfreq: float, sampling frequency in Hz
    :param window_sizes: list of int, window sizes in seconds to analyze
    :return: dict, keys are window sizes and values are lists of durations (in seconds) of zeros
    """
    # Ensure window_sizes is a list of individual sizes, not a list of lists
    window_sizes = np.array(window_sizes).flatten()

    # Convert window sizes from seconds to samples
    window_sizes_samples = [int(size * sfreq) for size in window_sizes]

    # Dictionary to store the results for each window size
    durations_per_window_size = {}

    # Iterate over each window size
    for window_size in window_sizes_samples:
        durations = []  # List to store durations of zeros for this window size

        # Slide window over the signal
        for start in range(0, len(signal) - window_size + 1, window_size):
            window = signal[start:start + window_size]
            
            # Count zeros and set to a proportion
            zero_count = np.count_nonzero(window != 0)/window.size
            
            durations.append(zero_count)

        # Store the durations list in the dictionary under the key of the current window size (in seconds)
        durations_per_window_size[window_size / sfreq] = durations

    return durations_per_window_size

from scipy.signal import find_peaks, butter, filtfilt

def butter_lowpass_filter(data, cutoff, fs, order=5):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def compute_pulse_pressure(art_data, sfreq, window_sizes, cutoff_freq=2.0, order=5):
    """
    Computes the pulse pressure (difference between systolic and diastolic pressure)
    over specified window sizes for a given ART data.

    :param art_data: numpy array, the ART data
    :param sfreq: float, sampling frequency in Hz
    :param window_sizes: list of int, window sizes in seconds to analyze
    :param cutoff_freq: float, cutoff frequency for low-pass filter (default 2.0 Hz)
    :param order: int, order of the low-pass filter (default 5)
    :return: dict, keys are window sizes and values are lists of pulse pressures
    """
    # Apply low-pass filter to remove high-frequency noise
    #filtered_data = butter_lowpass_filter(art_data, cutoff=cutoff_freq, fs=sfreq, order=order)
    filtered_data = art_data

    # Convert window sizes from seconds to samples
    window_sizes_samples = [int(size * sfreq) for size in window_sizes]

    # Dictionary to store the results for each window size
    pulse_pressures = {}

    # Iterate over each window size
    for window_size in window_sizes_samples:
        pulse_pressure_per_window = []

        # Slide window over the ART data
        for start in range(0, len(filtered_data) - window_size + 1, window_size):
            window = filtered_data[start:start + window_size]

            # Find peaks (systolic points)
            peaks_indices, _ = find_peaks(window, distance=sfreq*0.5)  # Minimum distance between peaks is 0.5 seconds
            troughs_indices, _ = find_peaks(-window, distance=sfreq*0.5)  # Invert data to find troughs (diastolic points)

            if len(peaks_indices) > 0 and len(troughs_indices) > 0:
                average_peaks = np.median(window[peaks_indices])
                average_troughs = np.median(window[troughs_indices])

                # Compute pulse pressure: average peak - average trough
                pulse_pressure = average_peaks - average_troughs
                pulse_pressure_per_window.append(pulse_pressure)
            else:
                # In case no peaks or troughs are found, append NaN
                pulse_pressure_per_window.append(np.nan)

        # Store the pulse pressures list in the dictionary under the key of the current window size (in seconds)
        pulse_pressures[window_size / sfreq] = pulse_pressure_per_window

    return pulse_pressures

def compute_pulse_pressure_with_peaks(art_data, sfreq, window_sizes):
    """
    Computes the pulse pressure over specified window sizes for ART data using peak detection.
    
    :param art_data: numpy array, the arterial pressure data
    :param sfreq: float, sampling frequency in Hz
    :param window_sizes: list of int, window sizes in seconds to analyze
    :return: dict, keys are window sizes and values are lists of computed pulse pressures
    """
    # Convert window sizes from seconds to samples
    window_sizes_samples = [int(size * sfreq) for size in window_sizes]

    # Dictionary to store the results for each window size
    pulse_pressures = {}

    # Iterate over each window size
    for window_size in window_sizes_samples:
        pulse_pressure_per_window = []

        # Slide window over the ART data
        for start in range(0, len(art_data) - window_size + 1, window_size):
            window = art_data[start:start + window_size]

            # Find peaks (systolic points)
            peaks_indices, _ = find_peaks(window)
            troughs_indices, _ = find_peaks(-window)  # Invert data to find troughs (diastolic points)

            # Calculate the average of the peaks and troughs
            if len(peaks_indices) > 0 and len(troughs_indices) > 0:
                average_peaks = np.median(window[peaks_indices])
                average_troughs = np.median(window[troughs_indices])

                # Compute pulse pressure: average peak - average trough
                pulse_pressure = average_peaks - average_troughs
                pulse_pressure_per_window.append(pulse_pressure)
            else:
                # In case no peaks or troughs are found, append a default value or NaN
                pulse_pressure_per_window.append(np.nan)

        # Store the pulse pressures list in the dictionary under the key of the current window size (in seconds)
        pulse_pressures[window_size / sfreq] = pulse_pressure_per_window

    return pulse_pressures

def compute_eeg_std(eeg_data, sfreq, window_sizes):
    """
    Computes the standard deviation over specified window sizes for a given EEG data.

    :param eeg_data: numpy array, the EEG data (channels x samples)
    :param sfreq: float, sampling frequency in Hz
    :param window_sizes: list of int, window sizes in seconds to analyze
    :return: dict, keys are window sizes and values are lists of standard deviations for each window
    """
    # Convert window sizes from seconds to samples
    window_sizes_samples = [int(size * sfreq) for size in window_sizes]

    # Dictionary to store the results for each window size
    std_per_window_size = {}

    # Iterate over each window size
    for window_size in window_sizes_samples:
        std_per_window = []

        # Slide window over the EEG data
        for start in range(0, eeg_data.shape[0] - window_size + 1, window_size):
            window = eeg_data[start:start + window_size]

            # Calculate mean and standard deviation
            mean = np.mean(window)
            std_dev = np.std(window)
    
            # Identify outliers: elements more than 4 standard deviations from the mean
            outliers = np.abs(window - mean) > 4 * std_dev
    
            # Replace outliers with NaN
            window_cleaned = np.where(outliers, np.nan, window)
    
            # Recalculate standard deviation with outliers removed
            std_dev_cleaned = np.nanstd(window_cleaned)
    
            std_per_window.append(std_dev_cleaned)

        # Store the cleaned standard deviations list in the dictionary under the key of the current window size (in seconds)
        std_per_window_size[window_size / sfreq] = std_per_window

    return std_per_window_size


full_name = ['Loss of brain blood flow','Loss of brain activity','Loss of systemic circulation']

# set file locations
patient_files = ['','','']
patient_ids = [2,3,4]
root = r"D:/MATTFILES/TCD_Clean/"

for pat_idx,file in enumerate(patient_files): 
    
    current_file = root+file
    # change directory 
    os.chdir(current_file)
    
    # load file
    Signals = loadmat("TCDClean.mat")
    
    EEG = Signals['TCDClean']['EEG'][0][0]
    ECG = Signals['TCDClean']['EEGECG'][0][0]
    TCD = Signals['TCDClean']['V'][0][0].T
    ART = Signals['TCDClean']['BP'][0][0].T
    
    # reduce to the smallest array
    eeg_shape = EEG.shape[0];  ecg_shape = ECG.shape[0]; 
    tcd_shape = TCD.shape[0];  art_shape = ART.shape[0]
    
    trim = np.min([eeg_shape,ecg_shape,tcd_shape,art_shape])
    EEG = EEG[0:trim,:]; ECG = ECG[0:trim,:]; 
    TCD = TCD[0:trim,:]; ART = ART[0:trim,:]; 
    
    
    # drop the clear outliers from the signal
    for chan in range(0,EEG.shape[1]):
        # get channel mean
        chan_mean = np.nanmean(EEG[:,chan])
        # replace outliers with mean
        EEG[(np.abs(zscore(EEG[:,chan],nan_policy = 'omit')) > 5),chan] = chan_mean
        
        # do it again
        chan_mean = np.nanmean(EEG[:,chan])
        EEG[(np.abs(zscore(EEG[:,chan],nan_policy = 'omit')) > 5),chan] = chan_mean
        
        
    # same with ECG 
    # drop the clear outliers from the signal
    for chan in range(0,ECG.shape[1]):
        # get channel mean
        chan_mean = np.nanmean(ECG[:,chan])
        # replace outliers with mean
        ECG[(np.abs(zscore(ECG[:,chan],nan_policy = 'omit')) > 5),chan] = chan_mean
        
        # do it again
        chan_mean = np.nanmean(ECG[:,chan])
        ECG[(np.abs(zscore(ECG[:,chan],nan_policy = 'omit')) > 5),chan] = chan_mean
    
    # apply some filtering
    #TCD = mne.filter.filter_data(TCD,sfreq = 100, l_freq = None, h_freq = 10,n_jobs = -1) # see if you could smooth out the signal for plotting
    #ART = mne.filter.filter_data(ART,sfreq = 100, l_freq = None, h_freq = 10,n_jobs = -1) # see if you could smooth out the signal for plotting
    
    
    info = mne.create_info(ch_names = ['F3-P3','F4-P4','Fz-Pz','F7-T5','F8-T6','A1-A2','ECGL','ECGR','TCD','ART'],sfreq = 100, 
                           ch_types = ['eeg','eeg','eeg','eeg','eeg','eeg','ecg','ecg','bio','bio'])
    
  
    Raw = mne.io.RawArray(np.concatenate([EEG,ECG,TCD,ART],1).T, info = info)
    # create event structure

    
    events = np.zeros((3,3),dtype = 'int')
    events[0,0] = Signals['TCDClean']['CCA_Matt'][0][0][0][0]; events[0,2] = 0;
    events[1,0] = Signals['TCDClean']['ISO'][0][0][0][0]; events[1,2] = 1;
    events[2,0] = Signals['TCDClean']['PP'][0][0][0][0] + 1; events[2,2] = 2;
    
    event_id = {'CCA':0,'CEA':1,'CA':2}
    
    
    # Assuming 'Raw' is already created and 'events' array is set up
    
    # Convert event sample indices to times in seconds
    event_times = events[:, 0] / 100  # Convert from samples to seconds (sfreq=100 Hz
    # Plotting
    fig, axs = plt.subplots(3, 1, figsize=(6.28, 3.14), sharex=True)
    
    # Assume you have already trimmed the data arrays to 'trim'
    times = np.arange(crit_times[pat_idx][0], crit_times[pat_idx][1], 1 / 100)  # Total duration of the trimmed data in seconds
    samples = np.arange(crit_times[pat_idx][0] * fs, crit_times[pat_idx][1] * fs, 1)
    
    #TCD_win = TCD[samples]
    #ART_win = ART[samples]
    
    # Plot each modality on a different subplot
    axs[0].plot(times, TCD[samples], label='TCD', color='blue')
    axs[1].plot(times, ART[samples], label='ART', color='red')
    axs[2].plot(times, np.mean(EEG[samples,0:5], axis=1) * 10e5, label='GFP', color='k')  # Assuming you want to plot mean EEG
    #axs[3].plot(times, ECG[samples, 0], label='ECG', color='red')  # Assuming first column is the ECG channel of interest
    
    
    colour = ['b','k','r']
    # Adding event lines and labels
    for ax in axs:
        ylim = ax.get_ylim()[1]  # Get the current y-limit of the plot
        for event_idx, event_time in enumerate(event_times):
            ax.axvline(x=event_time, color=colour[event_idx], linestyle='--', linewidth=1)
 
    # Set labels and titles
    axs[0].set_ylabel('BBFv (cm/s)')
    axs[1].set_ylabel('ABP (mmHg)')
    axs[2].set_ylabel('EEG (µV)')
    #axs[3].set_ylabel('ECG Signal (µV)')
    axs[2].set_xlabel('Time from recording start (seconds)')
    

    plt.tight_layout()
    plt.show()
      
    # lets create a similar plot but by aggregrating this information over a window


    # Assuming ZEROFLOW, PULSEPRESSURE, and EEGVAR are dictionaries with keys as window sizes and values as lists of computed value
    
    # compute proportion of time with zero flow
    ZEROFLOW = compute_zero_durations(TCD[samples],100,[1,10,30,60])
    
    PULSEPRESSURE = compute_pulse_pressure(ART[samples].reshape(-1,), 100, [1,10,30,60])
    
    PPP = compute_pulse_pressure_with_peaks(ART[samples].reshape(-1,), 100, [1,10,30,60])
    EEGVAR = compute_eeg_std(np.median(EEG[samples,0:5],1,keepdims=True), 100, [1,10,30,60])
    

    
    # Compute time vector for each window size
    time_vectors = {size: np.arange(times[0], times[0] + len(ZEROFLOW[size]) * size, size) for size in window_sizes}
    
        # Choose a window size to plot
    window_size_to_plot = 60
    
    # Create a figure and axis object
    fig,ax = plt.subplots(nrows = 3, figsize=(6.28, 6.28))
    
    # Plot zero flow proportion
    ax[0].plot(time_vectors[window_size_to_plot], np.array(ZEROFLOW[window_size_to_plot]) * 100, label='Zero Flow Proportion', \
             marker='o', linestyle='-', color='blue')
    
    ax[0].set_ylabel("Non-Zero Flow \n (%)")
    # Plot pulse pressure
    ax[1].plot(time_vectors[window_size_to_plot], PULSEPRESSURE[window_size_to_plot], label='Pulse Pressure',\
             marker='s', linestyle='-', color='red')
    
    ax[1].set_ylabel("Pulse Pressure \n (mm/Hg)")
    # Plot pulse pressure
    #ax[2].plot(time_vectors[window_size_to_plot], PPP[window_size_to_plot], label='Pulse Pressure',\
    #         marker='s', linestyle='-', color='red')
        
    # Plot EEG standard deviation
    ax[2].plot(time_vectors[window_size_to_plot], np.array(EEGVAR[window_size_to_plot]) * 10e5, label='GFP Std Deviation',\
               marker='^', linestyle='-', color='k')
    ax[2].set_ylabel("EEG \n Standard Deviation (uV)")
    
    # Adding labels and title
    plt.xlabel('Time from recording start (seconds)')
    #plt.title(f'Combined Metrics Over Time for Window Size {window_size_to_plot}s')
    # Adding event lines and labels
    for ax in ax:
        for event_idx, event_time in enumerate(event_times):
            ax.axvline(x=event_time, color=colour[event_idx], linestyle='--', linewidth=1)
            
    #plt.legend()
    plt.grid(True)
    plt.show()
    
    
    # Compute time vector for each window size
    time_vectors = {size: np.arange(times[0], times[0] + len(ZEROFLOW[size]) * size, size) for size in window_sizes}
    
    # Choose a window size to plot
    window_size_to_plot = 60
    
    # Create a figure and axis object
    fig, axs = plt.subplots(3, 1, figsize=(6.28, 6.28), sharex=True)
    
    # Plot the original time series data
    axs[0].plot(times, TCD[samples], label='TCD', color='blue')
    axs[1].plot(times, ART[samples], label='ART', color='red')
    axs[2].plot(times, np.mean(EEG[samples, 0:5], axis=1) * 10e5, label='GFP', color='k')
    
    # Create secondary y-axes for aggregated data
    axs2 = [ax.twinx() for ax in axs]
    
    # Plot the aggregated data
    axs2[0].plot(time_vectors[window_size_to_plot], np.array(ZEROFLOW[window_size_to_plot]) * 100, label='NonZero Flow Percentage', marker='o', linestyle='-', color='cyan')
    axs2[1].plot(time_vectors[window_size_to_plot], PULSEPRESSURE[window_size_to_plot], label='Pulse Pressure', marker='s', linestyle='-', color='orange')
    axs2[2].plot(time_vectors[window_size_to_plot], np.array(EEGVAR[window_size_to_plot]) * 10e5, label='EEG Std Deviation', marker='^', linestyle='-', color='brown')
    
    # Set labels and titles
    axs[0].set_ylabel('BBFv (cm/s)')
    axs[1].set_ylabel('ABP (mmHg)')
    axs[2].set_ylabel('EEG (µV)')
    axs[2].set_xlabel('Time from recording start (seconds)')
    axs2[0].set_ylabel('NonZero Flow (%)')
    axs2[1].set_ylabel('Pulse Pressure (mmHg)')
    axs2[2].set_ylabel('EEG Std Deviation (µV)')
    
    # Add event lines
    colour = ['b', 'k', 'r']
    for ax in axs:
        ylim = ax.get_ylim()[1]  # Get the current y-limit of the plot
        for event_idx, event_time in enumerate(event_times):
            ax.axvline(x=event_time, color=colour[event_idx], linestyle='--', linewidth=1)
    
    # # Adding labels and legends
    # for ax in axs:
    #     ax.legend(loc='upper left')
    #     ax.grid(True)
    # for ax in axs2:
    #     ax.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()
    
    # sort into epochs
    RawEpochs = mne.Epochs(Raw,events,event_id, tmin = -10, tmax = 30, baseline = None,preload=True, verbose = True)
    
    # for each event
    
    
    # Calculate the time vector for the epochs
    tmin, tmax, sfreq = -10, 30, 100  # start time, end time in seconds, and sampling frequency in Hz
    times = np.linspace(tmin, tmax, int((tmax - tmin) * sfreq + 1), endpoint=True)
    
    # Now update the plotting section
    for event_name, eid in event_id.items():
        fig, axs = plt.subplots(5, 1, figsize=(3.14, 5))
        data = RawEpochs[event_name].get_data(picks=None)  # All data for this event
        
        # Assuming eeg_idx, ecg_idx, tcd_idx, art_idx are defined somewhere in your code
        eeg_data = np.median(data[:,eeg_idx,:], axis=1)  # Average EEG data across channels
        ecg_data = np.median(data[:,ecg_idx,:],axis=1)  # Average ECG data across channels
        tcd_data = data[:,tcd_idx,:].squeeze()  # Assuming only one TCD channel
        art_data = data[:,art_idx,:].squeeze()  # Assuming only one ART channel
        nc_data = data[:,nc_idx,:].squeeze()
        
        # Plot the data with time axis
        axs[0].plot(times, tcd_data.T, color='blue', label='TCD')
        axs[0].set_ylabel('BBFv (cm/s)')
        axs[0].axvline(0, linestyle = "--", color = 'grey', alpha = 0.5)
        #axs[0].legend()
        
        axs[1].plot(times, art_data.T, color='red', label='ART')
        axs[1].set_ylabel('ART (mm/Hg)')
        axs[1].axvline(0, linestyle = "--", color = 'grey', alpha = 0.5)
        #axs[1].legend()
        
        axs[2].plot(times, eeg_data.T * 1e6, color='k', label='EEG')  # Adjust scaling if necessary
        axs[2].set_ylabel('EEG (μV)')
        axs[2].axvline(0, linestyle = "--", color = 'grey', alpha = 0.5)
        #axs[2].legend()
        axs[2].set_ylim([-10, 10])
        
        axs[3].plot(times, ecg_data.T * 1e6, color='g', label='ECG')  # Adjust scaling if necessary
        axs[3].set_ylabel('ECG (μV)')
        axs[3].axvline(0, linestyle = "--", color = 'grey', alpha = 0.5)
        axs[3].set_xlabel('Time (seconds)')
        #axs[3].legend()
        
        axs[4].plot(times, nc_data.T * 1e6, color='g', label='Non-Cephalic EEG')  # Adjust scaling if necessary
        axs[4].set_ylabel('Non-Cephalic \n EEG (μV)')
        axs[4].axvline(0, linestyle = "--", color = 'grey', alpha = 0.5)
        axs[4].set_xlabel('Time (seconds)')
        #axs[3].legend()
    
        axs[0].set_title(f'{full_name[eid]}')
        
        #plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room for the figure title
        plt.show()
    
        
        
    