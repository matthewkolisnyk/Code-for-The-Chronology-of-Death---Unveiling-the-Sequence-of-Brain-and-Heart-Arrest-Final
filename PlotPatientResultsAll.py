# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:45:10 2024

@author: mkolisnyk
"""

# plot TCD, ECG, and Phys results
import os
from scipy.io import loadmat
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 

# go to current directory
os.chdir(r"D:\MATTFILES\TCD_Clean")

# Function to update x-axis to show time in seconds or minutes
def update_xaxis_to_time(ax, start_sample, end_sample, sampling_rate=100, use_minutes=False):
    num_samples = end_sample - start_sample
    if use_minutes:
        # Convert the range to minutes
        max_time = num_samples / sampling_rate / 60  # Total duration in minutes
        step = max_time / 5  # Example step size for ticks, adjust as needed
        labels = [f"{i*step:.1f}" for i in range(6)]  # Generate labels
        ax.set_xticks(np.linspace(0, num_samples, 6))
        ax.set_xticklabels(labels)
        ax.set_xlabel('Time (minutes)')
    else:
        max_time = num_samples / sampling_rate  # Convert to seconds
        step = max_time / 5  # Example step size for ticks, adjust as needed
        labels = [f"{i*step:.1f}" for i in range(6)]  # Generate labels
        ax.set_xticks(np.linspace(0, num_samples, 6))
        ax.set_xticklabels(labels)
        ax.set_xlabel('Time (seconds)')

# Function to convert sample indices to time in minutes or seconds
def samples_to_time(sample_index, sampling_rate=100, use_minutes=False):
    if use_minutes:
        return sample_index / (sampling_rate * 60)  # Convert to minutes
    else:
        return sample_index / sampling_rate  # Convert to seconds


def convert_samples_to_plot_minutes(sample_index, start_sample, sampling_rate=100):
    # Convert both the sample index and the start sample to minutes
    sample_minutes = sample_index / sampling_rate / 60
    
    # Return the relative position in minutes
    return sample_minutes


# load AMA-typical plotting format

# Set the global font size, figure size, and line width for plots
plt.rcParams.update({
    'font.size': 12,  # Adjust based on your journal's requirements
    'figure.figsize': (12, 8),  # Typical figure size, but adjust as needed
    'figure.dpi': 400,
    'lines.linewidth': 2,
    'axes.labelsize': 12,  # Adjust label size as necessary
    'font.family': 'Arial',
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 12,
    'figure.autolayout': True,  # Adjust layout automatically
    'axes.spines.top': False,  # Remove top and right spines
    'axes.spines.right': False,
    'axes.grid': True,  # Enable grid (optional, based on your preference)
    'grid.color': 'gray',  # Grid color
    'grid.linestyle': 'None',  # Grid linestyle
    'grid.linewidth': 0.5,
    'grid.alpha': 0.5  # Grid line transparency
})

# Use a simple color palette
sns.set_palette('grey')  # For Seaborn plots, or choose a palette appropriate for color-blind readers if using color

# plot timeseries about the critical points
buffer = 3000; # 60 seconds before and after

# plot four panel version MARAT likes
clean_time = [(5000,6000),(0,1000),(15000,16000),(0,1000),(1000,2000),(590000,591000),(100000,101000)]
descent_time = [(10000,100000),(700000,800000),(1600000,1700000),(450000,550000),
                (150000,250000),(600000,700000),(650000,770000)]

# store figures 
FIGURES = []
# for each patient
patient_files = ["","","","","","",""]

# text jitter
text_jitter = [10,10,20,20,30,35,10]

for patient_idx, file in enumerate(patient_files):

    patient_file = patient_files[patient_idx]
    patient_id_anon = "Patient " + str(patient_idx + 1)
    patient_id = patient_file    
    
    # load in patient file
    data = loadmat(patient_file+"TCDClean.mat")['TCDClean']
    
    # break down MATLAB structure
    ABP = data['BP'][0][0][0].reshape(1,-1)
    V = data['V'][0][0][0].reshape(1,-1)
    ECG = data['ECG'][0][0].reshape(1,-1)
    
    # grab critical times
    CA_Time = data['PP'][0][0][0][0]
    CCA_Time = data['CCA_Matt'][0][0][0][0]
    
    # start with CA
    fig, ax = plt.subplots(nrows = 4, ncols = 1)
    
    # plot normal waves
    ax[0].plot(V[0,clean_time[patient_idx][0]:clean_time[patient_idx][1]],color = 'blue', \
               label = "Brain Blood Flow Velocity (BBFv)",alpha = 0.5); 
    ax2 = ax[0].twinx();
    #ax2.set_ylim(0,100)
    ax2.plot(ABP[0,clean_time[patient_idx][0]:clean_time[patient_idx][1]],color = 'red', \
             label = "Arterial Blood Pressure (ABP)",alpha = 0.5); 
    update_xaxis_to_time(ax[0], 0, 1000, sampling_rate=100, use_minutes=False)
    ax[0].set_ylabel("BBFv \n (cm/s)")
    ax2.set_ylabel("ABP \n (mmHg)")
    
    
    # Adding a single legend for both lines
    lines, labels = ax[0].get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='upper right',bbox_to_anchor=(1, 2.0))
    
    # add title
    ax[0].set_title(patient_id_anon + "\n", loc = 'center')
    
    # Plot descent
    start_sec = samples_to_time(descent_time[patient_idx][0], 100, use_minutes=False)
    end_sec = samples_to_time(descent_time[patient_idx][1], 100, use_minutes=False)
    ax[1].plot(V[0, descent_time[patient_idx][0]:descent_time[patient_idx][1]], color='blue',alpha = 0.5)
    ax2 = ax[1].twinx()
    #ax2.set_ylim(0,100)
    ax2.plot(ABP[0, descent_time[patient_idx][0]:descent_time[patient_idx][1]], color='red',alpha = 0.5)
    # Add vertical lines at the converted times
    ax2.axvline(CCA_Time- descent_time[patient_idx][0],0,\
                np.nanmax(ABP[0,descent_time[patient_idx][0]:descent_time[patient_idx][1]]), \
                    linestyle="--", color='black')
    ax2.axvline(CA_Time- descent_time[patient_idx][0],0,\
                np.nanmax(ABP[0,descent_time[patient_idx][0]:descent_time[patient_idx][1]]),\
                    linestyle="--", color='black')
    
    ax2.text(y=text_jitter[patient_idx] + np.nanmax(ABP[0,descent_time[patient_idx][0]:descent_time[patient_idx][1]]),  # X-coordinate for the text position
            x=CCA_Time- descent_time[patient_idx][0],  # Y-coordinate to align with the horizontal line
            s='loss of brain blood flow',  # Text content
            va='bottom',  # Vertical alignment
            ha='center',  # Horizontal alignment
            color='blue')  # Text color
    
    ax2.text(y=np.nanmax(ABP[0,descent_time[patient_idx][0]:descent_time[patient_idx][1]]),  # X-coordinate for the text position
            x=CA_Time- descent_time[patient_idx][0],  # Y-coordinate to align with the horizontal line
            s='loss of systemic circulation',  # Text content
            va='bottom',  # Vertical alignment
            ha='center',  # Horizontal alignment
            color='red')  # Text color
    
    
    # Manually setting x-axis labels for descent time to reflect absolute start time in minutes
    minutes = np.arange(start_sec / 60, end_sec / 60, (end_sec - start_sec) / (60 * 10))  # 10 labels
    ax[1].set_xticks(np.linspace(0, descent_time[patient_idx][1] - descent_time[patient_idx][0], len(minutes)))
    ax[1].set_xticklabels([f"{minute:.2f}" for minute in minutes])
    ax[1].set_xlabel("\n Time (minutes) from recording start \n")
    ax[1].set_ylabel("BBFv \n (cm/s)")
    ax2.set_ylabel("ABP \n (mmHg)")
    
    # plot CCA_Time
    ax[2].plot(V[0,CCA_Time - buffer:CCA_Time + buffer],color = 'blue',alpha = 0.5); 
    ax2 = ax[2].twinx()
    #ax2.set_ylim(0,100)
    ax2.plot(ABP[0,CCA_Time - buffer:CCA_Time + buffer],color = 'red',alpha = 0.5); 
    ax[2].set_ylabel("BBFv \n (cm/s)")
    ax2.set_ylabel("ABP \n (mmHg)")
    
    # change x axis to seconds
    ax[2].set_xticks(ticks = [0, 1000, 2000, 3000, 4000, 5000, 6000], \
                     labels = [-30, -20, -10, 0, 10, 20, 30])
    ax[2].set_xlabel("Seconds from loss of brain blood flow")
    
    # plot CA Time
    ax[3].plot(V[0,CA_Time - buffer:CA_Time + buffer],color = 'blue',alpha = 0.5); 
    ax2 = ax[3].twinx()
    #ax2.set_ylim(0,100)
    ax2.plot(ABP[0,CA_Time - buffer:CA_Time + buffer],color = 'red',alpha = 0.5); 
    ax[3].set_ylabel("BBFv \n (cm/s)")
    ax2.set_ylabel("ABP \n (mmHg)")
    # change x axis to seconds
    ax[3].set_xticks(ticks = [0, 1000, 2000, 3000, 4000, 5000, 6000], \
                     labels = [-30, -20, -10, 0, 10, 20, 30])
    ax[3].set_xlabel("Seconds from loss of systemic circulation")
    
    # plot line denoting CCA and CA
    ax[2].axvline(buffer, linestyle = "--", color = 'grey', alpha = 0.5)
    ax[3].axvline(buffer, linestyle = "--", color = 'black', alpha = 0.5)
    
    # optionally set the ylimit
    #for ax_idx in range(len(ax)):
    #    ax[ax_idx].set_ylim(0,100)
    
    plt.tight_layout()
    plt.show()
    FIGURES.append(fig)


# save figures as pdf
from matplotlib.backends.backend_pdf import PdfPages
import datetime

# Create a PDF file to save the figures
with PdfPages('PatientResults.pdf') as pdf:
    for fig in FIGURES:
        pdf.savefig(fig)  # Saves the current figure into a page of the PDF
        plt.close(fig)  # Close the figure to free memory
    
    
    # Optionally, you can add metadata to your PDF file
    d = pdf.infodict()
    d['Title'] = 'Patient Results - Depart'
    d['Author'] = 'Matthew Kolisnyk'
    d['Subject'] = 'Depart Patient-Specific Analyis'
    d['CreationDate'] = datetime.datetime.today()
    d['ModDate'] = datetime.datetime.today()





