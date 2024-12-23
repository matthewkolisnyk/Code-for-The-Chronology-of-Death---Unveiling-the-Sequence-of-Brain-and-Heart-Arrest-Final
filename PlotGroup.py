# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:47:57 2024

@author: mkolisnyk
"""

# plot TCD, ECG, and Phys results
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# go to current directory
os.chdir(r"D:\MATTFILES\TCD_Clean")

# load in excel with timing information 
file = 'Results-Table-Clean-v2.csv'

# load in patient file
data = pd.read_csv(file)

# store figures to save to a pdf
FIGURES = []

# load AMA-typical plotting format

# Set the global font size, figure size, and line width for plots
plt.rcParams.update({
    'font.size': 8,  # Adjust based on your journal's requirements
    'figure.figsize': (4, 2.2),  # Typical figure size, but adjust as needed
    'figure.dpi': 400,
    'lines.linewidth': 2,
    'axes.labelsize': 8,  # Adjust label size as necessary
    'font.family': 'Arial',
    'axes.titlesize': 8,
    'xtick.labelsize': 8,
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


# Use a simple color palette
sns.set_palette('Grays')  # For Seaborn plots, or choose a palette appropriate for color-blind readers if using color


# drop some of the ISO stuff
#data = data.drop(['ISO - CA', 'ISO - CCA',
#'CA Time (Samples).1', 'CA SP.1', 'CA DP.1', 'CA PP.1', 'CA HR.1',
#'CA MAP.1', 'CCA Time (Samples).1', 'CCA SP.1', 'CCA DP.1', 'CCA PP.1',
#'CCA HR.1', 'CCA MAP.1', 'ISO Time','Descent - CCA','Descent - CA'],axis = 1)

# gather some simple information about the data
nPat = data.shape[0]
cmap = sns.color_palette("PuBuGn", n_colors=nPat)

import matplotlib.patches as mpatches

# produce a simple boxplot of timing results
fig2, ax = plt.subplots(1, 1)
sns.boxplot(data=data, x="CCA Time (Seconds)",width = .4)
sns.swarmplot(data=data, x='CCA Time (Seconds)', size=5, color='gray')
plt.xlabel("Time from Loss of Systemic Circulation \n(Seconds)")
plt.axvline(0, linestyle='--', color='red',alpha = 0.5)  # Reference line at zero
plt.axvline(300, linestyle='--', color='black',alpha = 0.15)

# Shading regions based on conditions
ax.axvspan(-500, 0, color='lightblue', alpha=0.3)  # Time of loss of brain blood flow > time of loss of systemic circulation
ax.axvspan(0, 500, color='lightcoral', alpha=0.3)  # Time of loss of brain blood flow < time of loss of systemic circulation

# Create legend patches
#blue_patch = mpatches.Patch(color='lightblue', alpha=0.3, label='Loss of brain blood flow > Loss of systemic circulation')
#red_patch = mpatches.Patch(color='lightcoral', alpha=0.3, label='Loss of systemic circulation > Loss of brain blood flow ')

ax.text(-400,0.40,'Loss of Brain Blood Flow before \n Loss of Systemic Circulation', {'family': 'Arial','size':6,'color':'midnightblue'})
ax.text(100,0.40,'Loss of Brain Blood Flow after \n Loss of Systemic Circulation', {'family': 'Arial','size':6,'color':'darkred'})
ax.text(200,-0.58,'Current Guidelines for \n Circulatory Death', {'family': 'Arial','size':7})
ax.text(-100,-0.58,'Loss of Systemic \n     Circulation', {'family': 'Arial','size':7})
# Add legend to the plot
#plt.legend(handles=[blue_patch, red_patch], loc=[0,1], frameon=True)

# Customize remaining plot features
ax.get_yaxis().set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim([-500, 500])
plt.tight_layout()
plt.show()


# find ways to plot ISO timing

# create iso dataframe
dfISO = pd.concat([data.dropna()['ISO - CCA'],data.dropna()['ISO - CA'],data.dropna()['ID']],axis = 1)


# Reshape data to a long format
dfISO_long = pd.melt(dfISO, id_vars=['ID'], value_vars=['ISO - CCA', 'ISO - CA'],
                     var_name='Event', value_name='Value')

# adjust value to swithc ISO-CCA to CCA-ISO
dfISO_long.loc[dfISO_long["Event"] == 'ISO - CCA',"Value"] = dfISO_long.loc[dfISO_long["Event"] == 'ISO - CCA']['Value'] * -1

gray_palette = sns.color_palette("Greys", n_colors=len(dfISO_long['ID'].unique()) + 1)[1:]  # Exclude the lightest color

fig2, ax = plt.subplots(1, 1)

# Plot using stripplot, setting 'Event' as the y-axis to separate the two conditions
sns.swarmplot(data=dfISO_long, x='Value', y='Event', size=5, hue="ID", palette=gray_palette, alpha=0.8, legend=True)

# Add labels and customizations
plt.xlabel("Seconds")
plt.ylabel("")
plt.axvline(0, linestyle='--', color='red',alpha = 0.5)  # Reference line at zero
plt.axvline(300, linestyle='--', color='black',alpha = 0.5)  # Reference line at zero

# Shading regions based on conditions
ax.axvspan(-500, 0, color='lightblue', alpha=0.3)  # Time of loss of brain blood flow > time of loss of systemic circulation
ax.axvspan(0, 500, color='lightcoral', alpha=0.3)  # Time of loss of brain blood flow < time of loss of systemic circulation

# Customize remaining plot features
ax.spines['left'].set_visible(True)
ax.set_xlim([-500, 500])

ax.set_yticklabels(["Loss of \n Brain Blood Flow", "Loss of \n Systemic Circulation"])

plt.tight_layout()
plt.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position if necessary
plt.show()

FIGURES.append(fig2)

fig3, ax = plt.subplots(1, 1)

# Plot using stripplot, setting 'Event' as the y-axis to separate the two conditions
sns.swarmplot(data=dfISO_long[dfISO_long['Event'] == 'ISO - CCA'], x='Value', size=5, palette=gray_palette, alpha=0.8, legend=True)

# Add labels and customizations
plt.xlabel("Time from Loss of Brain Activity \n (Seconds)")
plt.ylabel("")
plt.yticks()
plt.axvline(0, linestyle='--', color='red',alpha = 0.5)  # Reference line at zero
#plt.axvline(300, linestyle='--', color='black',alpha = 0.5)  # Reference line at zero

# Shading regions based on conditions
ax.axvspan(-500, 0, color='lightblue', alpha=0.3)  # Time of loss of brain blood flow > time of loss of systemic circulation
ax.axvspan(0, 500, color='lightcoral', alpha=0.3)  # Time of loss of brain blood flow < time of loss of systemic circulation


ax.text(-460,0.40,'Loss of Brain Blood Flow before \n Loss of Brain Activity', {'family': 'Arial','size':6})
ax.text(40,0.40,'Loss of Brain Blood Flow after \n Loss of Brain Activity', {'family': 'Arial','size':6})
#ax.text(200,-0.58,'Current Guidelines for \n Death Determination', {'family': 'Arial','size':6})
ax.text(-120,-0.58,'Loss of Brain Activity', {'family': 'Arial','size':6})

# Customize remaining plot features
ax.spines['left'].set_visible(True)
ax.set_xlim([-500, 500])
plt.tick_params(left = False) 

ax.set_yticklabels([], minor=True)

plt.tight_layout()
#plt.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position if necessary
plt.show()


fig4, ax = plt.subplots(1, 1)

# Plot using stripplot, setting 'Event' as the y-axis to separate the two conditions
sns.swarmplot(data=dfISO_long[dfISO_long['Event'] == 'ISO - CA'], x='Value', size=5, palette=gray_palette, alpha=0.8, legend=True)

# Add labels and customizations
plt.xlabel("Time from Loss of Systemic Circulation \n (Seconds)")
plt.ylabel("")
plt.yticks()
plt.axvline(0, linestyle='--', color='red',alpha = 0.5)  # Reference line at zero
plt.axvline(300, linestyle='--', color='black',alpha = 0.5)  # Reference line at zero

# Shading regions based on conditions
ax.axvspan(-500, 0, color='lightblue', alpha=0.3)  # Time of loss of brain blood flow > time of loss of systemic circulation
ax.axvspan(0, 500, color='lightcoral', alpha=0.3)  # Time of loss of brain blood flow < time of loss of systemic circulation

ax.text(-460,0.40,'Loss of Brain Activity before \n Loss of Systemic Circulation', {'family': 'Arial','size':6})
ax.text(40,0.40,'Loss of Brain Activity after \n Loss of Systemic Circulation', {'family': 'Arial','size':6})
ax.text(200,-0.58,'Current Guidelines for \n Death Declaration', {'family': 'Arial','size':6})
ax.text(-180,-0.58,'Loss of Systemic Circulation', {'family': 'Arial','size':6})

# Customize remaining plot features
ax.spines['left'].set_visible(True)
ax.set_xlim([-500, 500])
plt.tick_params(left = False) 

ax.set_yticklabels([], minor=True)

plt.tight_layout()
#plt.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc='upper left')  # Adjust legend position if necessary
plt.show()


# drop some of the ISO stuff
data = data.drop(['ISO - CA', 'ISO - CCA',
'CA Time (Samples).1', 'CA SP.1', 'CA DP.1', 'CA PP.1', 'CA HR.1',
'CA MAP.1', 'CCA Time (Samples).1', 'CCA SP.1', 'CCA DP.1', 'CCA PP.1',
'CCA HR.1', 'CCA MAP.1', 'ISO Time','Descent - CCA','Descent - CA'],axis = 1)

# fold dataset 
CCA = data.filter(regex = "CCA")
CCA['Event'] = "CCA"
CCA['Patient'] = np.arange(1,nPat + 1)
CCA['Patient'] = CCA['Patient'].apply(lambda x: "Patient " + str(x))
CCA = CCA.values
CA = data.filter(regex = "(?<!C)CA")
CA['Event'] = "CA"
CA['Patient'] = np.arange(1,nPat + 1)
CA['Patient'] = CA['Patient'].apply(lambda x: "Patient " + str(x))
CA = CA.values 

dataEvent = pd.DataFrame(np.concatenate([CCA,CA],axis = 0),columns = ["Seconds", "Samples","Systolic Pressure",\
                                                                      "Diastolic Pressure","Pulse Pressure",\
                                                                          "Heart Rate", "Mean Arterial Pressure",\
                                                                              "Event","Patient"])

# Use a simple color palette                                                                              
#sns.set_palette('colorblind')  # For Seaborn plots, or choose a palette appropriate for color-blind readers if using color
print(list(sns.palettes.SEABORN_PALETTES.keys()))
cmap = sns.color_palette("PuBuGn", n_colors=nPat)

# Set the global font size, figure size, and line width for plots
plt.rcParams.update({
    'font.size': 8,  # Adjust based on your journal's requirements
    'figure.figsize': (2, 1.8),  # Typical figure size, but adjust as needed
    'figure.dpi': 400,
    'lines.linewidth': 2,
    'axes.labelsize': 8,  # Adjust label size as necessary
    'font.family': 'Arial',
    'axes.titlesize': 8,
    'xtick.labelsize': 8,
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


import matplotlib.patches as mpatches

fig3, ax = plt.subplots(nrows=1, ncols=2)

# Plotting each subplot
sns.boxplot(data=dataEvent, y="Pulse Pressure", x="Event", ax=ax[0], palette="Greys",width = 0.3)
sns.pointplot(data=dataEvent, y='Pulse Pressure', x="Event", ax=ax[0], hue="Patient", markersize=2, alpha=0.30, color='k')
ax[0].get_xaxis().set_visible(False)  # Remove x-axis ticks and labels
ax[0].legend([], frameon=False)

sns.boxplot(data=dataEvent, y="Mean Arterial Pressure", x="Event", ax=ax[1], palette="Greys",width = 0.3)
sns.pointplot(data=dataEvent, y='Mean Arterial Pressure', x="Event", ax=ax[1], hue="Patient", markersize=2, alpha=0.30, color='k')
ax[1].get_xaxis().set_visible(False)  # Remove x-axis ticks and labels
ax[1].legend([], frameon=False)

# Custom legend for x-axis categories
#loss_brain_flow_patch = mpatches.Patch(color="grey", label="Loss of Systemic Circulation")
#loss_sys_circ_patch = mpatches.Patch(color="darkgrey", label="Loss of Brain Blood Flow")

# Add custom legend above both plots
#plt.legend(handles=[loss_brain_flow_patch, loss_sys_circ_patch], loc='upper center', bbox_to_anchor=(0, 1.25), frameon=True)

# Adjust subplot layout to prevent squishing
fig3.subplots_adjust(top=0.85)  # Add space at the top

plt.tight_layout()
plt.show()

FIGURES.append(fig3)

# save figures as pdf
from matplotlib.backends.backend_pdf import PdfPages
import datetime

# Create a PDF file to save the figures
with PdfPages('GroupResults.pdf') as pdf:
    for fig in FIGURES:
        pdf.savefig(fig)  # Saves the current figure into a page of the PDF
        plt.close(fig)  # Close the figure to free memory
    
    
    # Optionally, you can add metadata to your PDF file
    d = pdf.infodict()
    d['Title'] = 'Group Results - Depart'
    d['Author'] = 'Matthew Kolisnyk'
    d['Subject'] = 'Depart Group Analyis'
    d['CreationDate'] = datetime.datetime.today()
    d['ModDate'] = datetime.datetime.today()

