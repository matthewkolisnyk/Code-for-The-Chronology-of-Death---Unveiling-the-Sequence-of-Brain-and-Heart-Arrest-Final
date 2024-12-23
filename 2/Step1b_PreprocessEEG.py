# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:33:42 2024

@author: mkolisnyk
"""

# try the same with EEG
# import library
from scipy.io import loadmat,savemat
import neurokit2 as nk
import os
import numpy as np
import mne
import matplotlib.pyplot as plt

# participant file name
file_name = ""
Raw = mne.io.read_raw_edf(file_name,  preload=True, exclude = ['Event','LOC',
'ROC',
'CHIN1',
'CHIN2',
'LAT1',
'LAT2',
'RAT1',
'RAT2',
'CHEST',
'ABD',
'FLOW',
'SNORE',
'DIF5',
'DIF6',
'POS',
'DC2',
'DC3',
'DC4',
'DC5',
'DC6',
'DC7',
'DC8',
'DC9',
'DC10',
'OSAT',
'PR'])

# creat montage file
ch_pos = {

'F3': np.array([-0.0502438,  0.0531112,  0.042192 ]),
'F4': np.array([0.0518362, 0.0543048, 0.040814 ]),
'F7': np.array([-0.0702629,  0.0424743, -0.01142  ]),
'F8': np.array([ 0.0730431,  0.0444217, -0.012    ]),
'Fz': np.array([0.0003122, 0.058512 , 0.066462 ]),
'P3': np.array([-0.0530073, -0.0787878,  0.05594  ]),
'P4': np.array([ 0.0556667, -0.0785602,  0.056561 ]),
'Pz': np.array([ 0.0003247, -0.081115 ,  0.082615 ]),
'C3': np.array([-0.0653581, -0.0116317,  0.064358 ]),
'C4': np.array([ 0.0671179, -0.0109003,  0.06358  ]),
'O1': np.array([-0.0294134, -0.112449 ,  0.008839 ]),
'O2': np.array([ 0.0298426, -0.112156 ,  0.0088   ]),
'Cz': np.array([ 0.0004009, -0.009167 ,  0.100244 ]),
'Fp1': np.array([-0.0294367,  0.0839171, -0.00699  ]),
'Fpz': np.array([ 0.0001123,  0.088247 , -0.001713 ]),
'Fp2': np.array([ 0.0298723,  0.0848959, -0.00708  ]),
'T3': np.array([-0.0841611, -0.0160187, -0.009346 ]),
'T5': np.array([-0.0724343, -0.0734527, -0.002487 ]),
'T4': np.array([ 0.0850799, -0.0150203, -0.00949  ]),
'T6': np.array([ 0.0730557, -0.0730683, -0.00254  ]),
'A1': np.array([-0.0860761, -0.0249897, -0.067986 ]),
'A2': np.array([ 0.0857939, -0.0250093, -0.068031 ])
}

coord_frame = 'mri'
nasion  = np.array([ 0.00146763,  0.08506715, -0.03483611], dtype='>f4')
lpa =  np.array([-0.08061612, -0.02908875, -0.04131077], dtype='>f4')
rpa = np.array([ 0.08436285, -0.02850276, -0.04127743], dtype='>f4')
hsp = None
hpi = None


dig = mne.channels.make_dig_montage(ch_pos=ch_pos, nasion=nasion, lpa=lpa, \
                              rpa=rpa, hsp=hsp, hpi=hpi, coord_frame=coord_frame)
    
Raw.set_channel_types({"F3":"eeg","F4":"eeg","F7":"eeg","F8":"eeg","Fz":"eeg","P3":"eeg","P4":"eeg",\
                       "Pz":"eeg","C3":"eeg","C4":"eeg","O1":"eeg","O2":"eeg","Cz":"eeg","Fp1":"eeg","Fpz":"eeg","Fp2":"eeg","T3":"eeg",\
                           "T5":"eeg","T4":"eeg","T6":"eeg","A1":"eeg","A2":"eeg","ECGL":"ecg","ECGR":"ecg"})
Raw.set_montage(dig)    

plt.plot(Raw.get_data()[5,:])
plt.show()

#eeg_raw.plot(proj = False, n_channels = len(eeg_raw.ch_names), remove_dc = False)

spectrum = Raw.compute_psd()
spectrum.plot(average=True, exclude="bads", picks = ['all'], amplitude=True)

# filter the data
eeg_raw = Raw.filter(l_freq = 0.1, h_freq = None, picks = ['eeg'])

# filter the ecg
eeg_raw = Raw.filter(l_freq = 0.5, h_freq = None, picks = ['ecg'])

eeg_raw.compute_psd(picks = ['all']).plot(picks = ['all'],amplitude = True)

eeg_raw = eeg_raw.filter(l_freq = None, h_freq = 55, picks = ['all'])

eeg_raw.compute_psd(picks = ['all']).plot(picks = ['all'],amplitude = True)

eeg_raw.plot(proj = False, n_channels = len(eeg_raw.ch_names), remove_dc = False)
plt.show()

# next notch filter
eeg_raw = eeg_raw.notch_filter(freqs = np.arange(60, 121, 60),picks = ['all'])

eeg_raw.compute_psd(picks = ['all']).plot(picks = ['all'],amplitude = True)
plt.show()

plt.plot(eeg_raw.get_data()[4,:])
plt.show()


# create bipolar montage
eeg_bipolar = eeg_raw.copy()
b_montage = ['F3 - P3','F4 - P4','Fz - Pz','F7 - T5','F8 - T6','A1 - A2']
anode = ['F3','F4','Fz','F7','F8','A1']; cathode = ['P3',"P4","Pz","T5","T6",'A2']
eeg_bipolar = mne.set_bipolar_reference(eeg_bipolar,anode = anode, cathode = cathode)

eeg_bipolar.plot()

plt.plot(eeg_bipolar.get_data()[4,:])
plt.show()

# save as mat file for later preprocessing .
savemat(file_name = "",mdict = {"data":eeg_bipolar.get_data()})


## now run the same preprocesing for the ECG file that came from the physiological recorder
import pandas as pd
ECG = pd.read_csv("", sep = "\t")['ECG1'].values

# create RAW object for ECG
ecg_info = mne.create_info(ch_names = ['ECG'],sfreq = 300, ch_types = ['ecg'],verbose = True)
ecg_raw = mne.io.RawArray(ECG.reshape(1,-1),info = ecg_info,first_samp = 0, copy = 'auto', verbose = True)

# plot raw data
#ecg_raw.plot(duration = 360, proj = False, n_channels = len(ecg_raw.ch_names), remove_dc = False)

spectrum = ecg_raw.compute_psd(picks=['ecg'])
spectrum.plot(average=True, picks=['ecg'], exclude="bads", amplitude=True)

# find the high amplitude components
np.argsort(-spectrum.get_data(picks = 'ecg'))[0][0:100]

# filter the data
ecg_raw = ecg_raw.filter(l_freq = 0.5, h_freq = 55, picks = ['ecg'])

ecg_raw.compute_psd(picks = ['ecg']).plot(picks = ['ecg'],amplitude = True)

# notch filter
ecg_raw = ecg_raw.notch_filter(freqs = np.arange(60, 121, 60),picks = ['ecg'])

#ecg_raw = ecg_raw.notch_filter(freqs = None,picks = ['ecg'],method='spectrum_fit')


ecg_raw.compute_psd(picks = ['ecg']).plot(picks = ['ecg'],amplitude = True)
plt.show()

# resample ecg
ecg_raw = ecg_raw.resample(sfreq = 256)

# save as mat file for later analysis.
savemat(file_name = "",mdict = {"data":ecg_raw.get_data()})