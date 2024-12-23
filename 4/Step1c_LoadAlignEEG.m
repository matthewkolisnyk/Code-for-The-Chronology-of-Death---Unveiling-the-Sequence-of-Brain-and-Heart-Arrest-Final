%% validate EEG 
clear
close all
clc

EEG = load("").data';
% load in channel info 
fs = 200;
ch_names = {'Oz',
 'Fp1',
 'T3',
 'O1',
 'C3',
 'Fpz',
 'Cz',
 'Fp2',
 'T4',
 'O2',
 'C4',
 'X1',
 'X2',
 'F3-P3',
 'F4-P4',
 'Fz-Pz',
 'F7-T5',
 'F8-T6',
 'A1-A2'};


% load in ECG
ECG = load("").data';

% assume recording time is the same
EEG_start = "";
EEG_end = "";
Phys_start = "";
Phys_end = "";
start_diff_sec = "";
end_diff_sec = "";

% convert to sample 
start_diff_sample = fs * start_diff_sec;
end_diff_sample = fs * end_diff_sec;

% same with ECG 
ECGt = ECG(start_diff_sample:end);

% get size info
[nsEEG, nChan] = size(EEG);
[nsECG,~] = size(ECGt);


%% align by final QRS

% final QRS for both EEG and ECG 
Last_QRS_EEG = 2221720; 
Last_QRS_ECG = 2211270;

shift = Last_QRS_EEG - Last_QRS_ECG;

if shift < 0
    EEGt = [NaN(abs(shift), length(ch_names)); EEG(:,:)];
elseif shift > 0
    EEGt = [EEG(abs(shift)+1:end,:); NaN(abs(shift), length(ch_names))];
else
    EEGt = EEG;
end

% plot to verify shift
figure
plot(ECGt)
hold on
plot(EEGt(:,12) * 10e3)
hold off



%% look at alignment relative to final QRS
Phys_start = "";
Last_QRS = "";

% relative to BP?
Last_QRS_sec = Last_QRS - Phys_start;
Last_QRS_sample = Last_QRS_sec * (fs);

Last_QRS_sample = Last_QRS_sample - 3989;

figure
plot(ECGt)
hold on
plot(EEGt(:,14)*10e4)
xline(Last_QRS_sample)
hold off

buffer = 5000;
figure
hold on
plot(ECGt(Last_QRS_sample - buffer:Last_QRS_sample + buffer)*10e1)
plot(EEGt(Last_QRS_sample - buffer:Last_QRS_sample + buffer,14)*10e4)
hold off

figure
hold on
plot(ECGt(Last_QRS_sample - buffer:Last_QRS_sample + buffer))
plot(EEGt(Last_QRS_sample - buffer:Last_QRS_sample + buffer,16:end)*10e4)
hold off


%% now look relative to last ABP
Last_PP = "";

Last_PP_sec = Last_PP - Phys_start;
Last_PP_sample = Last_PP_sec * (fs);

figure
plot(ECGt)
hold on
plot(EEGt(:,14))
xline(Last_PP_sample)
hold off

buffer = 5000;
figure
hold on
plot(ECGt(Last_PP_sample - buffer:Last_PP_sample + buffer))
plot(EEGt(Last_PP_sample - buffer:Last_PP_sample + buffer,13)*10e4)
hold off

buffer = 5000;
figure
hold on
plot(ECGt(Last_PP_sample - buffer:Last_PP_sample + buffer))
plot(EEGt(Last_PP_sample - buffer:Last_PP_sample + buffer,16:end)*10e4)
hold off


%% proposed time for iso eeg
Last_EEG = "";

Last_EEG_sec = Last_EEG - Phys_start;
Last_EEG_sample = Last_EEG_sec * (fs);

figure
plot(ECGt)
hold on
plot(EEGt(:,14))
xline(Last_EEG_sample)
hold off

buffer = 5000;
figure
hold on
plot(ECGt(Last_EEG_sample - buffer:Last_EEG_sample + buffer))
plot(EEGt(Last_EEG_sample - buffer:Last_EEG_sample + buffer,13)*10e4)
hold off

buffer = 5000;
figure
hold on
plot(ECGt(Last_EEG_sample - buffer:Last_EEG_sample + buffer))
plot(EEGt(Last_EEG_sample - buffer:Last_EEG_sample + buffer,16:end)*10e4)
hold off


%% store EEG with rest of TCD clean variables
% downsample to match rest
tfs = 100; bfs = fs;
EEGt = resample(EEGt,tfs,bfs);

TCD_start = "";
Phys_start = "";
start_diff_sec = TCD_start - EEG_start;

% convert to sample 
start_diff_sample = tfs * start_diff_sec;


% set 0s to nans at the beginning of the phys recording
start_buffer = 0;
noise_buffer = 1000;

EEGt = EEGt(start_diff_sample-start_buffer:end,:);

% load TCD
load("TCDClean.mat")
TCDClean.EEG = EEGt(:,end-5:end);
TCDClean.EEGECG = EEGt(:,12:13);
save("TCDClean.mat",'TCDClean')