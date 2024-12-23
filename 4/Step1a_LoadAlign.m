%% Script to load and align data for patient 4
% file names
TCD_filename = "";
Phys_filename = "";

% load TCD 
TCD = readtable(TCD_filename);
tfs = 100;  % sampling frequency for TCD

% load Phys
Phys = readtable(Phys_filename);
bfs = 300; % sampling frequency for physiological data

% extract exact measures of interest (velocity, arterial pulse pressure,
% and ecg
V = TCD{:,"Gate1ToProbeEnv_"}; BP = Phys{:,"P1"}; ECG = Phys{:,"ECG1"};
[nTCD,~] = size(TCD); [nBP,~] = size(BP);

% downsample BP data to TCD sampling frequency
BP = resample(BP,tfs,bfs);

% downsample ECG data to TCD sampling frequency
ECG = resample(ECG,tfs,bfs);

%% now that we are resampled, trim time series 
% note we only have second level resolution with phys data

% what are the start times for each of the devices
TCD_start = "";
Phys_start = "";


% what is the difference in start times
start_diff_sec = round(TCD_start) - Phys_start;

% convert to sample space
start_diff_sample = tfs * start_diff_sec;

% set 0s to nans at the beginning of the phys recording
start_buffer = 0; noise_buffer = 1000;
BPt = BP(start_diff_sample-start_buffer:end);

% BP noise index
BP_noise = find(BPt(1:noise_buffer) == 0);
BPt(BP_noise) = nan; 

% same for ECG
ECGt = ECG(start_diff_sample-start_buffer:end);

%% align time series
% align the arterial blood pressure and TCD signal using good quality data
% prior to cessation (where ART and TCD begin to look very different). 

lag_window = 100000; clean_window = 300000; skip_start = 1;
% use cross correlation to measure whether this alignment works 
[r,lags] = xcorr(V(skip_start:clean_window),BPt(skip_start:clean_window),lag_window);
stem(lags,r);
[~,ccidx] = max(r);

% shift results
shift = lag_window - ccidx; % note negative values mean BP lags V, 
% whereas positive values means V lags BP

% slightly adjusting alignment to match correspondence with PVCs
shift = shift - 116;

% Shift the V signal
if shift > 0
    Vs = [NaN(shift, 1); V(1:end-shift)];
elseif shift < 0
    Vs = [V(abs(shift)+1:end); NaN(abs(shift), 1)];
else
    Vs = V;
end

%% store relevant timing information (see Table 2 in the main manuscript)
Last_PP_sample = []; Last_QRS_sample = []; 


%% save adjusted signals
Vtime = 1:length(Vs);
TCDClean.V = Vs'; TCDClean.Time = Vtime; TCDClean.BP = BPt';
TCDClean.ECG = ECGt';
TCDClean.PP = Last_PP_sample; TCDClean.CriticalPoints = [Last_PP_sample,Last_QRS_sample];
TCDClean.QRS = Last_QRS_sample; TCDClean.ISO = nan;
save("TCDClean.mat","TCDClean")