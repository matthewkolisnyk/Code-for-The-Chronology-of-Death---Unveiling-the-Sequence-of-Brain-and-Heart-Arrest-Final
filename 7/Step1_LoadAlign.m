%% Script to load and align data for patient 7
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


%% now that we are resampled, trim time series 
% note we only have second level resolution with phys data

TCD_start = "";
Phys_start = "";
start_diff_sec = TCD_start - Phys_start;

% convert to sample 
start_diff_sample = tfs * start_diff_sec;

% note that TCD started before ABP

% set 0s to nans at the beginning of the phys recording
start_buffer = 0;
BPt = BP;
noise_buffer = 200;
BPt(1:noise_buffer) = nan;
ECGt(1:noise_buffer) = nan;
V = V(abs(start_diff_sample):end);

%% align time series
% align the arterial blood pressure and TCD signal using good quality data
% prior to cessation (where ART and TCD begin to look very different). 
lag_window = 100000; clean_window = 730000; skip_start = 200;
% use cross correlation to measure whether this alignment works 
[r,lags] = xcorr(V(skip_start:clean_window),BPt(skip_start:clean_window),lag_window);
[~,ccidx] = max(r);

% shift results. Very small (adjustment from max for slightly better alignment 
shift = lag_window - ccidx + 18; 

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


%% save results
Vtime = 1:length(Vs);
TCDClean.V = Vs'; TCDClean.Time = Vtime; TCDClean.BP = BPt'; TCDClean.ECG = ECG';
TCDClean.PP = Last_PP_sample; TCDClean.CriticalPoints = [Last_PP_sample,Last_QRS_sample];
TCDClean.QRS = Last_QRS_sample; TCDClean.ISO = nan;
save("TCDClean.mat","TCDClean")

