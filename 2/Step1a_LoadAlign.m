%% Script to load and align data for patient 1
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
V = TCD{:,"Gate2ToProbeEnv_"}; BP = Phys{:,"P1"}; ECG = Phys{:,"ECG1"};
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

% set time before TCD start to nan
start_buffer = 0;
BPt = BP(start_diff_sample-start_buffer:end);

% same with ECG 
ECGt = ECG(start_diff_sample-start_buffer:end);

%% align time series
% align the arterial blood pressure and TCD signal using good quality data
% prior to cessation (where ART and TCD begin to look very different). 

% how many samples do we want to compute correlation over?
lag_window = 800000; 

% how long into the do we have good quality signal
clean_window = 800000; 

% where does the good quality signal start? 
skip_start = 1;

% compute cross correlation for TCD and ART over these times
[r,lags] = xcorr(V(skip_start:clean_window),BPt(skip_start:clean_window),lag_window);

% grab the maximum cross correlation
[~,ccidx] = max(r);

% shift the velocity time series by the maximum cross correlation
shift = lag_window - ccidx; 

% two maxima in the cross correlation plot. Visualize inspection (confirmed by MS) prefers the second choice. 
% large shift is very likely due to corrupt TCD starting time in the file. 
shift = -67430;

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

