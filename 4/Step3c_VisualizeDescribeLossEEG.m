%% investigate EEG signals 

% load data
load("TCDClean.mat")
EEG = TCDClean.EEG(:,1:5);

ECG = TCDClean.ECG;
EEGECG = TCDClean.EEGECG;
% drop after QRS
buffer = 20000;

EEG = EEG(1:TCDClean.QRS+buffer,:);

% confirm QRS timing after shift
figure
plot(EEGECG(TCDClean.QRS-buffer:TCDClean.QRS+buffer)*10e4)
hold on
plot(ECG(TCDClean.QRS-buffer:TCDClean.QRS+buffer))


% try looking how the EEG timeseries changes over time with a rolling
% average
tfs = 100;
EEGmov = movmean(EEG,tfs,1,'omit');

% or just a regular average across channels
EEGmean = median(EEG,2);

% plot results

% also look after larger artifact
clean_window = 372440;

figure
plot(EEGmean)

% find the scale...
figure
plot(EEGmean(clean_window:end) * 10e1)

figure
plot(EEGmean(clean_window:end) * 10e2)

figure
plot(EEGmean(clean_window:end) * 10e3)

figure
plot(EEGmean(clean_window:end) * 10e4)

figure
plot(EEGmean(clean_window:end) * 10e5)


% collapse over fs

% Parameters
windowLengthSec = 60; % Window length in seconds
windowLengthSamples = windowLengthSec * tfs; % Convert window length to samples

% Number of complete windows in the signal
numWindows = floor(length(EEGmean) / windowLengthSamples);

% Initialize array to store average values
windowAverages = zeros(1, numWindows);
windowStandardDeviation = zeros(1, numWindows);

% Calculate average for each window
for k = 1:numWindows
    startIndex = (k-1) * windowLengthSamples + 1;
    endIndex = k * windowLengthSamples;
    windowAverages(k) = mean(EEGmean(startIndex:endIndex)* 10e5);
    windowStandardDeviation(k) = std(EEGmean(startIndex:endIndex)* 10e5);
end


% Plot the original moving average and the window averages
figure;
hold on;
plot(windowAverages, 'ro-');
plot(windowStandardDeviation, 'ko-');
xlabel('Time (s)');
ylabel('Amplitude');
title('Window Average Amplitude and Deviation');

% isolate iso eeg time
IsoWindow = 95;

isoStartIndex = (IsoWindow-1) * windowLengthSamples + 1;
isoEndIndex = IsoWindow * windowLengthSamples;


%% plot old ISO value and new ISO value
figure
plot(EEGmean(isoStartIndex:isoEndIndex))
hold on
plot(EEGmean(TCDClean.ISO:TCDClean.ISO+windowLengthSamples))
hold off

% store your isoeeg value
TCDClean.ISO = isoStartIndex;

save("TCDClean.mat","TCDClean")


