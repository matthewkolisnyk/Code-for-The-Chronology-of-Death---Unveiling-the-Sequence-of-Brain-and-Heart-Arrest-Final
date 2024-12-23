%% get key information from TCD and ABP signals at the time of CCA
clear; close all;
% load in timing information and signals
load("TCDClean.mat")

% what is in the structure
TCDClean

% get the parts of the data you will need
BloodVelocity = TCDClean.V;
BloodPressure = TCDClean.BP;
ECG = TCDClean.ECG;
EEG = mean(TCDClean.EEG * 10e4,2);
CCATime = TCDClean.CCA_Matt;

% shift CCA Time to be at the beginning of the final cycle
CCATime = CCATime + 5;

% plot the time series about the last time pulse pressure was less than
% 5mmHg
buffer = 6000; % 60 seconds

% extract 60 seconds pre and post CA
BloodVelocityCCA = BloodVelocity(CCATime - buffer: CCATime + buffer);
BloodPressureCCA = BloodPressure(CCATime - buffer: CCATime + buffer);
EEGCCA = EEG(CCATime - buffer: CCATime + buffer);
ECGCCA = ECG(CCATime - buffer: CCATime + buffer);

figure(1)
hold on
plot(BloodVelocityCCA,'b');
plot(BloodPressureCCA,'r');
plot(ECGCCA,'g');
plot(EEGCCA,'k');
xline(buffer)
%xticklabels(-buffer/100:(1/(buffer*2)):buffer/100)
%xlabel("Seconds")
hold off

% lets find the peaks
fs = 100;
WR = fs;
[peakValues, peakIndices] = findpeaks(BloodPressureCCA,'MinPeakDistance',WR,'MinPeakHeight',mean(BloodPressureCCA));

% Find troughs by negating the signal
[troughValues, troughIndices] = findpeaks(-BloodPressureCCA,'MinPeakDistance',WR);
troughValues = -troughValues; % Negate back to get the actual trough values

% Visualization
figure(3);
plot(BloodPressureCCA);
hold on;
plot(peakIndices, peakValues, 'ro'); % Peaks as red circles
plot(troughIndices, troughValues, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off


% do the same for ECG

[peakValuesECG, peakIndicesECG] = findpeaks(ECGCCA,'MinPeakDistance',WR,'MinPeakHeight',mean(ECGCCA));

% Find troughs by negating the signal
[troughValuesECG, troughIndicesECG] = findpeaks(-ECGCCA,'MinPeakDistance',WR);
troughValuesECG = -troughValuesECG; % Negate back to get the actual trough values

% Visualization
figure(4);
plot(ECGCCA);
hold on;
plot(peakIndicesECG, peakValuesECG, 'ro'); % Peaks as red circles
plot(troughIndicesECG, troughValuesECG, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off


% calculate systolic pulse pressure
% find index of blood pressure that precedes the CCA time (as defined by
% TCD)
peakIdx = find(sign(peakIndices - buffer) == -1);
troughIdx = find(sign(troughIndices - buffer) == -1);

peakIdxECG = find(sign(peakIndicesECG - buffer) == -1);
troughIdxECG = find(sign(troughIndicesECG - buffer) == -1);

% extract final cycle
last_cycle_time = peakIndices(peakIdx(end -1)):peakIndices(peakIdx(end) + 1);

% compute derivative of BP over this duration
last_cycle_dt = diff(BloodPressureCCA(last_cycle_time));
last_cycle_dtdt = diff(last_cycle_dt);

% plot final cycle
figure(8)
plot(BloodPressureCCA(last_cycle_time))
hold on
plot(last_cycle_dt)
plot(last_cycle_dtdt)
hold off

% trim the other cycles off the final cycle
% how many steps back from the peak do we need to go before having a
% negative derivative
last_peak_idx = find(BloodPressureCCA(last_cycle_time) == peakValues(peakIdx(end)));
start_trough = find(sign(last_cycle_dt(1:last_peak_idx-1)) == -1);
start_trough = start_trough(end) + 1;

% same for end
last_trough_idx = find(BloodPressureCCA(last_cycle_time) == peakValues(peakIdx(end)+1));
last_trough = find(sign(last_cycle_dt(1:last_trough_idx-1)) == -1);
last_trough = last_trough(end) + 1;

% now to calculate values over the final cycle
last_cycle = BloodPressureCCA(last_cycle_time(start_trough:last_trough));
plot(last_cycle);

SystolicPressure = max(last_cycle);
DiastolicPressure = last_cycle(end);
PulsePressure = SystolicPressure - DiastolicPressure;
MAP = trapz(last_cycle);
MAP = MAP/length(last_cycle);
r_r =  peakIndices(peakIdx(end)) - peakIndices(peakIdx(end -1));
% find r_r peak using ecg
last_peak_ecg = find(sign(peakIndicesECG - peakIndices(peakIdx(end)))==-1);
last_peak_ecg_time = peakIndicesECG(last_peak_ecg(end)+1) - peakIndicesECG(last_peak_ecg(end));
heart_rate = 60/(r_r/fs);
heart_rate_ecg = 60/(last_peak_ecg_time/fs);


% display results
fprintf("Systolic Pressure is: %f \n",SystolicPressure)
fprintf("Diastolic Pressure is: %f \n",DiastolicPressure)
fprintf("Pulse Pressure is: %f \n",PulsePressure)
fprintf("Heart Rate is: %f \n",heart_rate)
fprintf("Heart Rate ECG is: %f \n",heart_rate_ecg)
fprintf("MAP (Integral) is: %f \n",MAP)



%% now for CA
clear; close all;

% load in timing information and signals
load("TCDClean.mat")

% what is in the structure
TCDClean

% get the parts of the data you will need
BloodVelocity = TCDClean.V;
BloodPressure = TCDClean.BP;
ECG = TCDClean.ECG;
CATime = TCDClean.PP;

% bump CATime to the beginning of the final cycle
CATime = CATime - 39;


% plot the time series about the last time pulse pressure was less than
% 5mmHg
buffer = 6000; % 60 seconds

% extract 60 seconds pre and post CA
BloodVelocityCA = BloodVelocity(CATime - buffer: CATime + buffer);
BloodPressureCA = BloodPressure(CATime - buffer: CATime + buffer);
ECGCA = ECG(CATime - buffer: CATime + buffer);

figure(1)
hold on
plot(BloodVelocityCA,'b');
plot(BloodPressureCA,'r');
plot(ECGCA,'g');
xline(buffer)
%xticklabels(-buffer/100:(1/(buffer*2)):buffer/100)
%xlabel("Seconds")
hold off

% lets find the peaks
fs = 100;
WR = fs/1.25;
[peakValues, peakIndices] = findpeaks(BloodPressureCA,'MinPeakDistance',WR,'MinPeakHeight',mean(BloodPressureCA(1:buffer)));

% Find troughs by negating the signal
[troughValues, troughIndices] = findpeaks(-BloodPressureCA,'MinPeakDistance',WR);
troughValues = -troughValues; % Negate back to get the actual trough values

% Visualization
figure(3);
plot(BloodPressureCA);
hold on;
plot(peakIndices, peakValues, 'ro'); % Peaks as red circles
plot(troughIndices, troughValues, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off

% do the same for ECG

[peakValuesECG, peakIndicesECG] = findpeaks(ECGCA,'MinPeakDistance',WR,'MinPeakHeight',mean(ECGCA));

% Find troughs by negating the signal
[troughValuesECG, troughIndicesECG] = findpeaks(-ECGCA,'MinPeakDistance',WR);
troughValuesECG = -troughValuesECG; % Negate back to get the actual trough values

% Visualization
figure(4);
plot(ECGCA);
hold on;
plot(peakIndicesECG, peakValuesECG, 'ro'); % Peaks as red circles
plot(troughIndicesECG, troughValuesECG, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off

% calculate systolic pulse pressure
% find index of blood pressure that precedes the CCA time (as defined by
% TCD)
peakIdx = find(sign(peakIndices - buffer) == 1);
troughIdx = find(sign(troughIndices - buffer) == 1);

% extract final cycle
last_cycle_time = peakIndices(peakIdx(1)-1):peakIndices(peakIdx(1) + 1);

% compute derivative of BP over this duration
last_cycle_dt = diff(BloodPressureCA(last_cycle_time));
last_cycle_dtdt = diff(last_cycle_dt);

% plot final cycle
figure(8)
plot(BloodPressureCA(last_cycle_time))
hold on
plot(last_cycle_dt)
plot(last_cycle_dtdt)
hold off

% trim the other cycles off the final cycle
% how many steps back from the peak do we need to go before having a
% negative derivative
last_peak_idx = find(BloodPressureCA(last_cycle_time) == peakValues(peakIdx(1)));
start_trough = find(sign(last_cycle_dt(1:last_peak_idx-1)) == -1);
start_trough = start_trough(end) + 1;

% same for end
last_trough_idx = find(BloodPressureCA(last_cycle_time) == peakValues(peakIdx(1)+1));
last_trough = find(sign(last_cycle_dt(1:last_trough_idx-1)) == -1);
last_trough = last_trough(end) + 1;

% now to calculate values over the final cycle
last_cycle = BloodPressureCA(last_cycle_time(start_trough:last_trough));
plot(last_cycle);

SystolicPressure = max(last_cycle);
DiastolicPressure = last_cycle(end);
PulsePressure = SystolicPressure - DiastolicPressure;
MAP = trapz(last_cycle);
MAP = MAP/length(last_cycle);
r_r =  peakIndices(peakIdx(1)) - peakIndices(peakIdx(1) -1);
heart_rate = 60/(r_r/fs);

% find r_r peak using ecg
last_peak_ecg = find(sign(peakIndicesECG - peakIndices(peakIdx(1)))==-1);
last_peak_ecg_time = peakIndicesECG(last_peak_ecg(end)+1) - peakIndicesECG(last_peak_ecg(end));
heart_rate_ecg = 60/(last_peak_ecg_time/fs);

% display results
fprintf("Systolic Pressure is: %f \n",SystolicPressure)
fprintf("Diastolic Pressure is: %f \n",DiastolicPressure)
fprintf("Pulse Pressure is: %f \n",PulsePressure)
fprintf("Heart Rate ECG is: %f \n",heart_rate_ecg)
fprintf("Heart Rate is: %f \n",heart_rate)
fprintf("MAP (Integral) is: %f \n",MAP)


%% now for ISO
clear; close all;

% load in timing information and signals
load("TCDClean.mat")

% what is in the structure
TCDClean

% get the parts of the data you will need
BloodVelocity = TCDClean.V;
BloodPressure = TCDClean.BP;
ECG = TCDClean.ECG;
EEG = mean(TCDClean.EEG * 10e5,2);
ISOTime = TCDClean.ISO;

% bump ISOTime to the beginning of the final cycle
ISOTime = ISOTime + 80 + 88;


% plot the time series about the last time pulse pressure was less than
% 5mmHg
buffer = 6000; % 60 seconds

% extract 60 seconds pre and post CA
BloodVelocityISO = BloodVelocity(ISOTime - buffer: ISOTime + buffer);
BloodPressureISO = BloodPressure(ISOTime - buffer: ISOTime + buffer);
ECGISO = ECG(ISOTime - buffer: ISOTime + buffer);
EEGISO = EEG(ISOTime - buffer: ISOTime + buffer);

figure(1)
hold on
plot(BloodVelocityISO,'b');
plot(BloodPressureISO,'r');
plot(ECGISO,'g');
plot(EEGISO,'k');
xline(buffer)
%xticklabels(-buffer/100:(1/(buffer*2)):buffer/100)
%xlabel("Seconds")
hold off
% lets find the peaks
fs = 100/1.4;
WR = fs;
[peakValues, peakIndices] = findpeaks(BloodPressureISO,'MinPeakDistance',WR,'MinPeakHeight',60);

% Find troughs by negating the signal
[troughValues, troughIndices] = findpeaks(-BloodPressureISO,'MinPeakDistance',WR, 'MinPeakHeight',-45);
troughValues = -troughValues; % Negate back to get the actual trough values

% Visualization
figure(3);
plot(BloodPressureISO);
hold on;
plot(peakIndices, peakValues, 'ro'); % Peaks as red circles
plot(troughIndices, troughValues, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off

% do the same for ECG

[peakValuesECG, peakIndicesECG] = findpeaks(ECGISO,'MinPeakDistance',WR,'MinPeakHeight',mean(ECGISO));

% Find troughs by negating the signal
[troughValuesECG, troughIndicesECG] = findpeaks(-ECGISO,'MinPeakDistance',WR);
troughValuesECG = -troughValuesECG; % Negate back to get the actual trough values

% Visualization
figure(4);
plot(ECGISO);
hold on;
plot(peakIndicesECG, peakValuesECG, 'ro'); % Peaks as red circles
plot(troughIndicesECG, troughValuesECG, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off

% calculate systolic pulse pressure
% find index of blood pressure that precedes the CCA time (as defined by
% TCD)
peakIdx = find(sign(peakIndices - buffer) == 1);
troughIdx = find(sign(troughIndices - buffer) == 1);

% extract final cycle
last_cycle_time = peakIndices(peakIdx(1)-1):peakIndices(peakIdx(1) + 1);

% compute derivative of BP over this duration
last_cycle_dt = diff(BloodPressureISO(last_cycle_time));
last_cycle_dtdt = diff(last_cycle_dt);

% plot final cycle
figure(8)
plot(BloodPressureISO(last_cycle_time))
hold on
plot(last_cycle_dt)
plot(last_cycle_dtdt)
hold off

% trim the other cycles off the final cycle
% how many steps back from the peak do we need to go before having a
% negative derivative
last_peak_idx = find(BloodPressureISO(last_cycle_time) == peakValues(peakIdx(1)));
start_trough = find(sign(last_cycle_dt(1:last_peak_idx-1)) == -1);
start_trough = start_trough(end) + 1;

% same for end
last_trough_idx = find(BloodPressureISO(last_cycle_time) == peakValues(peakIdx(1)+1));
last_trough = find(sign(last_cycle_dt(1:last_trough_idx-1)) == -1);
last_trough = last_trough(end) + 1;

% now to calculate values over the final cycle
last_cycle = BloodPressureISO(last_cycle_time(start_trough:last_trough));
plot(last_cycle);

SystolicPressure = max(last_cycle);
DiastolicPressure = last_cycle(end);
PulsePressure = SystolicPressure - DiastolicPressure;
MAP = trapz(last_cycle);
MAP = MAP/length(last_cycle);
r_r =  peakIndices(peakIdx(1)) - peakIndices(peakIdx(1) -1);
heart_rate = 60/(r_r/fs);

% find r_r peak using ecg
last_peak_ecg = find(sign(peakIndicesECG - peakIndices(peakIdx(1)))==-1);
last_peak_ecg_time = peakIndicesECG(last_peak_ecg(end)+1) - peakIndicesECG(last_peak_ecg(end));
heart_rate_ecg = 60/(last_peak_ecg_time/fs);

% display results
fprintf("Systolic Pressure is: %f \n",SystolicPressure)
fprintf("Diastolic Pressure is: %f \n",DiastolicPressure)
fprintf("Pulse Pressure is: %f \n",PulsePressure)
fprintf("Heart Rate is: %f \n",heart_rate)
fprintf("Heart Rate ECG is: %f \n",heart_rate_ecg)
fprintf("MAP (Integral) is: %f \n",MAP)


