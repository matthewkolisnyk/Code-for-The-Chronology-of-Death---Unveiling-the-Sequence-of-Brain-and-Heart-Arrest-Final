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
CCATime = TCDClean.CCA_Matt;

% shift CCA Time to be at the beginning of the final cycle
CCATime = CCATime + 1;

% plot the time series about the last time pulse pressure was less than
% 5mmHg
buffer = 6000; % 60 seconds

% extract 60 seconds pre and post CA
BloodVelocityCCA = BloodVelocity(CCATime - buffer: CCATime + buffer);
BloodPressureCCA = BloodPressure(CCATime - buffer: CCATime + buffer);
ECGCCA = ECG(CCATime - buffer: CCATime + buffer);

figure(1)
hold on
plot(BloodVelocityCCA,'b');
plot(BloodPressureCCA,'r');
plot(ECGCCA,'g');
xline(buffer)
%xticklabels(-buffer/100:(1/(buffer*2)):buffer/100)
%xlabel("Seconds")
hold off

% lets find the peaks
fs = 100;
WR = fs * 1.1;
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

% calculate systolic pulse pressure
% find index of blood pressure that precedes the CCA time (as defined by
% TCD)
peakIdx = find(sign(peakIndices - buffer) == -1);
troughIdx = find(sign(troughIndices - buffer) == -1);

afterPeak = find(sign(peakIndices - buffer) == 1);
afterTrough = find(sign(troughIndices - buffer) == 1);

peakIdx = [peakIdx, afterPeak(1)]; 
troughIdx = [troughIdx, afterTrough(1)];

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
heart_rate = 60/(r_r/fs);

% display results
fprintf("Systolic Pressure is: %f \n",SystolicPressure)
fprintf("Diastolic Pressure is: %f \n",DiastolicPressure)
fprintf("Pulse Pressure is: %f \n",PulsePressure)
fprintf("Heart Rate is: %f \n",heart_rate)
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
CATime = CATime + 72;


% plot the time series about the last time pulse pressure was less than
% 5mmHg
buffer = 6000; % 60 seconds

% extract 60 seconds pre and post CA
BloodVelocityCCA = BloodVelocity(CATime - buffer: CATime + buffer);
BloodPressureCCA = BloodPressure(CATime - buffer: CATime + buffer);
ECGCCA = ECG(CATime - buffer: CATime + buffer);

figure(1)
hold on
plot(BloodVelocityCCA,'b');
plot(BloodPressureCCA,'r');
plot(ECGCCA,'g');
xline(buffer)
%xticklabels(-buffer/100:(1/(buffer*2)):buffer/100)
%xlabel("Seconds")
hold off

% lets find the peaks
fs = 100;
WR = fs/(1.25);
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

% calculate systolic pulse pressure
% find index of blood pressure that precedes the CCA time (as defined by
% TCD)
peakIdx = find(sign(peakIndices - buffer) == 1);
troughIdx = find(sign(troughIndices - buffer) == 1);

% extract final cycle
last_cycle_time = peakIndices(peakIdx(1)-1):peakIndices(peakIdx(1) + 1);

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
last_peak_idx = find(BloodPressureCCA(last_cycle_time) == peakValues(peakIdx(1)));
start_trough = find(sign(last_cycle_dt(1:last_peak_idx-1)) == -1);
start_trough = start_trough(end) + 1;

% same for end
last_trough_idx = find(BloodPressureCCA(last_cycle_time) == peakValues(peakIdx(1)+1));
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
r_r =  peakIndices(peakIdx(1)) - peakIndices(peakIdx(1) -1);
heart_rate = 60/(r_r/fs);

% display results
fprintf("Systolic Pressure is: %f \n",SystolicPressure)
fprintf("Diastolic Pressure is: %f \n",DiastolicPressure)
fprintf("Pulse Pressure is: %f \n",PulsePressure)
fprintf("Heart Rate is: %f \n",heart_rate)
fprintf("MAP (Integral) is: %f \n",MAP)


