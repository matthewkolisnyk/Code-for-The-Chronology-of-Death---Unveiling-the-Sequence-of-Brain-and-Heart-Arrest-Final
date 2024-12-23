%% visualize and describe timeseries during baseline period 
% note CCA = loss of brain blood and CA = loss of systemic circulation
% load in the aligned and timing data
data = load("TCDClean.mat").TCDClean;
last_QRS = data.QRS; last_PP = data.PP;

% clean data 

% start with identifying the first time there is zero flow
V = data.V;
BP = data.BP;

Time = data.Time;
ZeroFlow = find(V <= 0);
plot(V); 
hold on
plot(BP)
hold off


% calculate measure from a typical part of the baseline period
clean_index = 1000:6000;

figure
plot(V(clean_index))
hold on
plot(BP(clean_index))
hold off

baselineBP = BP(clean_index);
baselineV = V(clean_index);

fs = 100;
WR = fs * 0.45;
[peakValues, peakIndices] = findpeaks(baselineBP,'MinPeakDistance',WR,'MinPeakHeight',mean(baselineBP)+10);

% Find troughs by negating the signal
[troughValues, troughIndices] = findpeaks(-baselineBP,'MinPeakDistance',WR, 'MinPeakHeight',-mean(baselineBP));
troughValues = -troughValues; % Negate back to get the actual trough values

% Visualization
figure(3);
plot(baselineBP);
hold on;
plot(peakIndices, peakValues, 'ro'); % Peaks as red circles
plot(troughIndices, troughValues, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off

% for each cycle
% drop the first peak
%peakValues = peakValues(2:end); peakIndices = peakIndices(2:end);
% drop the last trough
troughValues = troughValues(1:end-1); troughIndices = troughIndices(1:end-1);

% number of peaks
numPeaks = length(peakValues);

bMAP = []; bPP = []
% for each peak find the cycle
for np = 1:numPeaks-1 % skip first and last cycle
    % extract final cycle
    cycle_time = troughIndices(np):peakIndices(np + 1);

    % compute derivative of BP over this duration
    cycle_dt = diff(baselineBP(cycle_time));
    cycle_dtdt = diff(cycle_dt);

    % how many steps back from the peak do we need to go before having a
    % negative derivative
    start_trough = find(sign(cycle_dt) == -1);
    start_trough = start_trough(end);

    % compute map on this cycle
    cycle = baselineBP(cycle_time(1:start_trough));
    MAP = trapz(cycle);
    MAP = MAP/length(cycle);
    bMAP = [bMAP,MAP];

    % compute pulse pressure on this cycle
    bPP = [bPP max(cycle) - min(cycle)];
end

fprintf("MAP (Integral) is: %f \n",median(bMAP))
fprintf("Pulse Pressure is: %f \n",median(bPP))


%% same for V
fs = 100;
WR = fs * 0.3;
[peakValues, peakIndices] = findpeaks(baselineV,'MinPeakDistance',WR,'MinPeakHeight',mean(baselineV)+30);

% Find troughs by negating the signal
[troughValues, troughIndices] = findpeaks(-baselineV,'MinPeakDistance',WR,'MinPeakHeight',-mean(baselineV));
troughValues = -troughValues; % Negate back to get the actual trough values

% Visualization
figure();
plot(baselineV);
hold on;
plot(peakIndices, peakValues, 'ro'); % Peaks as red circles
plot(troughIndices, troughValues, 'go'); % Troughs as green circles
legend('Data', 'Peaks', 'Troughs');
title('Peak and Trough Detection');
hold off

% for each cycle
% drop the first peak
peakValues = peakValues(2:end); peakIndices = peakIndices(2:end);
% drop the last trough
troughValues = troughValues(1:end-2); troughIndices = troughIndices(1:end-2);

% number of peaks
numPeaks = length(peakValues);

bSV = []; bDV = [];
% for each peak find the cycle
for np = 1:numPeaks-1 % skip first and last cycle
    % extract final cycle
    cycle_time = troughIndices(np):peakIndices(np + 1);

    % compute derivative of BP over this duration
    cycle_dt = diff(baselineV(cycle_time));
    cycle_dtdt = diff(cycle_dt);

    % how many steps back from the peak do we need to go before having a
    % negative derivative
    start_trough = find(sign(cycle_dt) == -1);
    start_trough = start_trough(end);

    % compute map on this cycle
    cycle = baselineV(cycle_time(1:start_trough));
    bSV = [bSV,max(cycle)];

    % compute pulse pressure on this cycle
    bDV = [bDV  min(cycle)];
end



% display results
fprintf("Systolic Velocity is: %f \n",median(bSV))
fprintf("Diastolic Velocity is: %f \n",median(bDV))
