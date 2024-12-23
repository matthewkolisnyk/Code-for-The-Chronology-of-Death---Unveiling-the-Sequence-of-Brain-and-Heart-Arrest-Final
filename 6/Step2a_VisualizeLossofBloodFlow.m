%% Identify critical points within the signal
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


% drop time of loss of signal
loss_time = 525796;
ZeroFlow = ZeroFlow(ZeroFlow > loss_time);
% remove ZeroFlow values at the time of last QRS
ZeroFlow = ZeroFlow(ZeroFlow < last_QRS);

% plot from the first time there was zero flow
figure
plot(V(ZeroFlow(1):last_QRS))

% measure the difference in times between zero flow measures. Find when the
% points seem to converge
ZeroInterval = diff(ZeroFlow);
[IntervalIndex, IntervalResidual] = findchangepts(ZeroInterval,MaxNumChanges=10);
[IntervalIndexMax, IntervalResidualMax] = findchangepts(ZeroInterval);

% visualize the fit
figure
plot(V);
xline(ZeroFlow(1));
xline(ZeroFlow(IntervalIndexMax))
xline(last_PP)

% iterate over change points until residual not changing 
nIter = 50; RESIDUALS = []; INDEX = {};
for i=1:nIter
    [IntervalIndexIterated, IntervalResidualIterate] = findchangepts(ZeroInterval,'Statistic', 'mean',MaxNumChanges=i);
    RESIDUALS = [RESIDUALS,IntervalResidualIterate];
    INDEX{i} = IntervalIndexIterated;
end


% where dt stops changing 
plot(RESIDUALS)
stop_index = 9;
ZeroIntervalStop = INDEX{stop_index}(end);

% visualize the fit
figure
plot(V);
xline(ZeroFlow(1),'r');
xline(ZeroFlow(IntervalIndexMax),'g')
xline(ZeroFlow(ZeroIntervalStop),'b')
xline(last_PP)

% Let's add it to the TCDClean structure
data.ZeroFlow = ZeroFlow(IntervalIndexMax);
data.ZeroFlowResid = ZeroFlow(ZeroIntervalStop);

% Add your recorded time as well.
% 667936 is likely biphasic oscillatory flow but there are a few case with
% spikes above 50. 
% 671309 is unequivocally spike activity < 50
data.CCA_Matt = 671309;

TCDClean = data;
TCDClean.CriticalPoints = [TCDClean.CriticalPoints,TCDClean.ZeroFlowResid,data.CCA_Matt];

% resave the structure
save("TCDClean.mat","TCDClean")

