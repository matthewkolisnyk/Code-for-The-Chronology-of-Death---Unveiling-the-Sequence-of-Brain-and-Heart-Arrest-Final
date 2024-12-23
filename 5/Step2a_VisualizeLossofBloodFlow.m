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

% re-identify ZeroFlow after removing lost signal
clean_index = 217473;
ZeroFlow = ZeroFlow(ZeroFlow > clean_index);

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
stop_index = 7;
ZeroIntervalStop = INDEX{stop_index}(end);

% visualize the fit
figure
plot(V);
xline(ZeroFlow(1),'r');
xline(ZeroFlow(IntervalIndexMax),'g')
xline(ZeroFlow(ZeroIntervalStop),'blue')
xline(last_PP)

% Let's add it to the TCDClean structure
data.ZeroFlow = ZeroFlow(IntervalIndexMax);
data.ZeroFlowResid = ZeroFlow(ZeroIntervalStop);

% Add your recorded time as well.
data.CCA_Matt = 229500;

TCDClean = data;
TCDClean.CriticalPoints = [TCDClean.CriticalPoints,TCDClean.ZeroFlowResid,data.CCA_Matt];

% resave the structure
save("TCDClean.mat","TCDClean")

