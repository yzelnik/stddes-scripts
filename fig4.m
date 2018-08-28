%% Note: This script runs simulations to reproduce results similar to that shown in Fig. 4 (and also some figures in the appendix).
%  However, as these take a very long time, we demonstrate here a somewhat 
%  shorter version, that might take 1 hour to run (instead of days)
%  In particular, to make the run more similar to the results shown in Fig.4, you should:
%  - increase diststr to 0.1
%  - increase randnum to 100
%  - increase totime to either 10000 (for variability) or 100000 (for collapse probability)
%  - change freq to the appropiate frequency is question
%  - change Ps.LocFunc to @L_Log or @L_SR for the LG and SR models respectively (instead of the default AE model)

%% define model and simulation setup

modeldef; % run script to define model parameters

freq    = 0.15; % frequency of disturbances
tottime = 1000; % time of simulation
randnum = 20;   % how many randomizations
bfres   = 11;   % how many points along sigma-axis
diststr = 0.05; % disturbance strength and also its smallest extent
distmax = 2;    % ratio between smallest and largest disturbance extent

% setup the poisson process (for disturbances), run simulation, calculate variance and avg. biomass
Es.FuncList = {@U_SetupPoissonProcess,@runframes,@C_CalcVariance,@T_AvgVal};

Es.TsSize   = 0.1; % time step size of time integration
Es.Frames   = 0:Es.TsSize:tottime; % define simulation time-points
Es.RecurFunc= @M_CutVar;  % define recurrent function (making the disturbances)
Es.ModPrm   = [-diststr 1 diststr -1]; % [dist-strength variable-index dist-extent random-dist-location]
Es.CalcRange= 0.2; % Calculate temporal variance starting from 20% of elapsed time into the simulation
Es.ReachVal = Ps.K*-0.9*diststr; % value of biomass to reach that is considered recovered

% define the parameters to change for the simulation: the dist. extent (sigma) and the randomization seed
Es.BfPrm   = {'Es.ModPrm(3)','Es.RandSeed'};
Es.BfRange = {[diststr; diststr*distmax; bfres; 0],[1; randnum; randnum; 0]};

Vs = Ps.K; % initial biomass value

%% Run simulations of recovery from single disturbance - to calculate return time and the approx. of V=f*int(g)
retfuncs  = {@M_CutVar,@run2ss,@C_MultiCalc};	
retcalcs  = {@C_CalcMoments,@C_ReachVal};
tic;
[~,bf0]=runpar(Vs,Ps,Es,'Es.BfPrm',Es.BfPrm{1},'Es.BfRange',Es.BfRange{1},'Es.FuncList',retfuncs,'Es.CalcList',retcalcs);
toc;
%% Run the simulations of multiple disturbances.
tic;
% Note: You can set Es.Verbose=0 to not print out info while running.
[~,bf1]=runpar(Vs,Ps,Es,'Es.Verbose',1,'Es.PppPrm',1/freq);
toc;

%% plot out return time, variability and collapse probability
axsigma = bf0(:,1);
moment2 = bf0(:,3);
rettime = bf0(:,4);
vardata = reshape(bf1(:,3),bfres,randnum);

subplot(1,3,1);
% plot the return time
plot(axsigma,rettime,'-*r')
ylabel('return time','fontSize',20);
xlim(diststr*[1 distmax])

subplot(1,3,2);
% plot the average (and error-estimation) of the variability
errorbar(axsigma,mean(vardata'),std(vardata'));
hold on;
% add the approximation using the second moment of the response function g(t)
plot(axsigma,moment2*freq,'k')
hold off;
xlim(diststr*[1 distmax])
ylabel('variability','fontSize',20);
xlabel('disturbance extent','fontSize',20);

subplot(1,3,3);
% plot the collapse probability (note that here it us using the same time frame tau)
plot(axsigma,mean(reshape(bf11(:,4)<1e-2,bfres,randnum)'),'-*g')
ylabel('collapse probability','fontSize',20);
xlim(diststr*[1 distmax])
ylim([0 1])
