%% setup model with 2-dimensional space
modeldef;
% redefine the model as 2-d space (size 200x200)
Ps.Nx=100;
Ps.Ny=100;
Ps.Lx=200; 
Ps.Ly=200;
% define the local dynamics to be of the Slow Recovery model
Ps.LocFunc = @L_SR;
Ps.gamma=2;
% define the disturbance extent to be a bit larger then the strength
Es.ModPrm=[-0.1 1 0.11 -1];

%% run short sim (a single disturbance)

[frm1,bf1] = runframes(1,Ps,Es,'Es.RecurFrames',[1 1e5],'Es.TsSize',0.05,'Es.OlDraw',0,'Es.RecurFunc',@M_CutVar,'Es.Frames',0:0.1:60,'Es.ModPrm',[-0.1 1 0.11 0.5]);

%% run long sim (multiple disturbances)
randind = 5; 
tottime = 800;  % simulation time
dstfreq = 0.04; % frequency of disturbances
res     = 0.05; % time resolution
frmjmp  = 40;   % every how many frames (snapshots) to save one

% define a series of disturbance times (from an exponential distribution)
rng(randind);
tmptms = cumsum(-1/dstfreq*log(rand(round(2*tottime*dstfreq),1)));
disturbtms = unique([1; ceil(tmptms/res)]);

% run this simulation
rng(randind);
tic;
[stlong,bflong] = runframes(1,Ps,Es,'Es.RecurFrames',disturbtms,'Es.OlDraw',0,'Es.RecurFunc',@M_CutVar,'Es.Frames',0:res:tottime,'Es.FramesChoice',1:frmjmp:1e5);
toc;

%% a simple figure (Fig. 3)

pnts = [200 158]; % which times points to show?

% create time-series from the trajectory following a single disturbance
plotres = 2; 
tmpbf = [[(-20:plotres:-1)'; bf1(1:10*plotres:end,1);bf1(end,1)+(1:plotres:20)'] ,[ones(20/plotres,1); bf1(1:10*plotres:end,2);ones(20/plotres,1)]];
% add some noise to time-series
tmpbf(:,2) = tmpbf(:,2)+rand(size(tmpbf(:,2)))/100-0.005;

% plot trajectory
subplot(2,2,1);
rng(1)
plot(tmpbf(:,1) ,tmpbf(:,2),'b.-','lineWidth',1.25,'markerSize',16);

% add other info (red dot, etc...)
hold on;
patch([bf1(:,1)' fliplr(bf1(:,1)')], [ones(1,length(bf1)) fliplr(bf1(:,2)')], 'k','FaceAlpha',.15,'edgeAlpha',0)
plot(tmpbf(:,1) ,tmpbf(:,2),'b.-','lineWidth',1.25,'markerSize',16);
plot(bf1(pnts(1),1),bf1(pnts(1),2),'r.','markerSize',48);
hold off;

set(gca,'fontSize',14);
set(gca,'yTick',0:0.05:1);
ylim([0.86 1.025])
set(gca,'xTick',0:20:80,'yTick',0:0.05:2);
xlabel('time','fontSize',24); ylabel('total biomass','fontSize',24);
xlim([-20 70])

% plot a 2-d snapshot (for the single disturbance)
subplot(2,2,2);
plotst(M_ShiftSt(frm1(:,:,pnts(1)),Ps,Es,'Es.ShiftPrm',[0.49 0.]),Ps,Es)
set(gca,'xTick',[],'yTick',[]);
xlabel('space','fontSize',24); ylabel('space','fontSize',24);


plotres=8;
subplot(2,2,3);
rng(1);
% create time-series from the trajectory due to multiple disturbances
tmpbf = bflong(1:plotres*20:end,:);
% add some noise to time-series
tmpbf(:,2) = tmpbf(:,2)+rand(size(tmpbf(:,2)))/50-0.01;

% plot time-series
plot(tmpbf(:,1),tmpbf(:,2),'b.-','lineWidth',1.25,'markerSize',16);
set(gca,'fontSize',14);
set(gca,'yTick',0:0.2:1);
xlabel('time','fontSize',24); ylabel('total biomass','fontSize',24);
ylim([0.65 1.02])
set(gca,'xTick',0:200:800)

% add variability measure
place=175;
avg=mean(bflong(:,2));
wid=std(bflong(:,2));
hold on
patch([bflong([1 end],1)' fliplr(bflong([1 end],1)')]-bflong(1,1), [repmat(avg-wid,1,2) repmat(avg+wid,1,2)], 'k','FaceAlpha',.15,'edgeColor',[1 1 1],'edgeAlpha',0)
plot(tmpbf(:,1),tmpbf(:,2),'b.-','lineWidth',1.25,'markerSize',16);
plot(bflong((pnts(2)-1)*frmjmp,1),bflong((pnts(2)-1)*frmjmp,2),'r.','markerSize',48);
hold off

% plot a 2-d snapshot (for the multiple-disturbances)
subplot(2,2,4);
plotst(M_ShiftSt(stlong(:,:,pnts(2)),Ps,Es,'Es.ShiftPrm',[0.49 0.87]),Ps,Es)
set(gca,'xTick',[],'yTick',[]);
xlabel('space','fontSize',24); ylabel('space','fontSize',24);

