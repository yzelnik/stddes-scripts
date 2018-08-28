% run simulations for 3 models for different values of disturbance extent
% (this should take a few minutes)
modeldef;

prmres = 101;
diststr = 0.2;
Es.ModPrm=[-diststr 1 diststr 0.5];
Es.TestFunc = {@T_AvgVal,@T_MinMax,[2 2],@T_LowReg};
Es.SegThresh=0.99;
bifrange=[diststr; diststr*3.; prmres; 0];
tsz=0.02;

tic;
[~,bf1]=runpar(1,Ps,Es,'Es.TsSize',tsz,'Ps.LocFunc',@L_Log,'Es.InitFunc',@M_InitRndSt,'Es.BfRange',bifrange,'Es.BfPrm','Es.ModPrm(3)','Es.FuncList',{@M_CutVar,@run2ss});
toc;
tic;
[~,bf2]=runpar(1,Ps,Es,'Es.TsSize',tsz/2,'Ps.LocFunc',@L_SR,'Es.BfRange',bifrange,'Es.BfPrm','Es.ModPrm(3)','Es.FuncList',{@M_CutVar,@run2ss}); 
toc;

tic;
[~,bf3]=runpar(1,Ps,Es,'Es.TsSize',tsz,'Ps.LocFunc',@L_Allee,'Es.BfRange',bifrange,'Es.BfPrm','Es.ModPrm(3)','Es.FuncList',{@M_CutVar,@run2ss});
toc;

%% extract trajectories
bfs={bf1,bf2,bf3};
columns=4;
ver=[]; hor=[]; retm=[];

% run over different model runs
for modind=1:length(bfs)
    bftmp=bfs{modind};
    smvar=[]; smavg=[]; 
    
    dx = diff(bftmp(1,3:4));
    for ii=1:size(bftmp,1)
        trjlen = find(abs(diff(diff(bftmp(ii,2:end))))>1e-10,1)+1;
        trj    = reshape(bftmp(ii,1+(1:trjlen*columns)),trjlen,columns);
        xx=trj(:,1);
        yy=1-trj(:,2);
        
        retm(ii,modind)=sum(yy>diststr/100)*dx;
        ver{modind}{ii}=trj(:,3);
        
        [~,ind]=min(abs(trj(:,3)-Es.SegThresh)); 
        vecdif = diff(trj(ind-2:ind-1,4));
        hortmp = [trj(1:ind-1,4) ; max(0,(1:size(trj,1)-ind)'*vecdif+trj(ind-1,4)); 0];
       
        hor{modind}{ii}=hortmp;
        tim{modind}{ii}=xx;
        rec{modind}{ii}=yy;
    end;
end;
retax=bfs{1}(:,1);

%% calculate local potential
r1=1; r2=1; r3=1;
A=0.4;
K=1;
gam=4;

N = 0:0.001:1.4*K; % biomass variable

% local potentials for 3 models
ef1 = -r1*(N.^2)/2+r1*(N.^3)/(3*K);
ef2 = -r2*((N.^(gam+2))/(K.^gam)).*(1/(gam+2)-(N/K)/(gam+3));
ef3 = -(-r3*(N.^2)/2+r3*(K+A)*(N.^3)/(3*K*A)-r3*(N.^4)/(4*K*A));


%% plot it all out

efs = {ef1,ef2,ef3};
jmps = [10 100 200];
ylms = [-0.2 0.05 ; -0.03 0.02; -0.05 0.04];
rmax = [22 220 420];
yext = [0.014 0.016; 0.0020 0.0028 ;0.005 0.008];

trjnum=12; % number of trajectories to show

% choose the trajectories 
chs = round((0.5:trjnum)*(prmres/trjnum));
x=0:0.01:1;

% setup plot
clf;
ha = tight_subplot(3,3,[0.03 0.08],[0.09 0.04],[0.045 0.015]);

focus = [2 length(chs)-1];

% go over 3 models
for ii=1:3

%biomass values on trajectories
bvals(1,:) = ver{ii}{chs(focus(1))}([1 100]);
bvals(2,:) = ver{ii}{chs(focus(2))}([1 100]); 
uvals = bvals;
for jj=1:length(bvals(:))
    [~,minind]=min(abs(N-bvals(jj)));
    uvals(jj)=efs{ii}(minind)+yext(ii,mod((jj-1),2)+1);
end;

% plot local potential
axes(ha(ii*3-2));
plot(N,efs{ii},'k','LineWidth',2);
ylim(ylms(ii,:));
xlim([0 1.4]);
set(gca,'fontSize',14);

if(ii==3) xlabel('biomass value','fontSize',24); set(gca,'XTick',0:0.5:2); else set(gca,'xTickLabel',[]); end;

% plot circles on potentials
set(gca,'YTick',[]);
hold on;
plot(bvals(1,1),uvals(1,1),'.m','markerSize',36)
plot(bvals(2,1),uvals(2,1),'.g','markerSize',36)
hold off;


% plot recovery trajectories
% start with the boundaries
axes(ha(ii*3-1));
plot(x,1-diststr./x,'k--','lineWidth',2)  
hold on; 
patch([x fliplr(x)], [ max(0,1-diststr*0.1./x) ones(size(x))*1.02], 'k','FaceAlpha',.5,'edgeColor',[1 1 1],'edgeAlpha',0)
text(0.035,0.79,'Recovered','fontSize',11,'Rotation',30,'color','w');
set(gca,'fontSize',14)

% plot each trajectory
for jj=1:length(chs)
    plot(hor{ii}{chs(jj)},ver{ii}{chs(jj)},'b-','lineWidth',1.5); 
    axis([0 0.615 0 1.002]);
end; 
% plot circles on trajectory
plot(hor{ii}{chs(focus(1))}(1),bvals(1,1),'om','markerSize',16,'lineWidth',2)
plot(hor{ii}{chs(focus(2))}(1),bvals(2,1),'og','markerSize',16,'lineWidth',2)
hold off;

if(ii==3) xlabel('size of disturbed region','fontSize',24); set(gca,'xTick',(0:6)/5); else set(gca,'xTick',[]); end;
if(ii==2) ylabel('local biomass level','fontSize',24); end;
set(gca,'yTick',(0:2)/2,'xTickLabel',(0:6)/5,'yTickLabel',(0:2)/2)

% plot return time curves
axes(ha(ii*3));
plot(retax,retm(:,ii),'r',retax(chs),retm(chs,ii),'bx','lineWidth',2,'markerSize',10); %ylim([0 20]);
hold on;
plot(retax(chs(focus(1))),retm(chs(focus(1)),ii),'om','markerSize',16,'lineWidth',2);
plot(retax(chs(focus(2))),retm(chs(focus(2)),ii),'og','markerSize',16,'lineWidth',2);
hold off;
set(gca,'fontSize',14)
if(ii==3) xlabel('disturbance extent','fontSize',24); else set(gca,'xTickLabel',[]); end;

xlim([0.2 0.6])
ylim([0 rmax(ii)]);
set(gca,'yTick',(0:10)*jmps(ii),'xTick',(0:10)/5);

end;

% add titles and model names
ttls={'local potential','recovery trajectories','return time'};
models = {'Logistic Growth','Slow Recovery','Allee Effect'};    
for ii=1:length(ttls)
    axes(ha(ii));
    title(ttls{ii},'fontSize',22);
end;
for ii=1:length(models)
    axes(ha(1+(ii-1)*length(ttls)))
    ylabel(models{ii},'fontSize',24); 
end;

