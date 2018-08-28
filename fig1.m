%% show dynamics following a localized disturbance
modeldef;

tic;
frm1 = runflow(1,Ps,Es,'Es.TsSize',0.1,'Es.TsNum',10,'Es.OlDraw',1,'Es.FuncList',{@M_CutVar,@run2ss},'Es.ModPrm',[-diststr*2 1 diststr*2 -1]);
toc;

%%
distext = [0.1 1];
tottms= [120 3];
% run simulations for different extent of disturbances
for ii=1:length(distext)
    
    tmpmodprm=[diststr 1 distext(ii) 0.5];
    [frm1,bf1] = runflow(1,Ps,Es,'Es.TsSize',0.01,'Es.TsNum',2,'Es.OlDraw',0,'Es.FuncList',{@M_CutVar,@runframes},'Es.ModPrm',tmpmodprm,'Es.Frames',(0:120)*(tottms(ii)/120));
    frms{ii} = frm1;
    bfs{ii} = [[-tottms(ii)/5 1;-tottms(ii)/1e4 1] ; bf1]; % add info before disturbance
    rettime(ii)=C_ReachVal(bf1,Ps,Es); % return time
end;


%% plot overall trajectory with a series of snapshots

% time points (snapshots) to show
chs = {[0 40 80]+1,[0 30 60]+1};
snapnum = length(chs{1});
tjmp = [40 1];

figure('Position',[50 50 800 500])

clf;
% Setup panels
ha(1) = axes('Units','normalized', 'Position',[0.075 0.10 0.44 0.41]);
ha(2) = axes('Units','normalized', 'Position',[0.545 0.10 0.44 0.41]);
for ii=1:snapnum
    ha(ii+2)        =axes('Units','normalized', 'Position',[0.075+0.45*(ii-1)/snapnum 0.64 0.45/(snapnum+0.2) 0.3]);
    ha(ii+2+snapnum)=axes('Units','normalized', 'Position',[0.545+0.45*(ii-1)/snapnum 0.64 0.45/(snapnum+0.2) 0.3]);
end;

snapnum = length(chs{1});
ttls = {'Localized disturbance','Global disturbance'};

% go over different disturances response
for ii=1:2
    
    % plot snapshots
    for jj=1:snapnum 
        axes(ha((ii-1)*snapnum+jj+2));
        % plot state 
        plotst(frms{ii}(:,:,chs{ii}(jj)),Ps,Es,'Es.StLineWidth',1.5,'Es.St1Color',[0 0.7 0])
        hold on; plot([0 Ps.Lx],[1 1],'k--'); hold off; % add initial
        ylim([-0.01 1.1]);
        set(gca,'fontSize',14)
        if(jj==2) title(ttls{ii},'fontSize',24);  end;
        if(jj>1)||(ii>1) set(gca,'yTickLabel',[]); end;
        %if(jj==snapnum) xlabel('space','fontSize',16); else set(gca,'xTickLabel',[]); end;
        set(gca,'xTick',[250 500],'yTick',[0 0.5 1]);
        
        if(ii==1) numtxt = sprintf('t=%.0f',bfs{ii}(2+chs{ii}(jj),1));%round((chs{ii}(jj)-1)/10));
            else numtxt = sprintf('t=%.2f',bfs{ii}(2+chs{ii}(jj),1)); end;%round((chs{ii}(jj)-1)/10));
        text(Ps.Lx/20,0.1,numtxt,'fontSize',22,'color','r');
        if(jj==2) xlabel('space','fontSize',24); end;
    end;
    axes(ha(3));
    ylabel('biomass','fontSize',24);
    
    % plot overall biomass over time
    axes(ha(ii));
    plotbf(bfs{ii},Es,'Es.BfColor',[0 0 1],'Es.BfLineWidth',2);
    hold on; 
    plot(bfs{ii}([1 end],1),[1 1]*(1-0.1*diststr),'k:'); 
    plot(bfs{ii}(2+chs{ii},1),bfs{ii}(2+chs{ii},2),'or','markerSize',14);
    hold off;
    set(gca,'fontSize',14)
    xlabel('time','fontSize',24); 
    %title(ttls{ii},'fontSize',24);
    ylim([0.88 1.02]);
    xlim([-0.1 1]*tottms(ii))
    set(gca,'xTick',0:tjmp(ii):tottms(ii));
    if(ii==1) 
        set(gca,'yTick',0:0.05:1); 
        ylabel('total biomass','fontSize',24); 
    else
        set(gca,'yTick',[]);
    end;
   
end;

