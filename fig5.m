%% define model and run simulations

modeldef; % run script for defining model

% number of points in disturbance space to simulate
pointnum= 50;
% define random points in the disturbance parameter space 

% uniform distribution, without extrme values around 0 and 1
basernd = rand(pointnum,2,2)*0.95+0.025; 
extval = basernd(:,:,1); % extent
strval = basernd(:,:,1).*basernd(:,:,2); %strength

% define parameter space
Es.BfPrm={'Es.ModPrm(1)','Es.ModPrm(3)'};
% functions to run per simulation: disturbance, time-integration, check time of recovery
Es.FuncList={@M_CutVar,@run2ss,@C_ReachVal};

rng(1); % set randomization seed

% run simulations for the Slow Recovery model
tic;
[~,bf1]=runpar(1,Ps,Es,'Ps.LocFunc',@L_SR,'Es.BfRange',[-strval(:,1) extval(:,1)]);
toc;
% run simulations for the Allee Effect model
tic;
[~,bf2]=runpar(1,Ps,Es,'Ps.LocFunc',@L_Allee,'Es.BfRange',[-strval(:,2) extval(:,2)]);
toc;

% get out the return times
rettms = [bf1(:,3) bf2(:,3)];
%% plot out the reconstructed curves

modelnms = {'SR','AE'};
for ii=1:2
    % plot the points in the parameter space
    subplot(2,2,ii*2-1);
    plot(extval(:,ii),strval(:,ii),'*r',[0 1 1 0],[0 1 0 0],'k'); 
    axis([0 1 0 1])
    title(sprintf('parameter space: %s',modelnms{ii}),'fontSize',20);
    xlabel('disturbance extent')
    ylabel('disturbance strength')
    
    % plot out the points as a reconstructed curve
    subplot(2,2,ii*2);
    plot(1-strval(:,ii)./extval(:,ii),rettms(:,ii)./strval(:,ii),'*r'); 
    title(sprintf('reconstructed curve: %s',modelnms{ii}),'fontSize',20);
    xlabel('normalized disturbance extent')
    ylabel('normalized return time')
   %plot(strval(:,ii),rettms(:,ii),'*r'); 
end;

