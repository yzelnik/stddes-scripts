%% define model parameters and such
diststr = 0.1;

Ps=struct('LocFunc',@L_Allee,'SpaFunc',@S_RD,'IntegFunc',@I_FDSIMP,'r',1,'K',1,'A',0.4,'gamma',4,'Ds',1,'VarNum',1,'Lx',500,'Ly',1,'Nx',1000,'Ny',1);
Es=struct('TsSize',0.05,'TsNum',2,'TimeDst',200,'SsThresh',1e-7,'NonNeg',1,'StSmall',0.01,'VarInd',1,'StAxis',[0 1.2],'TestFunc',@T_AvgVal,'ReachVal',0.99,'InitFunc',@M_InitUnfSt);

Es.ModPrm=[-diststr 1 diststr -1];


%%