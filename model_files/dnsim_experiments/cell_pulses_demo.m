cd /home/jason/models/dnsim/SparsePING/base;
% ---------------------------------------------------------
% Parameters
nE=80; nI=20; % nE=40; nI=10;
tauI=10; gI=.1; gE=.1; stim=7.5; noise=4;
% ---------------------------------------------------------
spec=[];
spec.nodes(1).label = 'E';
spec.nodes(1).multiplicity = nE;
spec.nodes(1).dynamics = {'V''=(current)./Cm'};
spec.nodes(1).mechanisms = {'K','Na','leak'};
spec.nodes(1).parameters = {'Cm',1,'V_IC',-70,'stim',stim,'noise',noise};
spec.nodes(2).label = 'I';
spec.nodes(2).multiplicity = nI;
spec.nodes(2).dynamics = {'V''=(current)./Cm'};
spec.nodes(2).mechanisms = {'K','Na','leak'};
spec.nodes(2).parameters = {'Cm',1,'V_IC',-70,'stim',0,'noise',noise};

spec.connections(1,2).label = 'E-I';
spec.connections(1,2).mechanisms = {'AMPA'};
spec.connections(1,2).parameters = {'tauDx',2,'g_SYN',gE*(80/nE)};
spec.connections(2,1).label = 'I-E';
spec.connections(2,1).mechanisms = {'GABAa'};
spec.connections(2,1).parameters = {'tauDx',tauI,'g_SYN',gI*(20/nI)};

% [allX,stats] = cell_pulses(spec,'plot_flag',1,'minamp',-7,'maxamp',7,'damp',1,'pulsetime',200,'dt',.01);
[allX,stats] = cell_pulses(spec,'plot_flag',1,'minamp',-15,'maxamp',15,'damp',.25,'pulsetime',200,'dt',.01);
figure; cellind=1;
subplot(2,1,1); plot(stats.amps,stats.Veq(:,cellind),'.-'); xlabel('I'); ylabel('V'); title(sprintf('Reff=%g',stats.Reff(cellind)));
subplot(2,1,2); plot(stats.amps,stats.Rss(:,cellind),'.-'); xlabel('I'); ylabel('f'); title(sprintf('fIslope=%g',stats.fIslope(cellind)));

[allMUA,netstats] = net_drives(spec,'targets',{'E','I'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);

%%
% get network spec for task simulations
cd /home/jason/models/dnsim/aro/base
ruletask_ResponseLR_CB_competition;

% remove inputs and noise sources
exclude={'InputGenerator2','randn'};
for i=1:length(spec.nodes)
  spec.nodes(i).dynamics={'V''=(current)./Cm'};
  spec.nodes(i).mechanisms=setdiff(spec.nodes(i).mechanisms,exclude);
end

[allX,stats] = cell_pulses(spec,'minamp',-15,'maxamp',15,'damp',1,'pulsetime',200,'plot_flag',1);

[allMUA1,netstats1] = net_drives(spec,'targets',{'E1','Ifast1'},...
  'minamp',-15,'maxamp',15,'damp',1,'pulsetime',200,'plot_flag',1);

[allMUA2,netstats2] = net_tuning(spec,'targets',{'E1'},...
  'minamp',1,'maxamp',2,'damp',1,'minfreq',5,'maxfreq',60,'dfreq',2.5,...
  'pulsetime',1000,'baseline',50,'plot_flag',1);
  % todo: add filter metrics (fc,BW) and support for multiple targets

% probe network
[allMUA,netstats] = net_drives(spec,'targets',{'E1','Ifast1'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E1','Islow1'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E1','Ifast2'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E1','Islow2'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E2','Ifast2'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E2','Islow2'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E1','E2'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);
[allMUA,netstats] = net_drives(spec,'targets',{'E1','Ifast1','Islow1'},'maxamp',15,'damp',2,'pulsetime',200,'plot_flag',1);

[allX,stats] = cell_pulses(spec,'plot_flag',1,'minamp',-15,'maxamp',15,'damp',1,'pulsetime',200,'dt',.01);
figure; cellind=1;
subplot(2,1,1); plot(stats.amps,stats.Veq(:,cellind),'.-'); xlabel('I'); ylabel('V'); title(sprintf('Reff=%g',stats.Reff(cellind)));
subplot(2,1,2); plot(stats.amps,stats.Rss(:,cellind),'.-'); xlabel('I'); ylabel('f'); title(sprintf('fIslope=%g',stats.fIslope(cellind)));

%% Task model: cell_pulses()
cd /home/jason/models/dnsim/aro/base
ruletask_ResponseLR_CB_competition;
% remove inputs and noise sources
exclude={'InputGenerator2','randn'};
for i=1:length(spec.nodes)
  spec.nodes(i).dynamics={'V''=(current)./Cm'};
  spec.nodes(i).mechanisms=setdiff(spec.nodes(i).mechanisms,exclude);
end
% store modified specification
spec0=spec; % full model
% ----------------------

specEsup=spec;
specEsup.nodes=specEsup.nodes(1);
spec=specEsup;

% test: gpas=.04 (taum=25ms|Cm=1)

% parms.minamp=-1; 
% parms.maxamp=1; 
% parms.damp=.1;
% parms.pulsetime=1000;
% parms.baseline=500;
% parms.override={spec.nodes(1).label,'parameters','Cm',2};

[allX,stats] = cell_pulses(spec,'minamp',-1,'maxamp',1,'damp',.1,'pulsetime',1000,'baseline',500,'plot_flag',1);

cellind=1;

% Cell passive biophysics:
% effective membrane resistance from I/V curve
Reff=stats.Reff(cellind);
% effective membrane time constant from voltage relaxation to (1/e)Vmax
taum_eff=stats.taum_median(cellind);
% leak resistance from model specification
i=find(cellfun(@(x)isequal(x,'gpas'),stats.model.nodes(cellind).parameters));
Rleak=1/stats.model.nodes(cellind).parameters{i+1};
% membrane capacitance from model specification
i=find(cellfun(@(x)isequal(x,'Cm'),stats.model.nodes(cellind).parameters));
Cm=stats.model.nodes(cellind).parameters{i+1};
% calc model membrane time constant
taum_model=Rleak*Cm;

fprintf('Model:      Rpas=%3.3g, taum=%3.3g\n',Rleak,taum_model);
fprintf('Experiment: Rpas=%3.3g, taum=%3.3g\n',Reff,taum_eff);
%{
Cm=1
  Model:      Rpas=58.8, taum=58.8
  Experiment: Rpas=59.1, taum=58.6
Cm=2
  Model:      Rpas=58.8, taum=118
  Experiment: Rpas=58.9, taum=117
%}
figure;
subplot(3,1,1); plot(stats.amps,stats.Veq(:,cellind),'.-'); xlabel('I'); ylabel('V'); title(sprintf('Reff=%g',stats.Reff(cellind)));
subplot(3,1,2); plot(stats.amps,stats.Rss(:,cellind),'.-'); xlabel('I'); ylabel('f'); title(sprintf('fIslope=%g',stats.fIslope(cellind)));
subplot(3,1,3); plot(stats.amps,stats.taum(:,cellind),'.-'); xlabel('I'); ylabel('taum'); title(sprintf('median=%g',stats.taum_median(cellind)));

%% Task model: net_drives()

spec=spec0; 
pops=[1 3];
spec.nodes=spec.nodes(pops);
spec.connections=spec.connections(pops,pops);

% parms.minamp=-2;
% parms.maxamp=2; 
% parms.damp=1;
% parms.pulsetime=1000;
% parms.baseline=100;
% parms.override={spec.nodes(1).label,'parameters','Cm',1};
% parms.targets={'E1','Ifast1'};

[allMUA,netstats] = net_drives(spec,'targets',{'E1','Ifast1'},'plot_flag',1,...
  'minamp',-2,'maxamp',2,'damp',1,'pulsetime',1000,'baseline',100);

%% Task model: net_tuning()

% parms.minamp=1;
% parms.maxamp=2; 
% parms.damp=1;
% parms.minfreq=10;
% parms.maxfreq=60;
% parms.dfreq=2.5;
% parms.pulsetime=1000;
% parms.baseline=100;

override={spec.nodes(1).label,'parameters','Cm',1};

[allMUA2,netstats2] = net_tuning(spec,'targets',{'E1'},...
  'minamp',1,'maxamp',2,'damp',1,'minfreq',10,'maxfreq',60,'dfreq',2.5,...
  'pulsetime',1000,'baseline',100,'plot_flag',1,'override',override);

%%

cd /home/jason/models/dnsim/aro/base
ruletask_ResponseLR_CB_competition;
% remove inputs and noise sources
exclude={'InputGenerator2','randn'};
for i=1:length(spec.nodes)
  spec.nodes(i).dynamics={'V''=(current)./Cm'};
  spec.nodes(i).mechanisms=setdiff(spec.nodes(i).mechanisms,exclude);
end
% store modified specification
spec0=spec; % full model

% sub-model to process
spec=spec0; 
pops=[1 3]; % select E (E1) and FS (Ifast1) from superficial layers
pops=[2 4]; % select E (E2) and FS (Ifast2) from deep layers
spec.nodes=spec.nodes(pops);
spec.connections=spec.connections(pops,pops);

% biophysical search space
biophys_node = spec.nodes(1).label;%'E1';
biophys_param = 'Cm';
biophys_values = [1 2];

tstart=tic;
% loop over biophysical values
for biophys_index=1:length(biophys_values)

  % vary:
  override={biophys_node,'parameters',biophys_param,biophys_values(biophys_index)};

  % cell pulses
  [allX,cellstats] = cell_pulses(spec,'minamp',-1,'maxamp',1,'damp',.1,...
    'pulsetime',1000,'baseline',500,'plot_flag',1,'override',override);

  % network drives
  [allMUA1,netstats1] = net_drives(spec,'minamp',1,'maxamp',2,'damp',1,...
    'pulsetime',1000,'baseline',500,'plot_flag',1,'override',override,...
    'targets',{'E1'});

  % network tuning
  [allMUA2,netstats2] = net_tuning(spec,'minamp',1,'maxamp',2,'damp',1,...
    'pulsetime',1000,'baseline',500,'plot_flag',1,'override',override,...
    'targets',{'E1'},'minfreq',10,'maxfreq',60,'dfreq',2.5);

  stats.(biophys_param)(biophys_index).cell_pulses=cellstats;
  stats.(biophys_param)(biophys_index).net_drives=netstats1;
  stats.(biophys_param)(biophys_index).net_tuning=netstats2;
  data.(biophys_param)(biophys_index).cell_pulses=allX;
  data.(biophys_param)(biophys_index).net_drives=allMUA1;
  data.(biophys_param)(biophys_index).net_tuning=allMUA2;
end
toc(tstart)

% collect metrics=f(input|Cm)
cellind=1;
taum=arrayfun(@(x)x.cell_pulses.taum_median(cellind),stats.(biophys_param));
fMUA=arrayfun(@(x)x.net_drives.OscFreq(:,1,cellind),stats.(biophys_param),'uni',0);
fc=arrayfun(@(x)x.net_tuning.filter_fc(:,cellind),stats.(biophys_param),'uni',0);
amps=stats.(biophys_param)(biophys_index).net_tuning.amps;

fMUA=[fMUA{:}]; % amps x Cm
fc=[fc{:}];     % amps x Cm
% [taum] = 1 x Cm

figure
subplot(3,1,1); 
plot(biophys_values,taum,'.-'); 
xlabel(sprintf('%s.%s',biophys_node,biophys_param));
ylabel('taum [ms]'); title('cell pulses');
subplot(3,1,2); 
plot(biophys_values,fMUA,'.-'); 
xlabel(sprintf('%s.%s',biophys_node,biophys_param));
ylabel('fMUA [Hz]'); title('net drives');
legend(cellfun(@(x)sprintf('%s.amp=%g',biophys_node,x),num2cell(amps),'uni',0));%,'location','eastoutside');
subplot(3,1,3); 
plot(biophys_values,fc,'.-'); 
xlabel(sprintf('%s.%s',biophys_node,biophys_param));
ylabel('fc [Hz]'); title('net tuning');
legend(cellfun(@(x)sprintf('%s.amp=%g',biophys_node,x),num2cell(amps),'uni',0));%,'location','eastoutside');


stats_sup = stats;

% next: vary tauI and look for relationship (fMUA,fc) ~ (taum,tauI)

% --------------------------------------------
% todo: (after finishing the processing section of this function)
% 1. create new func from this func to sweep sin(2*pi*finp) | amp. return: fc, BW
% 2. prepare script that sweeps Cm and calcs Reff, fc
%     - plot (Cm,Reff) (should be linear), and (Reff,fc) (probably monotonic)
% 3. apply observations to competition and creating new funcs to
% characterize that... consider the implications for heterogeneity

% 2. net_tuning()
%   - add fc and BW calculation to net_tuning() (see code in hetero analysis)
% 	- then use sPING (varying Cm or gleak) to look for relationship b/w cell tau and net fc

