%% Search PFC cell model param space for fits to tallie's data (07-Dec-2014)
clear all  % need this to clear persistent variable in holdstim()
if ~any(strfind(path,'/home/jason/models/dnsim/Tallie-ACd-cells/dependencies'))
  addpath(genpath('/home/jason/models/dnsim/Tallie-ACd-cells/dependencies'));
end
%cd /home/jason/models/dnsim/Tallie-ACd-cells/distributions/base;
cd(fileparts(mfilename('fullpath')));

% simulation parameters
dsfact=10; solver='euler'; dt=.01; usecoder=1;
% tonic input
onset=50; tonicstim=0; 

% step protocol
bltime=500; nsections=4; nsteps=1; isi=1000; steptime=400; tonictime=15000; % dsfact=10
stepstim=7; ramprate=4; dt=.01; membranearea=3000; tspan=[0 20000]; 

bltime=500; nsections=4; nsteps=1; isi=1000; steptime=400; tonictime=25000; % dsfact=10
stepstim=7; ramprate=.75; dt=.01; membranearea=3000; tspan=[0 30000]; 
  
bltime=500; nsections=4; nsteps=1; isi=1000; steptime=400; tonictime=85000; % dsfact=10
stepstim=7; ramprate=2; dt=.01; membranearea=3000; tspan=[0 85000]; 

bltime=500; nsections=4; nsteps=1; isi=1000; steptime=400; tonictime=100000; % dsfact=10
stepstim=7; ramprate=1.5; dt=.01; membranearea=3000; tspan=[0 100000]; 
% 
% bltime=500; nsections=4; nsteps=1; isi=2000; steptime=1000; tonictime=100000; % dsfact=10
% stepstim=7; ramprate=.6; dt=.01; membranearea=3000; tspan=[0 150000]; 
% 
% bltime=500; nsections=4; nsteps=1; isi=1000; steptime=400; tonictime=100000; % dsfact=10
% stepstim=7; ramprate=.6; dt=.01; membranearea=3000; tspan=[0 150000]; 

% lfac=3;
% bltime=500; nsections=5; nsteps=1; isi=lfac*1000; steptime=lfac*400;  % dsfact=10
% stepstim=7; ramprate=.6; dt=.01; membranearea=3000; tspan=[0 140000]; tonictime=tspan(2); 

lfac=1; depol=1; stepreset=0;
bltime=500; nsections=4; nsteps=1; isi=lfac*1000; steptime=lfac*400;  % dsfact=10
stepstim=7; ramprate=.75; dt=.01; membranearea=3000; tspan=[0 160000]; tonictime=tspan(2); 

% biophysical parameters
ECa=126.1; EK=-80; ENa=55; Cm=1; ncells=1;

gNap = .01;%.015;
wbgNa = 0;%25; % 25
gnaf = 25;%75;%100;%100;%25; %100
gkdr = 6;%6;%8;
gAHP = .01;%.025;%.05;
gks = 0;%.285;%.25
gM = .01;
taurca = 8*28.5714; % 8*28.5714
gcan = .001;%.0056;
gkca = .0;%.25; %0;%.6;
taumin = 0;
gh = .01;%.005;%.002;
eh = -10;
gpas = .04; % vary this to get Vrest distribution
epas = -66;%-85;%-70;%-69;%-66;
noise=.2;

% tspan=[0 40000];

% Model specification (base model)
spec=[];
spec.nodes(1).label = 'PYs';
spec.nodes(1).multiplicity = ncells;
spec.nodes(1).dynamics = {'v''=current/c'};
spec.nodes(1).mechanisms = {'naf','E_Nap','wbNa','cadyn','can','kca','kdr','M','AHP','h','pas','stim','randn','iStepProtocol'};
spec.nodes(1).parameters = {...
  'c',Cm,'v_IC',-65,'stim',tonicstim,'noise',noise,'onset',onset,'depol',depol,'stepreset',stepreset,...
  'isi',isi,'nsteps',nsteps,'steptime',steptime,'stepsize',stepstim,'nsections',nsections,'tonictime',tonictime,'membranearea',membranearea,'bltime',bltime,'ramprate',ramprate,...
  'gnaf',gnaf,'taurca',taurca,'gcan',gcan,'gNap',gNap,'wbgNa',wbgNa,...
  'gkca',gkca,'taumin',taumin,'gkdr',gkdr,'gks',gks,'gAHP',gAHP,'gM',gM,...
  'gh',gh,'gpas',gpas,'epas',epas,'eh',eh,'ek',EK,'eca',ECa,'ena',ENa};

data = runsim(spec,'timelimits',tspan,'dt',dt,'dsfact',dsfact,'debug',1,'SOLVER',solver,'verbose',0,'coder',usecoder);
plotv(data,spec,'varlabel','v');

% My analysis
AnalyzeStepProtocol;
  % ephys_labels,ephys_metrics
  % [Idep;fdep]
  % [Idep;fdep0]
  % [Idep;fdepSS]
  % [I;v]

% Tallie analysis
% disp('tallie''s analysis procedure')
% tmp.sim_data=data; tmp.spec=buildmodel(spec,'verbose',0);
% proc=getcharacteristics(tmp,'sim',[],[],-1);
% plotcharacteristics(proc,tmp.spec);
% [simphysiol,labels] = getcellparams(proc{1});
% for i=1:length(labels), fprintf('%s: %g\n',labels{i},simphysiol(i)); end

%{
Approx mean experimental values:
AHPtime: 35
AHPsize: n/a
spikewidth: 1.5
spikefreq: 30
isidiff23: n/a
threshfreq: 1.5
Vrest: -70
%}

return


scope = {'PYs' , 'PYs' ,'PYs' ,'PYs' , 'PYs' , 'PYs' ,'PYs' };
variable = {'epas'                    , 'gkdr' ,'gnaf'          ,'gAHP'          ,'gks'        ,'gkca'             ,'gh' };
values = {'[-85 -80 -75 -70 -65 -60]','[4 6 8]','[25 50 75 100]','[0 .025 .05 .075 .1]' ,'[0 .25 .5 .75 1]','[0 .25 .5 .75 1]','[0 .005 .010 .015 .020]'};

%{

save_flag=1; cluster_flag=1;

vals='[.0001 .001 .01]';
scope    = {'PYs' , 'PYs'  ,'PYs'  ,'PYs' , 'PYs'  , 'PYs'        };
variable = {'gcat', 'gcan' ,'gkca' ,'gM'  , 'gAHP' , 'gpas'       };
values   = {vals  ,  vals  , vals  , vals ,  vals  ,'[.04 .1]'};

[~,~,outdir]=simstudy(spec,scope,variable,values,'coder',1,...
  'dt',dt,'SOLVER','euler','rootdir',pwd,'timelimits',tspan,'dsfact',10,'sim_cluster_flag',cluster_flag,...
  'savedata_flag',(save_flag||cluster_flag),'savepopavg_flag',0,'savespikes_flag',0,'saveplot_flag',(save_flag||cluster_flag),...
  'plotvars_flag',1,'plotrates_flag',0,'plotpower_flag',0,'addpath','~/dnsim');
outdir

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ON Local Machine 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shell:
% cd /home/jason/models/dnsim/Tallie-ACd-cells/distributions/base
% scp sherfey@scc2.bu.edu:/projectnb/crc-nak/sherfey/projects/heterogeneity/cells/PYs-gh-gks-gkca-gcat-gcan-gAHP-epas_simphysiol.mat .

% ----------------------------------------------------
res=load('PYs-gh-gks-gkca-gcat-gcan-gAHP-epas_simphysiol.mat');
vars={'gh','gks','gkca','gcat','gcan','gAHP','epas'};

%res=load('PYs-epas-gpas-c-noise-gkca-gAHP-gh_simphysiol.mat');
%res=load('PYs-gpas-gAHP-gM-gkca-gcan-gcat_simphysiol.mat'); % load results copied from scc
  % tallie/PYs-gpas-gAHP-gM-gkca-gcan-gcat/20141208
  spaces=res.spaces;
  %vars=res.vars;
  simphysiol=res.simphysiol_use;
  labels=res.labels;
% ----------------------------------------------------
  
vars{1},spaces{1}
simphysiol,labels
% ----------------------------------------------------

% experimental distributions of physiological measures
celldist = load('/home/jason/models/dnsim/Tallie-ACd-cells/distributions/tallie_distributions_1.mat');
cellids = celldist.ic;
paramlabels = celldist.head;
paramvalues = celldist.D;    % [cells x params]

% compare distributions of experimental paramvalues & simulation simphysiol
expparmi=[6 8 9 10 11]; % [5 6 8 9 10 11]
simparmi=[1 3 6 7 4];   % [2 1 3 6 7  4]
expparms=paramvalues(:,expparmi);
explabel=paramlabels(expparmi);
simparms=simphysiol(simparmi,:)';
simlabel=labels(simparmi);
numparms=length(simparmi);

% plot and compare exp vs sim distributions
figure('position',[100 180 1780 730]);
nbin=15; nr=2; nc=numparms; 
for i=1:numparms
  subplot(nr,nc,i);
  hist(expparms(:,i),nbin);
  xlabel(explabel{i}); title('exp');
  xlims=xlim;
  subplot(nr,nc,i+nc);
  ind=(simparms(:,i)<max(expparms(:,i))) & (simparms(:,i)~=0);
  hist(simparms(ind,i),nbin);
  xlabel(simlabel{i}); title('all sims');
  xlim(xlims);
end


% plot each unique variable in search space separately
numspaces=cellfun(@str2num,spaces);
uvars=unique(vars,'stable');
Rvals=nan(length(uvars),numparms); pvals=Rvals;
for j=1:length(uvars) % loop over model parameters in search space
  %sel=find(strcmp(uvars{j},vars));
  figure('position',[100 180 1780 730]);
  nbin=15; nr=3; nc=numparms; 
  for i=1:numparms % loop over physiological measures of cell responses
    sel=1:size(res.spaces,1);
    subplot(nr,nc,i);
    tmp=expparms(:,i); xmin=min(tmp(:)); xmax=max(tmp(:));
    hist(expparms(:,i),linspace(xmin,xmax,nbin));%nbin);
    xlabel(explabel{i}); title('exp');
    xlim([xmin xmax]); %xlims=xlim;
    subplot(nr,nc,i+nc); 
    tag='';
    % %
    if 1
      %ind=simparms(sel,i)<max(expparms(:,i));
      ind=(simparms(sel,i)<max(expparms(:,i))) & (simparms(sel,i)~=0);
      sel=sel(ind); tag='_inrange';
    end
    % %
    hist(simparms(sel,i),linspace(xmin,xmax,nbin));%nbin);
    xlabel(simlabel{i}); title(sprintf('all %s sims',uvars{j}));
    %axis tight;
    xlim([xmin xmax]); %xlim(xlims);
    xvals=numspaces(sel,j); %xvals=numspaces(sel,1);
    yvals=simparms(sel,i);
    tmpx=xvals(~isnan(yvals));
    tmpy=yvals(~isnan(yvals));
    if ~isempty(tmpy)
      [R,p]=corr(tmpx,tmpy);
    else
      R=0; p=inf;
    end
    N=length(find(~isnan(yvals)));
    subplot(nr,nc,i+2*nc);
    scatter(xvals,yvals); lsline
    %scatter(log10(xvals),yvals); lsline
    %scatter(log10(xvals),log10(yvals)); lsline
    xlabel(uvars{j}); ylabel(simlabel{i});
    title(sprintf('R=%g, p=%2.2g, N=%g',R,p,N));
    Rvals(j,i)=R; % model param x physiol measures
    pvals(j,i)=p;
  end  
  print(sprintf('PYs-gh-gks-gkca-gcat-gcan-gAHP-epas_%s-dep_simphysiol%s.jpg',uvars{j},tag),'-djpeg'); close  
  %print(sprintf('PYs-c-epas-eh-gh-gkca-gcan-taurca-gks-gAHP_%s-dep_simphysiol%s.jpg',uvars{j},tag),'-djpeg'); close  
  %print(sprintf('PYs-epas-gpas-c-noise-gkca-gAHP-gh_%s-dep_simphysiol%s.jpg',uvars{j},tag),'-djpeg'); close  
  %print(sprintf('PYs-noise-gcan-gcat-gkca-gAHP-gh_%s-dep_simphysiol.jpg',uvars{j}),'-djpeg'); close  
  %print(sprintf('mechanisms-log-plus-%s_simphysiol.jpg',uvars{j}),'-djpeg'); close
  %print(sprintf('mechanisms-log-plus-%s_simphysiol_log-scale.jpg',uvars{j}),'-djpeg'); close
  %print(sprintf('mechanisms-log-plus-%s_simphysiol_log-log-scale.jpg',uvars{j}),'-djpeg'); close
  %print(sprintf('mechanisms-log-plus-%s_simphysiol_axis-tight.jpg',uvars{j}),'-djpeg'); close
end

% plot correlation matrix
figure('position',[110 440 1800 530]);
subplot(1,2,1);
imagesc(Rvals); caxis([-1 1]); colorbar; title('all corrs');
set(gca,'ytick',1:length(uvars),'yticklabel',uvars);
set(gca,'xtick',1:length(simlabel),'xticklabel',simlabel);
subplot(1,2,2);
imagesc(Rvals.*(pvals<.001)); caxis([-1 1]); colorbar; title('R (p<.001)');
set(gca,'ytick',1:length(uvars),'yticklabel',uvars);
set(gca,'xtick',1:length(simlabel),'xticklabel',simlabel);
%print(sprintf('mechanisms-log-plus_simphysiol_correlations.jpg',uvars{j}),'-djpeg');

% boundaries on experimental distributions
expranges = [min(expparms,[],1);max(expparms,[],1)];

% find sims with phys measures within those ranges
% (look at simstudy params (from loadbatchmodels) producing simphysiol b/w ranges)
inrange = cellfun(@(i)simparms(:,i)>=expranges(1,i)&simparms(:,i)<=expranges(2,i),num2cell(1:size(expranges,2)),'uni',0);
inrange = cat(2,inrange{:});
sel = find(sum(inrange,2)==numparms);

% print search space info
fprintf('sim search space:\n');
for i=1:length(uscopes)
  fprintf('%s.%s [%s]\n',uscopes{i},uvars{i},uvals{i});
end
spaces{1}

% print info on search space params that produce physiol responses in exp range
disp('sim search space with responses in exp range:');
for i=1:numel(scopes{1})
  fprintf('(%s.%s) \t',scopes{1}{i},vars{1}{i});
end
fprintf('(%g sims)\n',size(spaces{1},1));
for i=1:numel(sel)
  for j=1:numel(scopes{1})
    fprintf('(%s) \t',spaces{1}{sel(i),j});
  end
  fprintf('(sim %g)\n',sel(i));
end

% plot and compare exp vs sim distributions in range
figure('position',[100 180 1780 730]);
nbin=15; nr=2; nc=numparms; 
for i=1:numparms
  subplot(nr,nc,i);
  hist(expparms(:,i),nbin);
  xlabel(explabel{i}); title('exp');
  xlims=xlim;
  subplot(nr,nc,i+nc);
  hist(simparms(sel,i),nbin);
  xlabel(simlabel{i}); title('sims (in range)');
  xlim(xlims);
end






