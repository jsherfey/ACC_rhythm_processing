addpath(genpath('/home/jason/models/dnsim/Tallie-ACd-cells/dependencies'));
cd /home/jason/models/dnsim/Tallie-ACd-cells/distributions/base;
% inputs
onset=50; tonicstim=1; % tonic
bltime=500; nsections=4; nsteps=1; isi=1000; steptime=400; stepstim=7; tonictime=5000; % step protocol
% biophysical parameters
ECa=126.1; EK=-80; ENa=55; Cm=1; noise=0; ncells=1;

gnaf = 100;
gnap = 0;%.0005;
taurca = 28.5714;
gcat = 0;%.001;
gcan = 0;%.0056;
gkca = 0;%.6;
taumin = 0;
gkdr = 6;%6;
gks = 0;%.285;
gM = 0;%.084;
gAHP = 0;%.054;
gh = .48;%.002;
gpas = .1;%.04;
epas = -66;
eh = -10;
membranearea = 1500; % vary this to see if the exp heterogeneity could be due to cell size
% -------------------------------------------
%stepstim=7; tonicstim=0;
stepstim=0; tonicstim=.02; noise=1; gh=1; eh=-10;
tspan=[0 1000]; dt=.05; solver='euler';
% -------------------------------------------
%{
theta STO (h-mediated damped oscillation perpetuated by noisy kicks): 
  stepstim=0; tonicstim=0; noise=.1; gh=1; eh=-10;
  where h.txt: rtau(v)=1.5*(exp(0.033*(v+75))./(0.011*(1+exp(0.083*(v+75)))));
  and (gpas=.1; gnap=gcat=gcan=gkca=gks=gM=gAHP=0)
%}
% Model specification (base model)
spec=[];
spec.nodes(1).label = 'E';
spec.nodes(1).multiplicity = ncells;
spec.nodes(1).dynamics = {'v''=current/c'};
spec.nodes(1).mechanisms = {'naf','nap','cadyn','cat','can','kca','kdr','iks','M','AHP','h','pas','stim','randn','iStepProtocol'};
spec.nodes(1).parameters = {...
  'c',Cm,'v_IC',-65,'stim',tonicstim,'noise',noise,'onset',onset,...
  'isi',isi,'nsteps',nsteps,'steptime',steptime,'stepsize',stepstim,'nsections',nsections,'tonictime',tonictime,'membranearea',membranearea,'bltime',bltime,...
  'gnaf',gnaf,'gnap',gnap,'taurca',taurca,'gcat',gcat,'gcan',gcan,...
  'gkca',gkca,'taumin',taumin,'gkdr',gkdr,'gks',gks,'gM',gM,'gAHP',gAHP,...
  'gh',gh,'gpas',gpas,'epas',epas,'eh',eh,'ek',EK,'eca',ECa,'ena',ENa};

% Process simulation
data = runsim(spec,'timelimits',tspan,'dt',dt,'dsfact',10,'SOLVER',solver);
plotv(data,spec,'varlabel','v');% ylim([-70 -50]);
plotpow(data,spec,'NFFT',4096,'varlabel','v'); % 16384, 10000

t=data.epochs.time*1000;
V=data.epochs.data(end,:);
[Psd, freq, PsdManSmooth, freqbin] = powerspectrum(t,V(5:end));
[psdmax, jnat]=max(Psd); fnat = freq(jnat); [fnat psdmax]
vline(fnat,'k'); title(sprintf('f=%g Hz',round(fnat*10)/10));
    