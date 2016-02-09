% Model description: ACC network with 1 E population and 2 I populations
% varying in their time constants of inhibition (1 fast and 1 slow).
% Two oscillatory inputs are delivered to mutually exclusive subsets of
% E-cells; responses are assessed in terms of ROI coactivity, synchrony, 
% and spectral content. The question being asked is: how does biophysical
% heterogeneity and the number of intrinsic network frequencies (i.e., # of 
% inhibitory time constants) affect the network response to multiple inputs 
% at different frequencies (representing different dimensions of a multimodal 
% signal input along parallel pathways or simultaneous inputs from any 
% arbitrary source regions with different natural or relayed frequencies)? 
% Furthermore, we are interested in the additional effects of weak 
% recurrent coupling and different values of inhibitory time constants 
% (i.e., producing different intrinsic/natural network frequencies).
% Any effects on network filter characteristics (pass-band and bandwidth)
% would also be interesting.

% Note: the only direct sensory inputs to ACC in macaque come from auditory
% cortex. ACC gets inputs from all regions of PFC, including OFC, which
% itself receives direct multimodal inputs from various sensory regions,
% and from LPFC, which receives segregated auditory and visual inputs as 
% well as motor inputs from PMdr. Finally, ACC receives inputs from 
% Hippocampus (whereas LPFC does not). (source: Barbas lab).

% Future models that will build on this:
% This modeling projects takes the perspective that the 2 interneuron
% population version is a "deep layer" ACC model, and a 1 population
% version is a "superficial layer" ACC model, roughly motivated by
% distributions of CB+ and PV+ interneurons in ACC, with both being present
% in deep layers and mainly CB+ being present in superficial layers. A
% future "multi-layer" version of this model will include two modules, each
% like the one in this script, coupled by interlaminar connections. Beyond
% that multilayer ACC model, a yet larger model will be constructed to
% incorporate connections with DLPFC (ACC L5/6 <-> DLPFC L2/3) and thalamus
% (VA and MD; core and matrix). This model is considered an "ACC" model 
% because it incorporates biophysical heterogeneity that produces 
% heterogeneous physiological properties similar to those observed
% in rat dorsal AC (paper submitted; data from LeBeau lab), which exceeds 
% heterogeneity observed in other prefrontal regions (e.g., PrL and IL).

project_dir='~/projects/ACC_rhythm_processing';
if isempty(strfind(path,[project_dir, pathsep]));
  addpath(genpath(project_dir));  
end
tstart=tic;

%% SIMULATOR CONTROLS (what to do)
tspan=[0 1000];  % [ms], time limits on numerical integration
dt=.01;         % [ms], time step for numerical integration
solver='rk2';   % {options; 'rk4', 'rk2', 'rk1'}, solver to use
compile_flag=1; % whether to compile the simulation
verbose_flag=1; % whether to display log information
solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag};
% -------------------------------------------------------------------------
%% Population definitions

% population sizes
nE=8; % # of E-cells per population
nI=2; % # of I-cells per population
prefix=sprintf('ACC-1Layer-2inputs_%gE%gIs%gIf',nE,nI,nI);
% -------------------------------------------------------------------------
% general input parameters (used for all inputs, whether experiment or other):
inpE=7; % E input amplitude [uA/cm2]
inpI=0;  % I input amplitude [uA/cm2]
KE1=1;KE2=1;KI1=0;KI2=0;
% -------------------------------------------------------------------------
% Connectivity kernels and parameters:

% time constants [ms]
tauE=2;   % E
tauIf=5; % If
tauIs=15; % Is

% max conductances [uS/cm2] (homogeneous or lognormal)
gEE=0;  % E->E
gEI=.2; % I->E
gIE=.2; % E->I
gII=0;  % I->I

% kernels
netconEE=zeros(nE,nE); % E->E, N_pre x N_post
netconEI=ones(nI,nE);  % I->E, N_pre x N_post
netconIE=ones(nE,nI);  % E->I, N_pre x N_post
netconII=zeros(nI,nI); % I->I, N_pre x N_post
% -------------------------------------------------------------------------
% Heterogeneous biophysics (distribution options: uniform or lognormal): E,I:{gi,wi}(currents)

% type of heterogeneity
heterogeneity_degree=0; % options: {0,1,2}
distribution='uniform'; % options: {'uniform','normal','lognormal'}

% degree of heterogeneity
switch heterogeneity_degree
  case 0          % homogeneous population
    sigma=0;
  case 1          % 10% spread (standard deviation)
    sigma=.1*mu;
  case 2          % 25% spread (standard deviation)
    sigma=.25*mu;
end

% distribution type
switch distribution
  case 'uniform'
    paramdist=@(mu,sigma,num_vals)unifrnd(mu-sigma,mu+sigma,[1 num_vals]);
  case 'normal'
    paramdist=@(mu,sigma,num_vals)normrnd(mu,sigma,[1 num_vals]);
  case 'lognormal'
    MU=log(mu^2/sqrt(sigma^2+mu^2));
    SIGMA=sqrt(log(sigma^2/mu^2+1));
    paramdist=@(mu,sigma,num_vals)lognrnd(MU,SIGMA,[1 num_vals]);
end

% draw parameter values from desired distribution
mu=1; % mean parameter value (mean gmax)
param_values=paramdist(mu,sigma,nE);

% gNa
% gK
% gM
% gH
% gNap
% Epas

% -------------------------------------------------------------------------
% Homogeneous biophysics
% Cm
% Rpas
% ENa
% EK
% Eh
% -------------------------------------------------------------------------
% Analysis parameters
bin_size=30;  % [ms], for spike rate calculation
bin_shift=10; % [ms], for spike rate calculation
% -------------------------------------------------------------------------
% Plotting parameters
plot_var='v'; % plot voltages [mV] for all populations
% -------------------------------------------------------------------------
% Network Specification
s=[];
s.pops(1)    =ACC_Ecell_specification;
s.pops(end).name='E';
s.pops(end).size=nE;
s.pops(end).parameters={'input',inpE,'K1',KE1,'K2',KE2};
s.pops(end+1)=ACC_Icell_specification;
s.pops(end).name='If';
s.pops(end).size=nI;
s.pops(end).parameters={'input',inpI,'K1',KI1,'K2',KI2};
s.pops(end+1)=ACC_Icell_specification;
s.pops(end).name='Is';
s.pops(end).size=nI;
s.pops(end).parameters={'input',inpI,'K1',KI1,'K2',KI2};
s.cons=[];
s.cons(1).direction    ='E->If';
s.cons(end).mechanism_list='iAMPA';
s.cons(end).parameters={'gSYN',gIE,'tauI',tauE,'netcon',netconIE};
s.cons(end+1).direction='E->Is';
s.cons(end).mechanism_list='iAMPA';
s.cons(end).parameters={'gSYN',gIE,'tauI',tauE,'netcon',netconIE};
s.cons(end+1).direction='E<-If';
s.cons(end).mechanism_list='iGABAa';
s.cons(end).parameters={'gSYN',gEI,'tauI',tauIf,'netcon',netconEI};
s.cons(end+1).direction='E<-Is';
s.cons(end).mechanism_list='iGABAa';
s.cons(end).parameters={'gSYN',gEI,'tauI',tauIs,'netcon',netconEI};

% set amplitude of network input
gINPUT=5;
plot_flag=1;
% todo: base on drive to achieve natural frequency near a target value

if 0
  % determine precise natural frequency of network with tonic injection
  s_=ApplyModifications(s,{'E','input',gINPUT;'Is','input',0;'If','input',0});
  data=SimulateModel(s_,solver_options{:}); % tonic drive
  data=CalcPower(data);
  f0=round(data.E_v_Power_MUA.PeakFreq);

  % prepare model for experiment (turn off inputs)
  s_=ApplyModifications(s,{'E','input',0;'Is','input',0;'If','input',0});

  % set parameters of multimodal rhythm processing experiment
  f1=f0+[-5 5];%[-5 -2.5 0 2.5 5];
  f2=f0+[-5 5];%[-5 -2.5 0 2.5 5];
  target='E';           % options: 'E','If','Is',{'E','If'},{'E','Is'},{'E','If','Is'}
  input_type='poisson'; % options: 'poisson','sin','noisy_sin','rectified_sin'
  gINPUT=5; dc=0; ac=1; baseline=.1; tau=2;
  num_repetitions=1;

  data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  PlotData(data,'variable',plot_var,'plot_type','waveform');
  PlotData(data(1).model.fixed_variables.E_s1,'plot_type','waveform');
  PlotData(data(1).model.fixed_variables.E_s2,'plot_type','waveform');
  PlotData(data(1).model.fixed_variables.E_s1,'plot_type','power','xlim',[0 100]);
  PlotData(data(1).model.fixed_variables.E_s2,'plot_type','power','xlim',[0 100]);
end
% simulation time (simulated 1 sec):
% rk1: 7.5 sec
% rk2: 8.5 sec
% rk4: 11.5 sec

% Plot the Poisson inputs from @ProbeTwoRhythms:
T=data(1).time;
K1=data(1).model.fixed_variables.E_K1;
K2=data(1).model.fixed_variables.E_K2;
s1=data(1).model.fixed_variables.E_s1;
s2=data(1).model.fixed_variables.E_s2;
Ginput=gINPUT*(repmat(K1,length(T),1).*s1+repmat(K2,length(T),1).*s2);
Iinput=-Ginput.*data(1).E_v;
PlotData(Ginput,'plot_type','waveform');
PlotData(Iinput,'plot_type','waveform');


% ----------------------------------------------
% WORKSHEET (05-Feb-2016)

% Q1: is there a regime where one population suppresses the other?

% S01: turn off tonic inputs and vary I->E strength
experiment=@ProbeTwoRhythms;
modifications={'E','input',0;'Is','input',0;'If','input',0;
               'E->If','gSYN',.2; 'E->Is','gSYN',.2};
s_=ApplyModifications(s,modifications);
vary={'If->E','gSYN',[0 .2];'Is->E','gSYN',[0 .2]};
data=SimulateModel(s_,'vary',vary,'experiment',experiment,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram','xlim',[0 500])

% S02: turn off tonic inputs and vary E->I strength
experiment=@ProbeTwoRhythms;
modifications={'E','input',0;'Is','input',0;'If','input',0;
               'If->E','gSYN',.2; 'Is->E','gSYN',.2};
s_=ApplyModifications(s,modifications);
vary={'E->If','gSYN',[0 .2];'E->Is','gSYN',[0 .2]};
data=SimulateModel(s_,'vary',vary,'experiment',experiment,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram')
PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100])

% S03: vary (f1,f2) over f0+/-20Hz where f0=40Hz | gINPUT=5
f1=f0+(-20:5:20); f2=f1; target='E'; input_type='poisson';
gINPUT=5; dc=0; ac=1; baseline=.1; tau=2;
num_repetitions=1;
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram');
plottype='rastergram'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');
PlotData(data,'variable','E_v','plot_type','power','xlim',[0 100]);
plottype='power'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');
clear data;

% S04: turn off tonic inputs and vary gINPUT with f1=20-60Hz and f2=40Hz
f1=20:5:60; f2=40; gINPUT=1:9; dc=0; ac=1; baseline=.1; tau=2;
num_repetitions=1; target='E'; input_type='poisson';
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram');
plottype='rastergram'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');
PlotData(data,'variable','E_v','plot_type','power','xlim',[0 100]);
plottype='power'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');

% S05: manually explore effects of different strength inputs
compile_flag=0; solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag};
f1=55; f2=40; gINPUT=.05; dc=0; ac=1; baseline=.1; tau=2;
num_repetitions=1; target='E'; input_type='poisson';
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
PlotData(data,'variable','E_v','plot_type','power','xlim',[0 100]);
PlotData(data,'variable',plot_var,'plot_type','rastergram');
% observation: much weaker inputs are sufficient to drive E-cells (Na,K) w/ Poisson input

% S06: repeat S04 with weaker inputs
compile_flag=1; solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag};
f1=20:5:60; f2=40; gINPUT=.01:.01:.09; dc=0; ac=1; baseline=.1; tau=2;
num_repetitions=1; target='E'; input_type='poisson';
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram');
plottype='rastergram'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');
PlotData(data,'variable','E_v','plot_type','power','xlim',[0 100]);
plottype='power'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');
% observation:
  % gINPUT<.04: population driven at freq closer to natural freq dominates very significantly
  % gINPUT=.05: ~balanced
  % gINPUT .06-.09: faster has slight advantage
  % gINPUT >1: slower has slight advantage

% space of interest:
% f1=[20 40 60]
% gINPUT: .01-.08

% S07: same as S06 except focused on space of interest
% note: may need to do multiple realizations; start w/ 1 and compare to S06
gINPUT=[.01 .02 .08]; f1=[20 40 60]; f2=40;
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram');
plottype='rastergram'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');
PlotData(data,'variable','E_v','plot_type','power','xlim',[0 100]);
plottype='power'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');

gINPUT=.01:.01:.08; f1=[20 40 60]; f2=40; num_repetitions=3;
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
PlotData(data,'variable',plot_var,'plot_type','rastergram');
plottype='rastergram'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%gx_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),num_repetitions,plottype),'-djpeg');
PlotData(data,'variable','E_v','plot_type','power','xlim',[0 100]);
plottype='power'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%gx_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),num_repetitions,plottype),'-djpeg');

% S08: examine smaller differences f1-f2
gINPUT=[.01 .02 .08]; f1=35:2.5:45; f2=40;
% ...

% Test01: determine precise natural frequency of network with homogeneous poisson inputs
modifications={'E','input',0;'Is','input',0;'If','input',0};
s_=ApplyModifications(s,modifications);
target='E'; input_type='poisson';
gINPUT=.01:.01:.05; f=10:10:80; num_repetitions=3; dc=0; ac=1; baseline=.1; tau=2;
[data,stats]=ProbeResonance(s_,'f',f,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});

% TODO:
% 1. create @ProbeResonance (from ProbeTwoRhythms) to vary nonhomogeneous 
%    poisson frequency of a single input. should return data with fMUA and
%    area power (from CalcPower) as well as fc for each sweep of 
%    frequencies (based on population FR from CalcFR).
% 2. using @ProbeResonance: vary gINPUT (gINPUT) (w/num_repetitions=5):
%     - with f=0 (homogeneous poisson):          calc fMUA=f(gINPUT) [based on Pxx]
%     - with f=0-100Hz (nonhomogeneous poisson): calc fc=f(gINPUT)   [based on FR]
%     (average estimates over 5 realizations; for E, If, and Is)
%     (also plot Pxx with error shading over num_repetitions)
% 3. add cell ROIs to SelectData
% 4. using @ProbeTwoRhythms + SelectData & CalcFR: vary gINPUT (gINPUT), f1,
%    and f2 then generate plots revealing regimes of competition vs
%    cooperation between E1 and E2: examine |(r1/sum(I1))-(r2/sum(I2))|
%    and plots outlined on paper.

% - figure out how to pass parameters to experiment through SimulateModel;
%   could use to examine how tauI affects fMUA and fc | @ProbeResonance.

% after understanding what's going on in the simple homogeneous population
% of (Na,K) cells, study the effects of adding non-spiking currents, then
% heterogeneity of single parameters (note: will require adding a layer of
% model realizations on top of everything else). i.e., continue with steps
% outlined with Nancy.

% ------------------------------------
tstart=tic;

% turn off injected inputs
modifications={'E','input',0;'Is','input',0;'If','input',0};
s_=ApplyModifications(s,modifications);

target='E'; input_type='poisson';

% Experiment #1: homogeneous poisson probing natural frequency
gINPUT=0:.01:.1; f=0; num_repetitions=5; dc=0; ac=1; baseline=.1; tau=2;
%gINPUT=0:.01:.22; f=0; num_repetitions=5; dc=0; ac=1; baseline=.1; tau=2;
[~,stats_f0]=ProbeResonance(s_,'f',f,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
plottype='Exp1-ProbeResonance-HomogeneousPoisson'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f(1),f(end),plottype),'-djpeg');

f0=stats_f0.E_v_Power_MUA.repetition_sets.PeakFreq_mu;
f0e=stats_f0.E_v_Power_MUA.repetition_sets.PeakFreq_sd/sqrt(stats_f0.num_repetitions);
% figure; plot(gINPUT,f0,'o-'); xlabel('input strength'); ylabel('network natural frequency');

figure
var=stats_f0.MUA_variables{1};
x=stats_f0.(var).repetition_sets.parameters;
y=stats_f0.(var).repetition_sets.PeakFreq_mu;
e=stats_f0.(var).repetition_sets.PeakFreq_sd/sqrt(stats_f0.num_repetitions);
subplot(2,1,1); plot_CI(x,y,e,'b');
xlabel(strrep([stats_f0.(var).repetition_sets.varied{:}],'_','\_'))
ylabel('Peak Spectral Frequency [Hz]');
y=stats_f0.(var).repetition_sets.PeakArea_mu;
e=stats_f0.(var).repetition_sets.PeakArea_sd/sqrt(stats_f0.num_repetitions);
subplot(2,1,2); plot_CI(x,y,e,'b');
xlabel(strrep([stats_f0.(var).repetition_sets.varied{:}],'_','\_'))
ylabel('Area Power Around Peak');
plottype='Exp1-MUA-power'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f(1),f(end),plottype),'-djpeg');

% Experiment #2: nonhomogeneous poisson probing center frequency
gINPUT=.01:.01:.06; f=[20 40 60]; num_repetitions=4; dc=0; ac=1; baseline=.1; tau=2;
%gINPUT=.01:.01:.09; f=10:5:75; num_repetitions=1; dc=0; ac=1; baseline=.1; tau=2;
[~,stats_fc]=ProbeResonance(s_,'f',f,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,solver_options{:});
plottype='Exp2-ProbeResonance-RhythmicPoisson'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f(1),f(end),plottype),'-djpeg');

fc=stats_fc.E_v_FR.sweep_sets.repetition_sets.sweep_pop_FR_max_mu;
fce=stats_fc.E_v_FR.sweep_sets.repetition_sets.sweep_pop_FR_max_sd/sqrt(stats_fc.num_repetitions);
% figure; plot(gINPUT,fc,'o-'); xlabel('input strength'); ylabel('network center frequency');

figure
var=stats_fc.FR_variables{1};
x=stats_fc.(var).sweep_sets.repetition_sets.parameters;
y=stats_fc.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_mu;
e=stats_fc.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_sd/sqrt(stats_fc.num_repetitions);
plot_CI(x,y,e,'r');
xlabel(strrep(stats_fc.(var).sweep_sets.repetition_sets.varied{1},'_','\_'))
ylabel('Resonance Frequency [Hz]');
plottype='Exp2-Tuning-curve'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f(1),f(end),plottype),'-djpeg');

x=gINPUT;
y1=fc(1:length(gINPUT)); y1e=fce(1:length(gINPUT));
y2=f0(1:length(gINPUT)); y2e=f0e(1:length(gINPUT));
figure; 
subplot(2,1,1); plot_CI(x,y1,y1e,'r'); hold on; plot_CI(x,y2,y2e,'b'); 
xlabel('gINPUT'); ylabel('freq [Hz]'); legend('fc','fMUA');
subplot(2,1,2); plot(gINPUT,y1-y2,'k-o'); 
xlabel('gINPUT'); ylabel('fc-fMUA [Hz]'); ylim([0 40])
plottype='Exps1-2_fc-vs-fMUA'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f(1),f(end),plottype),'-djpeg');

% Experiment #3: two nonhomogeneous poisson inputs probing competition vs cooperation
gINPUT=[.01 .02 .03 .04 .05]; f1=[20 40 60]; f2=40; dc=0; ac=1; baseline=.1; tau=2;
num_repetitions=4; target='E'; input_type='poisson'; post_downsample_factor=2;
[data,stats]=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'num_repetitions',num_repetitions,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'baseline',baseline,'post_downsample_factor',post_downsample_factor,solver_options{:});

toc(tstart)

% Examine effects of varying (f1,f2) and gINPUT on E1 and E2
num_gINPUTs=length(gINPUT);
num_f1s=length(f1);
p_gINPUT=stats(1).E_v_FR.repetition_sets.parameters(1,:); % 1 x nsims(=num_gINPUTs*num_f1s)
p_f1=stats(1).E_v_FR.repetition_sets.parameters(2,:);   % 1 x nsims
E1_FR=mean(stats(1).E_v_FR.repetition_sets.pop_FR_mu,1); % 1 x nsims
E2_FR=mean(stats(2).E_v_FR.repetition_sets.pop_FR_mu,1); % 1 x nsims
E_fc=fc(1:num_gINPUTs); % 1 x num_gINPUTs
E_f0=f0(1:num_gINPUTs); % 1 x num_gINPUTs

% calculate mean input strengths (Imean)
data_gINPUT=[data.E_gINPUT];
data_f1=[data.E_f1];
E1_Imean=zeros(num_gINPUTs,num_f1s);
E2_Imean=zeros(num_gINPUTs,num_f1s);
for i=1:num_gINPUTs
  for j=1:num_f1s
    dat=data(data_gINPUT==gINPUT(i) & data_f1==f1(j));
    I1avg=0; I2avg=0;
    for k=1:length(dat)
      T=dat(k).time;
      K1=dat(k).model.fixed_variables.E_K1;
      K2=dat(k).model.fixed_variables.E_K2;
      s1=dat(k).model.fixed_variables.E_s1(1:post_downsample_factor:end,:);
      s2=dat(k).model.fixed_variables.E_s2(1:post_downsample_factor:end,:);
      G1=gINPUT(i)*(repmat(K1,length(T),1).*s1);
      G2=gINPUT(i)*(repmat(K2,length(T),1).*s2);
      Itmp=-(G1.*dat(k).E_v);                 % I1 current to all cells
      I1avg=I1avg+mean(mean(Itmp(:,K1==1)));  % average I1 current to E1 cells
      Itmp=-(G2.*dat(k).E_v);                 % I2 current to all cells
      I2avg=I2avg+mean(mean(Itmp(:,K2==1)));  % average I2 current to E2 cells
    end
    E1_Imean(i,j)=I1avg/length(dat);
    E2_Imean(i,j)=I2avg/length(dat);
  end
end
% plot the actual mean input delivered to E-cells
figure('position',[325 300 1320 550])
subplot(1,2,1); plot(f1,E1_Imean,'o-'); 
xlabel('input frequency to E1'); ylabel('mean input to E1'); 
legend(cellfun(@(x)['gINPUT=' num2str(x)],num2cell(gINPUT),'uni',0));
subplot(1,2,2); plot(f2*ones(1,length(f1)),E2_Imean,'o-'); 
xlabel('input frequency to E2'); ylabel('mean input to E2'); 
legend(cellfun(@(x)['gINPUT=' num2str(x)],num2cell(gINPUT),'uni',0));
plottype='mean-input-strength'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');

% Plot rastergrams for all experiment 3 sims
PlotData(data,'variable',plot_var,'plot_type','rastergram');
plottype='rastergram'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');

% for each gINPUT:
figure('position',[330 50 1200 910]); 
ncols=6; nrows=num_gINPUTs;
rlims=[min([E1_FR(:);E2_FR(:)]) max([E1_FR(:);E2_FR(:)])];
drlims=[min(E1_FR(:)-E2_FR(:)) max(E1_FR(:)-E2_FR(:))];
drlims=[-max(abs(drlims)) max(abs(drlims))];
for i=1:num_gINPUTs
  % this input strength
  gINPUT_=gINPUT(i);                      % scalar
  % indices to ProbeTwoRhythms sims with this input strength
  sel=find(p_gINPUT==gINPUT_);            % 1 x length(f1)
  % natural and resonance frequencies of network given this input strength
  fc_=E_fc(i);                        % scalar
  f0_=E_f0(i);                        % scalar
  % normalized input frequencies
  f1_dc=abs(f1-fc_);                        % 1 x length(f1)
  f2_dc=abs(f2-fc_);                        % scalar
  f2_dc=f2_dc*ones(size(f1_dc));            % 1 x length(f1)
  f1_d0=abs(f1-f0_);                        % 1 x length(f1)
  f2_d0=abs(f2-f0_);                        % scalar
  f2_d0=f2_d0*ones(size(f1_d0));            % 1 x length(f1)
  % collect FRs
  r1=E1_FR(:,sel);                    % 1 x length(f1)
  r2=E2_FR(:,sel);                    % 1 x length(f1)
  % collect mean input strengths
  I1=E1_Imean(i,:);                   % 1 x length(f1)
  I2=E2_Imean(i,:);                   % 1 x length(f1)
  % calculate competition index 
  r1_Inorm=r1./I1;                    % 1 x length(f1)
  r2_Inorm=r2./I2;                    % 1 x length(f1)
  competition=abs(r1_Inorm-r2_Inorm); % 1 x length(f1)
  % generate plots (one row per gINPUT)
  
  subplot(nrows,ncols,ncols*(i-1)+1) % (f,r)
  plot(f1,r1,'bo-',f1,r2,'ro-'); ylim(rlims);
  xlabel(sprintf('f1 (f2=%g)',f2)); ylabel('spike rate'); legend('E1','E2');
  xmin=min(xlim); xmax=max(xlim);
  ymin=min(ylim); ymax=max(ylim);
  text_xpos=xmin+.05*(xmax-xmin);
  text_ypos=ymin+.9*(ymax-ymin);
  text(text_xpos,text_ypos,sprintf('gINPUT=%g',gINPUT_));

  subplot(nrows,ncols,ncols*(i-1)+2) % (f1-f2,r1-r2)
  plot(f1-f2,r1-r2,'ko-'); ylim(drlims);
  line(xlim,[0 0],'linestyle','--','color','g'); 
  line([0 0],ylim,'linestyle','--','color','g');
  xlabel('f1-f2'); ylabel('rE1-rE2');

  subplot(nrows,ncols,ncols*(i-1)+3) % (dist(f1,f0)-dist(f2,f0),r1-r2)
  plot(f2_d0-f1_d0,r1-r2,'ko-'); ylim(drlims);
  line(xlim,[0 0],'linestyle','--','color','g'); 
  line([0 0],ylim,'linestyle','--','color','g');
  xlabel('|f2-fMUA|-|f1-fMUA|'); ylabel('rE1-rE2');
  
  subplot(nrows,ncols,ncols*(i-1)+4) % (dist(f1,fc)-dist(f2,fc),r1-r2)
  plot(f2_dc-f1_dc,r1-r2,'ko-'); ylim(drlims);
  line(xlim,[0 0],'linestyle','--','color','g'); 
  line([0 0],ylim,'linestyle','--','color','g');
  xlabel('|f2-fc|-|f1-fc|'); ylabel('rE1-rE2');
  
  subplot(nrows,ncols,ncols*(i-1)+5) % (f1-f2,r1-r2)
  plot(f1,competition,'ko-'); 
  xlim([min(f1) max(f1)]); ylim([0 10]);
  line(xlim,[0 0],'linestyle','--','color','g'); 
  line([0 0],ylim,'linestyle','--','color','g');
  xlabel(sprintf('f1 (f2=%g)',f2)); ylabel('competition');  
  
  subplot(nrows,ncols,ncols*(i-1)+6) % (I,r)
  plot(I1,r1,'bo-',I2,r2,'ro-'); ylim(rlims);
  xlabel('mean input'); ylabel('spike rate'); %legend('E1','E2');
  
end
plottype='Exps1-3_E1-vs-E2-summary'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');


% ------------------------------------

toc(tstart);
return

% todo: rename plot_CI to Plot_CI (update in ProbeResonance)

PlotData(data,'variable',plot_var,'plot_type','rastergram');
PlotData(data,'variable','v','plot_type','power','xlim',[0 100]);
PlotData(SelectData(data,'varied',{'E_realization',1}),'variable',plot_var,'plot_type','rastergram');


% -------------------------------------------------------------
% todo:
% - proceed with sequence outlined during meeting with Nancy
% -------------------------------------------------------------

%{
% SIN
input_type='sin'; % options: 'sin','noisy_sin','rectified_sin'
% run multimodal rhythm processing experiment
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});
if plot_flag
  % plot results over (f1,f2):
  % plot voltage-derived rastergram
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  plottype='rastergram'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  % plot voltage power spectrum
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  plottype='power'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  % plot voltage waveforms
  PlotData(data,'variable',plot_var,'plot_type','waveform');
  plottype='waveform'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  figure; plot(data(1).time,data(1).model.fixed_variables.E_I1(:,1)); ylabel('I1');
  plottype='input'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
end
data_sin=data;

% RECTIFIED_SIN
input_type='rectified_sin'; % options: 'sin','noisy_sin','rectified_sin'
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});
if plot_flag
  % plot results over (f1,f2):
  % plot voltage-derived rastergram
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  plottype='rastergram'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  % plot voltage power spectrum
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  plottype='power'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  % plot voltage waveforms
  PlotData(data,'variable',plot_var,'plot_type','waveform');
  plottype='waveform'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  figure; plot(data(1).time,data(1).model.fixed_variables.E_I1(:,1)); ylabel('I1');
  plottype='input'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
end
data_rectsin=data;

% NOISY RECTIFIED_SIN
input_type='noisy_rectified_sin'; % options: 'sin','noisy_sin','rectified_sin'
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});
if plot_flag
  % plot results over (f1,f2):
  % plot voltage-derived rastergram
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  plottype='rastergram'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  % plot voltage power spectrum
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  plottype='power'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  % plot voltage waveforms
  PlotData(data,'variable',plot_var,'plot_type','waveform');
  plottype='waveform'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
  figure; plot(data(1).time,data(1).model.fixed_variables.E_I1(:,1)); ylabel('I1');
  plottype='input'; print(gcf,sprintf('%s-injection-amp%g_%g-%gHz_%s.jpg',input_type,gINPUT,f1(1),f1(end),plottype),'-djpeg');
end

% Poisson input:
% homogeneous:
if 0
  input_type='poisson'; % options: 'poisson','sin','noisy_sin','rectified_sin'
  f1=0; f2=0; dc=0; ac=5; tau=.002; gINPUT=5;
  target='E';       % options: 'E','If','Is',{'E','If'},{'E','Is'},{'E','If','Is'}
  data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'num_repetitions',num_repetitions,solver_options{:});
  PlotData(data,'variable',plot_var,'plot_type','waveform');
  PlotData(data,'variable',plot_var,'plot_type','power');
  data=CalcPower(data);
  f0=round(data.E_v_Power_MUA.PeakFreq)
  figure; plot(data.time,data.model.fixed_variables.E_s1(:,1))
end

% nonhomogeneous:
if 0
  f1=f0+[-5 5];%[-5 0 5];%(-1:1);
  f2=f0+[-5 5];
  num_repetitions=1;
  input_type='poisson'; % options: 'poisson','sin','noisy_sin','rectified_sin'
  target='E';       % options: 'E','If','Is',{'E','If'},{'E','Is'},{'E','If','Is'}
  dc=0; ac=5; tau=2/1000; gINPUT=5;
  s_=ApplyModifications(s,{'E','input',0;'Is','input',0;'If','input',0});
  data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'num_repetitions',num_repetitions,solver_options{:});
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  PlotData(data,'variable',plot_var,'plot_type','waveform');
  figure; plot(data(1).time,data(1).E_v)
end

% POISSON: AC only
% f1=f0+[-5 -2.5 0 2.5 5];%[-5 0 5];%(-1:1);
% f2=f0+[-5 -2.5 0 2.5 5];
num_repetitions=1;
input_type='poisson'; % options: 'poisson','sin','noisy_sin','rectified_sin'
dc=0; ac=1; tau=2/1000; gINPUT=5;
s_=ApplyModifications(s,{'E','input',0;'Is','input',0;'If','input',0});
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'num_repetitions',num_repetitions,solver_options{:});
if plot_flag
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  plottype='rastergram'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  plottype='power'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
%   PlotData(data,'variable',plot_var,'plot_type','waveform');
%   plottype='waveform'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
  figure; plot(data(1).time,data(1).model.fixed_variables.E_s1(:,1)); ylabel('s1');
  plottype='input'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
end
data_ac=data;

% POISSON: AC only (lower lambda_ac)
input_type='poisson'; % options: 'poisson','sin','noisy_sin','rectified_sin'
dc=0; ac=.5; tau=2/1000; gINPUT=5;
s_=ApplyModifications(s,{'E','input',0;'Is','input',0;'If','input',0});
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'num_repetitions',num_repetitions,solver_options{:});
if plot_flag
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  plottype='rastergram'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  plottype='power'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
%   PlotData(data,'variable',plot_var,'plot_type','waveform');
%   plottype='waveform'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
  figure; plot(data(1).time,data(1).model.fixed_variables.E_s1(:,1)); ylabel('s1');
  plottype='input'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
end
data_ac=data;

% POISSON: AC + DC
% f1=f0+[-5 -2.5 0 2.5 5];%[-5 0 5];%(-1:1);
% f2=f0+[-5 -2.5 0 2.5 5];
num_repetitions=1;
input_type='poisson'; % options: 'poisson','sin','noisy_sin','rectified_sin'
dc=.01; ac=1; tau=2/1000; gINPUT=5;
s_=ApplyModifications(s,{'E','input',0;'Is','input',0;'If','input',0});
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'dc',dc,'ac',ac,'tau',tau,'num_repetitions',num_repetitions,solver_options{:});
if plot_flag
  PlotData(data,'variable',plot_var,'plot_type','rastergram');
  plottype='rastergram'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
  PlotData(data,'variable',plot_var,'plot_type','power','xlim',[0 100]);
  plottype='power'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
%   PlotData(data,'variable',plot_var,'plot_type','waveform');
%   plottype='waveform'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
  figure; plot(data(1).time,data(1).model.fixed_variables.E_s1(:,1)); ylabel('s1');
  plottype='input'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
end
data_acdc=data;
%}

% Questions for decisions group:
% Q: input? (current injection vs poisson; current poisson vs more complicated?)
% Q: does it make more sense to match spike rates by randomly selecting spikes from the inputs with faster modulation frequency?
% Q: should I determine f0 with a tonic DC input or a poisson DC input?
% - talk to Salva about poissrnd units

% note: poisson s1 produces I(t)=gINPUT*s1*V where injections are directly
% I(t)=gINPUT*I1.

f1=f0+[-5 0 5];
f2=f0+[-5 0 5];
num_repetitions=3;

input_type='rectified_sin';
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});

PlotData(data,'variable',plot_var,'plot_type','rastergram');
PlotData(SelectData(data,'varied',{'E_realization',1}),'variable',plot_var,'plot_type','rastergram');
PlotData(SelectData(data,'varied',{'E_realization',2}),'variable',plot_var,'plot_type','rastergram');
PlotData(SelectData(data,'varied',{'E_realization',3}),'variable',plot_var,'plot_type','rastergram');
%   PlotData(SelectData(data,'varied',{'E_realization',1}),'variable',plot_var,'plot_type','rastergram');
%   plottype='rastergram_realization1'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
%   PlotData(SelectData(data,'varied',{'E_realization',1}),'variable',plot_var,'plot_type','power','xlim',[0 100]);
%   plottype='power_realization1'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');
%   PlotData(SelectData(data,'varied',{'E_realization',1}),'variable',plot_var,'plot_type','waveform');
%   plottype='waveform_realization1'; print(gcf,sprintf('%s-amp%g_DC%g-AC%g-tau%g_%g-%gHz_%s.jpg',input_type,gINPUT,dc,ac,tau,f1(1),f1(end),plottype),'-djpeg');

%% test: Salva Poisson input

% homogeneous poisson input
f=0; dc=0; ac=1; tau=2; baseline=.1;
S=get_input('poisson',nE,T,f,dc,ac,tau,.5,baseline);
PlotData(S,'plot_type','waveform','variable','data');
PlotData(S,'plot_type','power','variable','data','xlim',[0 100]);

% nonhomogeneous poisson input
tspan=[0 1000];
on=0; off=1000; latency=.1; 
baseline=.1; dc=1; ac=5; tau=2; kick=1; interval=tspan(2);
f=40/1000; fspread=.03; phase=0; xc=.25; consigma=.001;
s=getGenExtPoissonTotalGating(on,off,latency,f,fspread,phase,consigma,baseline,dc,ac,tau,kick,Npop,interval,dt,xc);
PlotData(S','plot_type','waveform','variable','data');


Npop=8;
dt=.01;
tOn_pfcInp = 0;               % in ms
tOff_pfcInp = 1000;           % in ms
latency_pfcInp = .1;           % in ms (rise and decay time constant for the input)
freq_pfcInp = .04;              % kHz, sequence of frequencies, e.g. 0 kHz only DC component, 10 Hz modulation (alpha) or 25 Hz (beta) or 35 Hz (low gamma), or first alpha, then beta/gamma, etc.
normFreqSigma_pfcInp = 0.03;  % normalized Freq sigma.
  % increase to get broader spectrum (decrease to get more pure frequency component)
  % todo: change code: remove scaling by input frequency
phase_pfcInp = 0;             % useful to introduce a phase lag between different components
widthSigma_pfcInp = 0.001;    % 0.001 represents an abrupt connectivity transition (in contrast to 0.1; it only applies to dc+ac not the baseline)
rate_pfcInp_baseline = .1;     % kHz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
rate_pfcInp_dc = 1;           % kHz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
rate_pfcInp_ac = 5;
tau_pfcInp = 2;               % ms, exponential decay time constant (AMPA)
kick_pfcInp = 1;              % conductance increase after a spike
interval = 1000;              % interval (ms)
g_pfcInp = 0;                 % mS/cm^2, external conductance (rate normalized)
E_pfcInp = 0;                 % reversal potential (mV; AMPA)
x_c=.25;
freq_pfcInp = .04;              % kHz, sequence of frequencies, e.g. 0 kHz only DC component, 10 Hz modulation (alpha) or 25 Hz (beta) or 35 Hz (low gamma), or first alpha, then beta/gamma, etc.
S = getGenExtPoissonTotalGating(tOn_pfcInp,tOff_pfcInp,latency_pfcInp,freq_pfcInp,normFreqSigma_pfcInp,phase_pfcInp,widthSigma_pfcInp,rate_pfcInp_baseline,rate_pfcInp_dc,rate_pfcInp_ac,tau_pfcInp,kick_pfcInp,Npop,interval,dt,x_c);
PlotData(S','plot_type','waveform','variable','data');
x_c=.75;
freq_pfcInp = .02;              % kHz, sequence of frequencies, e.g. 0 kHz only DC component, 10 Hz modulation (alpha) or 25 Hz (beta) or 35 Hz (low gamma), or first alpha, then beta/gamma, etc.
S = getGenExtPoissonTotalGating(tOn_pfcInp,tOff_pfcInp,latency_pfcInp,freq_pfcInp,normFreqSigma_pfcInp,phase_pfcInp,widthSigma_pfcInp,rate_pfcInp_baseline,rate_pfcInp_dc,rate_pfcInp_ac,tau_pfcInp,kick_pfcInp,Npop,interval,dt,x_c);
PlotData(S','plot_type','waveform','variable','data');
PlotData(S','plot_type','power','variable','data');


% rows: simultaneous inputs
% cols: different epochs

Npop=5;
dt=.01;
tOn_pfcInp = 0;               % in ms
tOff_pfcInp = 1000;           % in ms
latency_pfcInp = [.1;.1];           % in ms (rise and decay time constant for the input)
freq_pfcInp = [.02;.04];              % kHz, sequence of frequencies, e.g. 0 kHz only DC component, 10 Hz modulation (alpha) or 25 Hz (beta) or 35 Hz (low gamma), or first alpha, then beta/gamma, etc.
normFreqSigma_pfcInp = [0.03;.03];  % normalized Freq sigma.
phase_pfcInp = [0;0];             % useful to introduce a phase lag between different components
widthSigma_pfcInp = 0.001;    % 0.001 represents an abrupt connectivity transition (in contrast to 0.1; it only applies to dc+ac not the baseline)
rate_pfcInp_baseline = .1;     % kHz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
rate_pfcInp_dc = [1;1];           % kHz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
rate_pfcInp_ac = [5;5];
tau_pfcInp = 2;               % ms, exponential decay time constant (AMPA)
kick_pfcInp = 1;              % conductance increase after a spike
interval = 2000;              % interval (ms)
g_pfcInp = 0;                 % mS/cm^2, external conductance (rate normalized)
E_pfcInp = 0;                 % reversal potential (mV; AMPA)
x_c=[.25;.75];
S = getGenExtPoissonTotalGating(tOn_pfcInp,tOff_pfcInp,latency_pfcInp,freq_pfcInp,normFreqSigma_pfcInp,phase_pfcInp,widthSigma_pfcInp,rate_pfcInp_baseline,rate_pfcInp_dc,rate_pfcInp_ac,tau_pfcInp,kick_pfcInp,Npop,interval,dt,x_c);
PlotData(S','plot_type','waveform','variable','data');
PlotData(S','plot_type','power','variable','data');


%% Izhikevich study of computational properties

% - create Izhikevich study and add to tutorial
%   based on: http://www.izhikevich.org/publications/izhikevich.m

eqns={
  'a=.02; b=.2; c=-65; d=6; I=14';
  'dv/dt=.04*v^2+5*v+140-u+I; v(0)=-70';
  'du/dt=a*(b*v-u); u(0)=-20';
  'if(v>=30)(v=c;u=u+d)';
  };
P='pop1'; % name of population
vary={
  {P,'a',.02; P,'b',.2 ; P,'c',-50; P,'d',2;  P,'I',15} % tonic bursting
  {P,'a',.01; P,'b',.2 ; P,'c',-65; P,'d',8;  P,'I',30} % spike frequency adaptation
  {P,'a',.02; P,'b',.2 ; P,'c',-65; P,'d',6;  P,'I',7}  % spike latency
  {P,'a',.03; P,'b',.25; P,'c',-52; P,'d',0;  P,'I',0}  % rebound burst
  {P,'a',1;   P,'b',1.5; P,'c',-60; P,'d',0;  P,'I',-65}%bistability
  {P,'a',.02; P,'b',1  ; P,'c',-55; P,'d',4;  P,'I',1}  % accomodation
  {P,'a',-.02;P,'b',-1 ; P,'c',-60; P,'d',8;  P,'I',80} % inhibition-induced spiking
  {P,'a',-.026;P,'b',-1; P,'c',-45; P,'d',0;  P,'I',70} % inhibition-induced bursting
  };
data=SimulateModel(eqns,'tspan',[0 250],'vary',vary);
PlotData(data);


%% test: different inputs with custom external function

f1=f0+[-5 -2.5 0 2.5 5];
f2=f0+[-5 -2.5 0 2.5 5];
num_repetitions=1;

input_type='sin';
data1=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});
PlotData(data1,'variable',plot_var,'plot_type','rastergram');
PlotData(data1,'variable',plot_var,'plot_type','power','xlim',[0 100]);
figure; plot(data1(1).time,data1(1).model.fixed_variables.E_I1)

input_type='rectified_sin';
data2=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});
PlotData(data2,'variable',plot_var,'plot_type','rastergram');
PlotData(data2,'variable',plot_var,'plot_type','power','xlim',[0 100]);
figure; plot(data2(1).time,data2(1).model.fixed_variables.E_I1)

d=SelectData(data1,'time_limits',[100 200]); [d(2).time(1),d(2).time(end)]
d=SelectData(data1,'varied',{'E_f1',f0+[0 5]})
[data1.E_f1],[d.E_f1]
[data1.E_f2],[d.E_f2]
PlotData(d,'variable',plot_var,'plot_type','rastergram');

%% test: multiple realizations and select data with varied subsets

f1=f0+[-5 0 5];
f2=f0+[-5 0 5];
num_repetitions=3;

input_type='sin';
data=ProbeTwoRhythms(s_,'f1',f1,'f2',f2,'gINPUT',gINPUT,'SOA',0,'target',target,'input_type',input_type,'num_repetitions',num_repetitions,solver_options{:});

PlotData(data,'variable',plot_var,'plot_type','rastergram');
PlotData(SelectData(data,'varied',{'E_realization',1}),'variable',plot_var,'plot_type','rastergram');
PlotData(SelectData(data,'varied',{'E_realization',2}),'variable',plot_var,'plot_type','rastergram');
PlotData(SelectData(data,'varied',{'E_realization',3}),'variable',plot_var,'plot_type','rastergram');

%% test: ProbeSpikeThreshold

s=[];
s.pops=ACC_Ecell_specification;

eqns='dv/dt=@current+input; input=0; {iNa,iK}';
[thresh_amp,thresh_v,thresh_FR] = ProbeSpikeThreshold(eqns)
% confirmation test:
data=SimulateModel(s,'vary',{'','input',0:thresh_amp+2});
% plot waveforms
PlotData(data)
% add line at voltage threshold
line(xlim,[thresh_v thresh_v])
% plot spike rates vs input amplitude
PlotFR(data) % threshold amp: 5
% add line at amplitude threshold
line([thresh_amp thresh_amp],ylim);


%% test: plotting waveforms, spectra, and rastergrams

spec=[];
spec.populations(1).name='E1';
spec.populations(1).equations='dv[5]/dt=@current+amp*(1+randn(1,Npop)); amp=10; {iNa,iK}';
spec.populations(2).name='E2';
spec.populations(2).equations='dv[2]/dt=@current; {iNa,iK}';
spec.connections(1).direction='E1->E2';
spec.connections(1).mechanism_list='iAMPA';
vary={'E1','amp',[0 10 20]; 'E1->E2','gSYN',[0 .05 .1]};
data=SimulateModel(spec,'vary',vary,'tspan',[0 800]);
% plot voltage waveforms
PlotData(data,'variable','v','plot_type','power');
% plot voltage power spectrum
PlotData(data,'variable','v','plot_type','waveform');
% plot voltage-derived rastergram
PlotData(data,'variable','v','plot_type','rastergram');

% todo: use the profiler to determine what is taking so long!
% todo: fix indexing for monitors when coder used and pop size=1

% edit(data(1).simulator_options.solve_file)
% p=load('solve/params.mat')

%% test: modifications with external functions, commas, and single quotes

eqns='dv/dt=-v+I(k,:); I=get_input(0,Npop,T,f); f=5';
data=SimulateModel(eqns,'tspan',[0 1000]);
figure; plot(data.time,data.model.fixed_variables.pop1_I)

eqns='dv/dt=-v+I(k,:); I=get_input(''rectified_sin'',Npop,T,f); f=5';
data=SimulateModel(eqns,'tspan',[0 1000]);
figure; plot(data.time,data.model.fixed_variables.pop1_I)

s_=ApplyModifications(s,{'E','equations','cat(ODE1,+I1(k,:)+I2(k,:)'});
s_.populations(1).equations

input_type='rectified_sin';
inp1_str=sprintf('I1=get_input(''%s'',Npop,T,f1)',input_type);
inp2_str=sprintf('I2=get_input(''%s'',Npop,T,f2)',input_type);
eqn_mods=sprintf('cat(ODE1,+I1(k,:)+I2(k,:); %s; %s)',inp1_str,inp2_str);
s_=ApplyModifications(s,{'E','equations',eqn_mods});
s_.populations(1).equations
vary={'E','f1',5;'E','f2',5};
data=SimulateModel(s_,'tspan',[0 1000],'vary',vary)
figure; plot(data.time,data.model.fixed_variables.E_I1)

modifications={...
'E','equations','cat(ODE1,+K1.*(@I1)+K2.*(@I2))';
'E','mechanism_list','+(pfcPoissonInput@I1)';
'E','mechanism_list','+(pfcPoissonInput@I2)'};
s_=ApplyModifications(s,modifications);
s_.populations(1).equations
s_.populations(1).mechanism_list

%% -------------------------------------------------------------------------
% Run Simulations
tic
data=SimulateModel(s,'tspan',tspan,'dt',dt,'solver',solver,'experiment',experiment,'compile_flag',compile_flag,'verbose_flag',verbose_flag);
toc
figure; 
subplot(3,1,1); plot(data.time,data.E_v); ylabel('E_v');
subplot(3,1,2); plot(data.time,data.If_v); ylabel('If_v');
subplot(3,1,3); plot(data.time,data.Is_v); ylabel('Is_v');
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Search space parameters
vary_gSYN_IfE=0:.1:.5;  % E->If:{gSYN}
vary_gSYN_IsE=0;        % E->Is:{gSYN}
% vary parameters
vary={'E->If','gSYN',vary_gSYN_IfE;'E->Is','gSYN',vary_gSYN_IsE};
tic
data=SimulateModel(s,'vary',vary,'compile_flag',compile_flag,'verbose_flag',verbose_flag);
toc
PlotFR(data,'variable',plot_var,'bin_size',bin_size,'bin_shift',bin_shift)
% -------------------------------------------------------------------------

eqns={
  'a=.02; b=.2; c=-65; d=6; I=14';
  'dv/dt=.04*v^2+5*v+140-u+I';
  'du/dt=a*(b*v-u)';
  'if(v>=30)(v=c;u=u+d)';
  };
P='pop1'; % name of population
vary={
  {P,'a',.02; P,'b',.2 ; P,'c',-65; P,'d',6;  P,'I',14} % tonic spiking
  {P,'a',.02; P,'b',.25; P,'c',-65; P,'d',6;  P,'I',.5} % phasic spiking
  {P,'a',.02; P,'b',.2 ; P,'c',-50; P,'d',2;  P,'I',15} % tonic bursting
  {P,'a',.02; P,'b',.25; P,'c',-65; P,'d',.05;P,'I',.6} % phasic bursting
  {P,'a',.02; P,'b',.2 ; P,'c',-65; P,'d',4;  P,'I',10} % mixed mode  
  {P,'a',.01; P,'b',.2 ; P,'c',-65; P,'d',8;  P,'I',30} % spike frequency adaptation
  {P,'a',.02; P,'b',-.1; P,'c',-55; P,'d',6;  P,'I',0}  % Class 1
  {P,'a',.2;  P,'b',.26; P,'c',-65; P,'d',0;  P,'I',0}  % Class 2
  {P,'a',.02; P,'b',.2 ; P,'c',-65; P,'d',6;  P,'I',7}  % spike latency
  {P,'a',.05; P,'b',.26; P,'c',-60; P,'d',0;  P,'I',0}  % subthreshold oscillations
  {P,'a',.1;  P,'b',.26; P,'c',-60; P,'d',-1; P,'I',0}  % resonator
  {P,'a',.02; P,'b',-.1; P,'c',-55; P,'d',6;  P,'I',0}  % integrator
  {P,'a',.03; P,'b',.25; P,'c',-60; P,'d',4;  P,'I',0}  % rebound spike
  {P,'a',.03; P,'b',.25; P,'c',-52; P,'d',0;  P,'I',0}  % rebound burst
  {P,'a',.03; P,'b',.25; P,'c',-60; P,'d',4;  P,'I',0}  % threshold variability
  {P,'a',1;   P,'b',1.5; P,'c',-60; P,'d',0;  P,'I',-65}%bistability
  {P,'a',1;   P,'b',.2 ; P,'c',-60; P,'d',-21;P,'I',0}  % DAP
  {P,'a',.02; P,'b',1  ; P,'c',-55; P,'d',4;  P,'I',0}  % accomodation
  {P,'a',-.02;P,'b',-1 ; P,'c',-60; P,'d',8;  P,'I',80} % inhibition-induced spiking
  {P,'a',-.026;P,'b',-1; P,'c',-45; P,'d',0;  P,'I',80} % inhibition-induced bursting
  };
data=SimulateModel(eqns,'tspan',[0 1000],'vary',vary);
PlotData(data);


%% References:

% MAX CONDUCTANCE DISTRIBUTIONS 

% modeling: draw gmax from uniform then fit some by searching normal or
% lognormal distributed values around candidates:
% O’Leary, T., Williams, A. H., Franci, A., & Marder, E. (2014). Cell types, network homeostasis, and pathological compensation from a biologically plausible ion channel expression model. Neuron, 82(4), 809-821.
% SI: http://www.cell.com/cms/attachment/2015989148/2036876781/mmc1.pdf

% experiment + modeling: complicated, choose uniform spanning realistic range
% O'Leary, T., Williams, A. H., Caplan, J. S., & Marder, E. (2013). Correlations in ion channel expression emerge from homeostatic tuning rules. Proceedings of the National Academy of Sciences, 110(28), E2645-E2654.
% http://www.pnas.org/content/110/28/E2645.full

% motivation for lognormal distribution of synaptic conductances:
% Sarid, Leora, et al. "Modeling a layer 4-to-layer 2/3 module of a single column in rat neocortex: interweaving in vitro and in vivo experimental observations." Proceedings of the National Academy of Sciences 104.41 (2007): 16353-16358.
% http://www.pnas.org/content/104/41/16353.full

% whole-cell currents through AMPARs were lognormally distributed:
% Smith, T. Caitlin, Lu-Yang Wang, and James R. Howe. "Heterogeneous conductance levels of native AMPA receptors." The Journal of Neuroscience 20.6 (2000): 2073-2085.
% http://www.jneurosci.org/content/20/6/2073.full



%{
% homogeneous parameters:   sigma=0
% heterogeneous parameters: sigma>0

% uniform distribution:
unifrnd(gmax-sigma,gmax+sigma)

% normal distribution:
normrnd(gmax,sigma)

% lognormal distribution:
variance=sigma^2;
MU = log(gmax^2/sqrt(variance+gmax^2));
SIGMA = sqrt(log(variance/gmax^2+1));
lognrnd(MU,SIGMA)

% lognormal distribution of synaptic weights
g_SYN=1.09
gvar = [1]
gbar = g_SYN / mean(sum(netcon,1))
MU = log(gbar^2 / sqrt(gvar+gbar^2))
SIGMA = sqrt(log(gvar/gbar^2 + 1))
gmax = lognrnd(MU,SIGMA,Npop,1)
ISYN(V,s) = gmax.*((kernel).*(netcon*s).*(V-E_SYN))

%}
