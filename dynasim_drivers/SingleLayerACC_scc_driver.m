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
tspan=[0 2000];  % [ms], time limits on numerical integration
dt=.01;         % [ms], time step for numerical integration
solver='rk2';   % {options; 'rk4', 'rk2', 'rk1'}, solver to use
compile_flag=1; % whether to compile the simulation
verbose_flag=1; % whether to display log information
downsample_factor=10;
solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag,'downsample_factor',downsample_factor};
% -------------------------------------------------------------------------
%% Population definitions

% population sizes
nE=8; % # of E-cells per population
nI=2; % # of I-cells per population
prefix=sprintf('ACC_1Layer_%gE%gIs%gIf_iNaK',nE,nI,nI);
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


% -------------------------------------------------------------------------
% turn off injected inputs for experiments
modifications={'E','input',0;'Is','input',0;'If','input',0};
s_=ApplyModifications(s,modifications);

% set experimental parameters
gINPUT=.01:.01:.08; fINPUT=10:10:60; num_repetitions=5;
dc=.1; ac=1; baseline=0; tau=2; target='E';

% set simulation and analysis parameters
cluster_flag=1; sims_per_job=25; mem_limit='8G';
post_downsample_factor=1;
proc_variables={'E_v','If_v','Is_v'}; % []

f=2; % Hs
dt=.01; % ms
t=0:(dt/1000):1/f; % one cycle
n=length(t);
dc_nonhomo=.1; 
ac=1;
I = max(0,ac*sin(2*pi*f*t)+dc_nonhomo);
Imean = sum( I ) / n;
dc_homo=Imean;

nf=length(fINPUT);
dc_homo = zeros(1,nf);
for i=1:nf
  t=0:(dt/1000):1/fINPUT(i); % one cycle
  I = max(0,ac*sin(2*pi*fINPUT(i)*t)+dc);
  dc_homo(i) = sum(I)/length(t);
end
dc_homo

f=1:5; num_repetitions=1; gINPUT=[.01 .02 .06 .1];
[model,vary]=PrepProbeResonance(s_,'f',0,'dc',dc_homo,'ac',ac,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'tau',tau,'baseline',baseline);
[dat1,studyinfo]=SimulateModel(model,'vary',vary,solver_options{:});
[model,vary]=PrepProbeResonance(s_,'f',f,'dc',dc_nonhomo,'ac',ac,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'tau',tau,'baseline',baseline);
[dat2,studyinfo]=SimulateModel(model,'vary',vary,solver_options{:});
for i=1:4,mean(dat1(i).model.fixed_variables.E_s(:)),end
for i=1:4,mean(dat2(i).model.fixed_variables.E_s(:)),end

% match the mean strengths of the homogeneous and nonhomogeneous
% poisson-based inputs by adjusting the dc component to match the means 
% of their lambda's:
baseline=0; dc_nonhomo=.1; ac=1;
t=0:1e-5:1; % one cycle
I=max(0,ac*sin(2*pi*t)+dc_nonhomo);
dc_homo=sum(I)/length(t);

% -------------------------------------------------------------------------
% EXPERIMENT #1: Homogeneous Poisson
study_dir_exp1=sprintf('dynasim_studies/%s_TonicPoisson',prefix);
[model,vary]=PrepProbeResonance(s_,'f',0,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'dc',dc_homo,'ac',ac,'tau',tau,'baseline',baseline);
[data1,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir_exp1,solver_options{:},'cluster_flag',cluster_flag,'sims_per_job',sims_per_job,'memory_limit',mem_limit,'prefix',prefix);
% -------------------------------------------------------------------------
% EXPERIMENT #2: Nonhomogeneous Poisson (rhythmic)
study_dir_exp2=sprintf('dynasim_studies/%s_OneRhythmPoisson',prefix);
[model,vary]=PrepProbeResonance(s_,'f',fINPUT,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'dc',dc_nonhomo,'ac',ac,'tau',tau,'baseline',baseline);
[data2,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir_exp2,solver_options{:},'cluster_flag',cluster_flag,'sims_per_job',sims_per_job,'memory_limit',mem_limit,'prefix',prefix);
% -------------------------------------------------------------------------
% EXPERIMENT #3: Two rhythmic Poisson inputs
study_dir_exp3=sprintf('dynasim_studies/%s_TwoRhythmPoisson',prefix);
[model,vary]=PrepProbeTwoRhythms(s_,'f1',fINPUT,'f2',fINPUT,'gINPUT',gINPUT,'target',target,'num_repetitions',num_repetitions,'dc',dc_nonhomo,'ac',ac,'tau',tau,'baseline',baseline);
[data3,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir_exp3,solver_options{:},'cluster_flag',cluster_flag,'sims_per_job',sims_per_job,'memory_limit',mem_limit,'prefix',prefix);
% -------------------------------------------------------------------------
% POST-PROCESSING
if cluster_flag
  data1=ImportData(study_dir_exp1);%,'variables',proc_variables);
  data2=ImportData(study_dir_exp2);%,'variables',proc_variables);
  data3=ImportData(study_dir_exp3);%,'variables',proc_variables);
end
% experiment 1: calc network fMUA=f(gINPUT)
[stats_f0,data1]=CalcPopulationStats(data1,'repetition_parameter',[target '_repetition']);
% experiment 2: calc network fc=f(gINPUT)
[stats_fc,data2]=CalcPopulationStats(data2,'sweep_parameter',[target '_f'],'repetition_parameter',[target '_repetition']);
% experiment 3: calc ROI responses
N=model.specification.populations(1).size;
variable=data(1).labels{1};
dat1=SelectData(data3,'roi',{variable,find(data3(1).model.fixed_variables.E_K1>0)});
dat2=SelectData(data3,'roi',{variable,find(data3(1).model.fixed_variables.E_K2>0)});
clear stats
stats(1)=CalcPopulationStats(dat1,'variable',variable,'repetition_parameter',[target '_repetition']);
stats(2)=CalcPopulationStats(dat2,'variable',variable,'repetition_parameter',[target '_repetition']);  
clear dat1 dat2
% downsample for plotting
if post_downsample_factor>1
  data1=DownsampleData(data1,post_downsample_factor);
  data2=DownsampleData(data2,post_downsample_factor);
  data3=DownsampleData(data3,post_downsample_factor);
end

% compute mean inputs for experiments 1 and 2
inds1=1:num_repetitions:length(data1);
inds2=1:(num_repetitions*length(fINPUT)):length(data2);
Smean_f0=zeros(1,length(gINPUT)); 
Smean_fc=zeros(1,length(gINPUT)); 
for i=1:length(gINPUT)
  Smean_f0(i)=gINPUT(i)*mean(data1(inds1(i)).model.fixed_variables.E_s(:));
  Smean_fc(i)=gINPUT(i)*mean(data2(inds2(i)).model.fixed_variables.E_s(:));
end
% compare natural freq f0=f(gINPUT) vs resonance freq fc=f(gINPUT):
f0=stats_f0.E_v_Power_MUA.repetition_sets.PeakFreq_mu;
f0e=stats_f0.E_v_Power_MUA.repetition_sets.PeakFreq_sd/sqrt(stats_f0.num_repetitions);
fc=stats_fc.E_v_FR.sweep_sets.repetition_sets.sweep_pop_FR_max_mu;
fce=stats_fc.E_v_FR.sweep_sets.repetition_sets.sweep_pop_FR_max_sd/sqrt(stats_fc.num_repetitions);

x=gINPUT;
y1=fc(1:length(gINPUT)); y1e=fce(1:length(gINPUT));
y2=f0(1:length(gINPUT)); y2e=f0e(1:length(gINPUT));
figure; 
subplot(2,2,1); plot_CI(x,y1,y1e,'r'); hold on; plot_CI(x,y2,y2e,'b'); 
xlabel('gINPUT'); ylabel('freq [Hz]'); legend('fc','fMUA');
subplot(2,2,3); plot(x,Smean_fc,'r',x,Smean_f0,'b','linewidth',3)
xlabel('gINPUT'); ylabel('<G(t)>'); legend('for fc','for fMUA');
subplot(2,2,[2 4]); plot(Smean_fc,fc,'r',Smean_f0,f0,'b','linewidth',3);
xlabel('<G(t)>'); ylabel('freq [Hz]'); legend('fc','fMUA');
plottype='Exps1-2_fc-vs-fMUA'; print(gcf,sprintf('%s_gINPUT%g-%g_f-%g-%gHz_%s.jpg',prefix,gINPUT(1),gINPUT(end),f(1),f(end),plottype),'-djpeg');

return

% set amplitude of network input
gINPUT=5;
plot_flag=1;
% todo: base on drive to achieve natural frequency near a target value

% observation:
  % gINPUT<.04: population driven at freq closer to natural freq dominates very significantly
  % gINPUT=.05: ~balanced
  % gINPUT .06-.09: faster has slight advantage
  % gINPUT >1: slower has slight advantage

% space of interest:
% f1=[20 40 60]
% gINPUT: .01-.08

% ------------------------------------
tstart=tic;

% turn off injected inputs
modifications={'E','input',0;'Is','input',0;'If','input',0};
s_=ApplyModifications(s,modifications);

target='E'; input_type='poisson';

% Experiment #1: homogeneous poisson probing natural frequency
gINPUT=0:.01:.08; f=0; num_repetitions=10; dc=0; ac=1; baseline=.1; tau=2;
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
gINPUT=.01:.01:.08; f=[20 40 60]; num_repetitions=10; dc=0; ac=1; baseline=.1; tau=2;
% gINPUT=.01:.01:.08; f=[10 20 30]; num_repetitions=10; dc=0; ac=1; baseline=.1; tau=2;
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
gINPUT=[.01:.01:.08]; f1=[20 40 60]; f2=40; dc=0; ac=1; baseline=.1; tau=2;
% gINPUT=[.01:.01:.08]; f1=[10 20 30]; f2=20; dc=0; ac=1; baseline=.1; tau=2;
num_repetitions=10; target='E'; input_type='poisson'; post_downsample_factor=2;
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
  xlabel(sprintf('f1 (f2=%g)',f2)); ylabel('spike rate'); %legend('E1','E2');
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
  xlabel('mean input'); ylabel('spike rate'); legend('E1','E2');
  
end
plottype='Exps1-3_E1-vs-E2-summary'; print(gcf,sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz_%s.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end),plottype),'-djpeg');

outfile=sprintf('%s_%s-gINPUT%g-%g_f1-%g-%gHz_f2-%g-%gHz.jpg',prefix,input_type,gINPUT(1),gINPUT(end),f1(1),f1(end),f2(1),f2(end))
save(outfile,'stats_f0','stats_fc','stats','E1_Imean','E2_Imean','f1','f2','gINPUT','-v7.3')

% then:
% 1. redo with: (look at frequencies around fMUA)
%   Exp2: gINPUT=.01:.01:.08; f=[10 20 30]; num_repetitions=10; dc=0; ac=1; baseline=.1; tau=2;
%   Exp3: gINPUT=[.01:.01:.08]; f1=[10 20 30]; f2=20; dc=0; ac=1; baseline=.1; tau=2;
% 2. compare rhythmic vs nonrhythmic competition using:
%   f=[0 10 20 30] and [0 20 40 60]

% note:
% freqs around fMUA: f=[10 20 30]
% freqs around fc:   f=[20 40 60]

% ------------------------------------

toc(tstart);