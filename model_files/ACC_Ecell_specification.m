function pop = ACC_Ecell_specification(type,nE,heterogeneous_params,heterogeneity_degree,distribution)
if nargin<5, distribution='normal'; end
if nargin<4, heterogeneity_degree=1; end
if nargin<3, heterogeneous_params={}; end
if nargin<2, nE=2; end
if nargin<1, type='ACd_Class1'; end
  
% ACC E-cell model (Pyramidal cell)

% default tonic injected current parameters
Iapp=0;     % uA/cm2, amplitude of injected current
onset=0;    % ms, start time of injected current
offset=inf; % ms, stop time of injected current

% default noise parameters
baseline_rate=.1; % kHz, baseline poisson rate
kick=1;     % mS/cm2, conductance kick per spike
gNOISE=.03;  % mS/cm2, noise scaling factor 
tauAMPA=2;  % ms, decay time constant
EAMPA=0;   % mV, reversal potential
% note: default parameters chosen s.t. baseline cell FR ~ 1-4Hz

% default cell parameters
Cm=1.2;       % uF/cm2, membrane capacitance [Durstewitz 2003]
v_IC=-65;     % mV, voltage initial condition
v_IC_noise=1; % mV, scale of normally distributed IC noise
              % IC: v(0)=v_IC+v_IC_noise*randn

% List of intrinsic mechanisms
biophysics={};
switch type
  case 'ACd_Class1' % ~ Group 2
    baseline_rate=3.5; % chosen s.t. baseline cell FR ~ 1-4Hz
    %mechanism_list={'iNaF','iKDR','ileak','ih','CaBuffer','iCaT','iKCa'};
    %biophysics={'gleak',.4,'gh',gh,'gCaT',2,'gKCa',.2,'Eleak',Eleak,'Cm',1.2}; % noise=5    
    mechanism_list={'iNaF','iKDR','ileak','ih','CaBuffer','iCaT','iKCa','iM'};
    if ismember('gh',heterogeneous_params)
      gh=get_heterogeneous_values(.025,1.4915,distribution,heterogeneity_degree,nE);
      gh=max(0,gh);
    else
      gh=.025;
    end
    if ismember('Eleak',heterogeneous_params)
      Eleak=get_heterogeneous_values(-75,.0997,distribution,heterogeneity_degree,nE);
    else
      Eleak=-75;
    end
    biophysics={'gleak',.15,'gh',gh,'gCaT',0,'gKCa',.2,'Eleak',Eleak,...
      'Cm',1.2,'akdr_scale',.75,'amnaf_scale',2,'gKDR',10,'gNaF',35,'gM',3};
    % Distribution across layers (note: similar to Group 2):
    %    superficial  middle    deep
    %       43%        33%      24%    
    % Model IPs:  AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
    % Class 1      11.6ms            -75       .6mV         f(Inp)    .91 (ARif=.5)
    % Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
    %            94mV       1.7ms      f(Inp)       8Hz
    % Heterogeneity:
    %   gh -> produces appropriately heterogeneous Ih_abssag:
    %       exp_sag_mu=0.885; exp_sag_sd=0.66; level=2*(exp_sag_sd/exp_sag_mu) (1.4915)
    %       gh_mu=.025; gh_sigma=level*gh_mu;
    %       gh=max(0,paramdist(gh_mu,gh_sigma,num_cells));
    %   Eleak -> produces appropriately heterogeneous RMP and Ih_abssag:
    %       exp_RMP_mu=-74.9; exp_RMP_sd=9.96; level=.75*(exp_RMP_sd/exp_RMP_mu); (0.0997)
    %       Eleak_mu=-75; Eleak_sigma=level*Eleak_mu;
    %       Eleak=paramdist(Eleak_mu,Eleak_sigma,num_cells);
  case 'ACd_Class2' % ~ Group 3
    mechanism_list={'iNaF','iKDR','ileak','cadyn','can','iAHP','M','ih'};
    if ismember('gh',heterogeneous_params)
      gh=get_heterogeneous_values(.001,.92,distribution,heterogeneity_degree,nE);
      gh=max(0,gh);
    else
      gh=.001;
    end
    if ismember('Eleak',heterogeneous_params)
      Eleak=get_heterogeneous_values(-70,.0641,distribution,heterogeneity_degree,nE);
    else
      Eleak=-70;
    end    
    biophysics={'hV1NaF',23.1,'hV2NaF',25.1,'gh',gh,'akdr_scale',1,'gAHP',.4,'Eleak',Eleak,'gleak',.017,'gM',.1,'gcan',.0056,'taurCa',1.93,'taudCa',200/7,'gKDR',6,'gNaF',50.25,'CaRest',5e-5,'aAHP_scale',1e6,'Cm',1};
    % Distribution across layers:
    %    superficial  middle    deep
    %       23%        23%      46%
    %   Note: Fiona Group 3 is not present in superficial layers...
    % Model IPs: [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
    % Class 2      26ms              -73mV     .97mV    14ms        .97 (ARif=.89)
    % Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
    %          84mV       1.45ms       f(Inp)     1-2Hz
    % Heterogeneity:
    %   gh -> produces appropriately heterogeneous Ih_abssag:
    %       exp_sag_mu=0.95; exp_sag_sd=0.437; level=2*(exp_sag_sd/exp_sag_mu); (.92)
    %       gh_mu=.001; gh_sigma=level*gh_mu;
    %       gh=max(0,paramdist(gh_mu,gh_sigma,num_cells));
    %   Eleak -> produces appropriately heterogeneous RMP and Ih_abssag:
    %       exp_RMP_mu=-75.4; exp_RMP_sd=6.44; level=.75*(exp_RMP_sd/exp_RMP_mu); (.0641)
    %       Eleak_mu=-70; Eleak_sigma=level*Eleak_mu;
    %       Eleak=paramdist(Eleak_mu,Eleak_sigma,num_cells);
  case 'ACd_Class3' % ~ Group 1
    mechanism_list={'RSiKDR','RSiNaF','iM','ileak','ih'};
    biophysics={'gNaF',200,'gKDR',20,'gleak',.4,'Eleak',-75,'Cm',1.2,'gM',20,'gh',1.5};    
    % Distribution across layers:
    %    superficial  middle    deep
    %       14%        71%      14%
    % Model IPs: [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
    % Class 3      9ms [x]          -70mV     2.5mV       37ms [x]     .89 (ARif=.25)
    % Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
    %          103mV       .42ms [x]    f(Inp)     1-5Hz
  case 'ACd_Class4' % ~ Group 4
    mechanism_list={'iNaF','iKDR','ileak','CaBuffer','iHVA','iAHP','ih'};
    biophysics={'gNaF',75,'gKDR',4,'gleak',.4,'Eleak',-65,'Cm',1.2,'gHVA',2,'gAHP',10,'gh',10};
    % Distribution across layers:
    %    superficial  middle    deep
    %       33%        33%      33%
    % Model IPs: [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
    % Class 4      4.5ms          -62mV     1.7mV       11ms [x]     (stops)
    % Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
    %          82mV       1.5ms       f(Inp)     1-2Hz
  case 'PY'
    % todo: use {iNa,iK,iCaT,iAHP,iNaP,iH} {iM?}
    % iNa,iK,iCa
    %mechanism_list={'iNaf','iKdr','iM','iCa','iCan','CaBuffer','ileak'};
    mechanism_list={'iNaf','iKdr','iM','ileak'};
    gNaf=100; gKdr=25; EK=-80; ENa=50; gleak=.3; Eleak=-54.4; gM=.75;
    biophysics={'gNaf',gNaf,'gKdr',gKdr,'gM',gM,'EK',EK,'ENa',ENa,...
      'gleak',gleak,'Eleak',Eleak,'Cm',1};
  case 'NaKM'
    mechanism_list={'iNa','iK','iM','ileak'};
  case 'NaKMCa'
    mechanism_list={'iNa','iK','iM','ileak','iCa'};
  case 'PFC_L23' % Durstewitz
    % 'iNa','iK','iCa','iNap','CaDyn','iCan','iL'
    % mechanism_list={'iNa','iK','iCa','iNap','CaDyn','iCan','iL'};
    mechanism_list={'iNaF','iKDR','iCa','iNaP','cadyn','iHVA','ileak'};
  otherwise    
    mechanism_list={'iNa','iK'};
end

% Experimentally determined intrinsic properties of ACd cell classes:
%            AHP duration              hyperpol-sag   <ISI>      adaptation
% TARGETS: [AHP_time2trough]   [RMP]   [Ih_abssag] [ISI_median]    [AR23]
%          'AHP(time2trough)'  'RMP'      'Ih'    'median(ISIs)' 'ISI2/ISI3'
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16          .7-1.1
% Class 2      11-32ms       -81 to -69  .5-1.4mV    13-30         .5-.7
% Class 3      30-70ms       -77 to -65  2.4-4.4mV   15-19         .7-1.1
% Class 4      2-17ms        -67 to -60  1.7-2.9mV   22-39        (stops)
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%         'SpikeAmp' 'SpikeWidth'  'ISI1'   'ThreshRate'
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz

% Poisson background activity
sNOISE=sprintf('sNOISE=get_input(''poisson'',Npop,T,0,0,0,tauAMPA,ones(1,Npop),baseline_rate,0,kick); baseline_rate=%g; tauAMPA=%g; kick=%g',baseline_rate,tauAMPA,kick);
iNOISE='noise(t)=-gNOISE.*sNOISE(k,:).*(v-EAMPA); monitor noise';

% population dynamics
eqns=sprintf('dv/dt=(@current+noise(t)+Iapp*(t>=onset&t<=offset))/Cm; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; Iapp=%g; onset=%g; offset=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,v_IC,v_IC_noise,Cm,Iapp,onset,offset,EAMPA,gNOISE,sNOISE);
%eqns=sprintf('dv/dt=(@current+noise(t)+Iapp*(t>=onset&t<=offset))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; Iapp=%g; onset=%g; offset=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,Cm,Iapp,onset,offset,EAMPA,gNOISE,sNOISE);
%eqns=sprintf('dv/dt=(@current+noise(t))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,Cm,EAMPA,gNOISE,sNOISE);

% add biophysical parameters to equation string
parameters={};
for i=1:2:length(biophysics)-1
  if length(biophysics{i+1})==1
    eqns=sprintf('%s; %s=%g',eqns,biophysics{i},biophysics{i+1});
  else
    parameters{end+1}=biophysics{i};
    parameters{end+1}=biophysics{i+1};
  end
  %eqns=sprintf('%s; %s=',eqns,biophysics{i});
  %eqns=[eqns '[' sprintf('%g ',biophysics{i+1}) ']'];
end

% collect model components in population specification
pop.equations=eqns;
pop.parameters=parameters;
pop.name='';
pop.size=nE;
pop.mechanism_list=mechanism_list;

% Test background noise:
% data=SimulateModel(eqns,'tspan',[0 1000])
% PlotData(data,'variable',{'v','noise'})
% PlotData(data.model.fixed_variables.pop1_sNOISE)

% Test initial condition noise:
% eqns_=sprintf('dv[10]/dt=(@current+noise(t))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=1; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,EAMPA,gNOISE,sNOISE);
% data=SimulateModel(eqns_,'tspan',[0 1000]); 
% data.pop1_v(1,:)
% PlotData(data)

function values=get_heterogeneous_values(mu,spread,distribution,degree,num_cells)  
% spread
sigma=abs(degree*spread*mu);

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
values=paramdist(mu,sigma,num_cells);

% % mean parameter value (mean gmax)
% mu=1; 
% 
% % degree of heterogeneity
% sigma=spread*degree;
% % switch degree
% %   case 0          % homogeneous population
% %     sigma=0;
% %   case 1          % 10% spread (standard deviation)
% %     sigma=.1*mu;
% %   case 2          % 25% spread (standard deviation)
% %     sigma=.25*mu;
% % end


%         dE=21; % soma diameter [um]. area=pi*(d/2)^2 dendritic area=4*soma area.
%         SE=pi*(dE/2)^2;
%         cm    =.8;
%         gm    =.54;  Eleak =-80;
%         gNa   =80;   ENa   =55;
%         gK    =36;   EK    =-80;
%         gKs   =25;   EKs   = EK;
%         gNap  = .6;  ENap  = ENa;  napshift=5;
%         gCan  =.07;  ECan  =-20;   tauca=80; alphacaf=.002;  % gCan=.1
%         gCa   =(29.4/SE)*100; ECa=150;
%         EsomaTemplate.mechanisms = {'iL','iNa','iK','iCa','iNap','CaDyn','iCan','randn',InputMechanism};
%         EsomaTemplate.parameters = {...
%           'Cm',Cm,'g_l',gm,'E_l',Eleak,'V_IC',-70,'noise',Enoise,'IC_noise',IC_noise,...
%           'gNa',gNa,'ENa',ENa,'gNap',gNap,'napshift',napshift,'ENap',ENap,'gKf',gK,'EKf',EK,'gCaf',gCa,'ECaf',ECa,... 
%           'gCan',gCan,'ECan',ECan,'tauca',tauca,'alphacaf',alphacaf,...
%           inputparams{:}};

