function pop = ACC_Icell_specification(type)
% Generic ACC I-cell model (CB+ or PV+ interneuron)

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
Cm=1;         % uF/cm2, membrane capacitance
v_IC=-65;     % mV, voltage initial condition
v_IC_noise=1; % mV, scale of normally distributed IC noise
              % IC: v(0)=v_IC+v_IC_noise*randn

% Poisson background activity
sNOISE=sprintf('sNOISE=get_input(''poisson'',Npop,T,0,0,0,tauAMPA,ones(1,Npop),baseline_rate,0,kick); baseline_rate=%g; tauAMPA=%g; kick=%g',baseline_rate,tauAMPA,kick);
iNOISE='noise(t)=-gNOISE.*sNOISE(k,:).*(v-EAMPA); monitor noise';

% List of intrinsic mechanisms
biophysics={};
switch type
  case 'HH'
    %mechanisms='{iNa,iK}';
    mechanisms='{iNa,iK,ileak}';
    biophysics={'gleak',.3,'Eleak',-65};
  otherwise
    error('only ''HH'' is supported for I-cells at this time.');
end

% population dynamics
eqns=sprintf('dv/dt=(@current+noise(t)+Iapp*(t>=onset&t<=offset))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; Iapp=%g; onset=%g; offset=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,Cm,Iapp,onset,offset,EAMPA,gNOISE,sNOISE);
%eqns=sprintf('dv/dt=(@current+noise(t))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,Cm,EAMPA,gNOISE,sNOISE);

% add biophysical parameters to equation string
for i=1:2:length(biophysics)-1
  eqns=sprintf('%s; %s=%g',eqns,biophysics{i},biophysics{i+1});
end

% collect model components in population specification
pop.equations=eqns;
pop.parameters={};
pop.name='';
pop.size=2;
pop.mechanism_list={};

% Test background noise:
% data=SimulateModel(eqns,'tspan',[0 1000])
% PlotData(data,'variable',{'v','noise'})
% PlotData(data.model.fixed_variables.pop1_sNOISE)

% Test initial condition noise:
% eqns_=sprintf('dv[10]/dt=(@current+noise(t))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=1; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,EAMPA,gNOISE,sNOISE);
% data=SimulateModel(eqns_,'tspan',[0 1000]); 
% data.pop1_v(1,:)
% PlotData(data)
