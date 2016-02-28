function pop = ACC_Ecell_specification(type)
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
Cm=1;         % uF/cm2, membrane capacitance
v_IC=-65;     % mV, voltage initial condition
v_IC_noise=1; % mV, scale of normally distributed IC noise
              % IC: v(0)=v_IC+v_IC_noise*randn

% Poisson background activity
sNOISE=sprintf('sNOISE=get_input(''poisson'',Npop,T,0,0,0,tauAMPA,ones(1,Npop),baseline_rate,0,kick); baseline_rate=%g; tauAMPA=%g; kick=%g',baseline_rate,tauAMPA,kick);
iNOISE='noise(t)=-gNOISE.*sNOISE(k,:).*(v-EAMPA); monitor noise';

% population dynamics
eqns=sprintf('dv/dt=(@current+noise(t)+Iapp*(t>=onset&t<=offset))/Cm; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; Iapp=%g; onset=%g; offset=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,v_IC,v_IC_noise,Cm,Iapp,onset,offset,EAMPA,gNOISE,sNOISE);
%eqns=sprintf('dv/dt=(@current+noise(t)+Iapp*(t>=onset&t<=offset))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; Iapp=%g; onset=%g; offset=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,Cm,Iapp,onset,offset,EAMPA,gNOISE,sNOISE);
%eqns=sprintf('dv/dt=(@current+noise(t))/Cm; %s; %s; v(0)=v_IC+v_IC_noise*randn(1,Npop); v_IC=%g; v_IC_noise=%g; Cm=%g; EAMPA=%g; gNOISE=%g; %s',iNOISE,mechanisms,v_IC,v_IC_noise,Cm,EAMPA,gNOISE,sNOISE);

% List of intrinsic mechanisms
biophysics={};
switch type
  case 'ACd_Class1'
    mechanism_list={'iNaF','iKDR','ileak','ih','CaBuffer','iCaT','iKCa'};
    biophysics={'gleak',.15,'gh',.2,'gCaT',1.5,'gKCa',.2,'Eleak',-75,'noise',5};
  case 'PY'
    % todo: use {iNa,iK,iCaT,iAHP,iNaP,iH} {iM?}
    % iNa,iK,iCa
    %mechanism_list={'iNaf','iKdr','iM','iCa','iCan','CaBuffer','ileak'};
    mechanism_list={'iNaf','iKdr','iM','ileak'};
    gNaf=100; gKdr=25; EK=-80; ENa=50; gleak=.3; Eleak=-54.4; gM=.75;
    biophysics={'gNaf',gNaf,'gKdr',gKdr,'gM',gM,'EK',EK,'ENa',ENa,...
      'gleak',gleak,'Eleak',Eleak};
  case 'NaKM'
    mechanism_list={'iNa','iK','iM','ileak'};
  case 'NaKMCa'
    mechanism_list={'iNa','iK','iM','ileak','iCa'};
  case 'PFC_L23' % Durstewitz
    % 'iNa','iK','iCa','iNap','CaDyn','iCan','iL'
    mechanism_list={'iNa','iK','iCa','iNap','CaDyn','iCan','iL'};
  otherwise    
    mechanism_list={'iNa','iK'};
end

% add biophysical parameters to equation string
for i=1:2:length(biophysics)-1
  eqns=sprintf('%s; %s=%g',eqns,biophysics{i},biophysics{i+1});
end

% collect model components in population specification
pop.equations=eqns;
pop.parameters={};
pop.name='';
pop.size=2;
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

