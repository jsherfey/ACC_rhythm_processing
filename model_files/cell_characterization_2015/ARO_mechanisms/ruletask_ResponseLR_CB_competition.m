% Model: E/I % Base model: /home/jason/models/dnsim/CoherenceEffects/base/sim_transfers.m
cd(fileparts(mfilename('fullpath')));
if ~any(strfind(path,pwd)), addpath(genpath(pwd)); end; cd mechanisms
% ---------------------------------------------------------
% simulation parameters
dt=.01; dsfact=10; usecoder=1; solver='euler'; % solver='rk2'
% tspan=[0 2000];
tspan=[0 7500];
% tspan=[0 3500]; %[450 5000];%[450 3500];
t=tspan(1):dt:tspan(2); nt=length(t); xlims=tspan;

% What to plot?
plot_series_flag=1; skip_plotv=1; summary_plots=1; buildonly=0; 
plotinputs=0; plotconnect=0; plotcurrents_flag=0; visible_status='on'; QC_flag=0; verbose=0;
rulesims_flag=1; plot_raster_only=1; plot_spectrogram_flag=0; plot_coh_flag=0;

noiseless_flag=0; force_all2all=1; 

%{
focus on effects of sync input when r(Left)=r(Right)
for exploring competition in the full model:
EsomaType=[2 9]
...
may need to increase nI2>1 and add fractional connections from response
groups to each I2 cell (in that case, set nECellsPerResponse=5)

for simulated superficial layer inputs:
EsomaType=[9]
nLayers=1
nE=10  (nECellsPerBit=5, nBits=2)
nI=2   (nICellsPerRule=2, nRules=1)
nIs=0  (nIsCellsPerRule=0)

add to template.dynamics: "InputMatrix(t).*(0-V)" for excitatory input matrix
and set:
inpscale=.02;
StimAmpE=1.5*inpscale;
StimAmpI=.75*inpscale; % feedforward input to middle layers (PV+ INs)
StimAmpIs=0*inpscale;
RuleAmpI=0*inpscale;
RuleAmpE=.2*inpscale; % .0; % drive to E-cells to maintain rE | inhibition
RuleAmpIs=17*inpscale;% 12; % feedback input to upper layers (CB+ INs)

%}

% Define task elements
nBits=2;          % 2 or 4 (eg, (ColorLeft,OrientationRight) or (ColorLeft,ColorRight,OrientationLeft,OrientationRight))
nECellsPerBit=5;
nResponses=2;     % 1 or 2 (eg, Left, Right)
nECellsPerResponse=5; % 2
nRules=2;         % 1 or 2 (eg, Color, Orientation)
nICellsPerRule=2;
nIsCellsPerRule=2;
nI2=2;                         % number of fast I-cells (downstream layers)
nIs2=2;                        % number of slow I-cells (downstream layers)
FFfraction=.4; % fraction of a Bit's cells connecting to a given Response cell
EI2fraction=1;%.8;
use_predef_topology=1;

include_thalamus_flag=0; % 0: none, 1: asleep (bursty TC), 2: awake (no burst)
nTCPerResponse=1;
nRE=2; % set to at least 2 for lateral inhibition (decreases rate, more realistic)

% Set population sizes
numLayers=2;                % number of layers
% input layer
nE=nBits*nECellsPerBit;     % number of E-cells (superficial layer)
nI=nRules*nICellsPerRule;   % number of fast I-cells (superficial layer)
nIs=nRules*nIsCellsPerRule; % number of slow I-cells (superficial layer)
% downstream layers
nE2=nResponses*nECellsPerResponse; % number of E-cells (downstream layers)
nTC=nResponses*nTCPerResponse;

% Define ROIs for inputs
CatA=1:nECellsPerBit;             % ColorLeft
CatB=nECellsPerBit+CatA;          % ColorRight
CatC=nECellsPerBit+CatB;          % OrientationLeft
CatD=nECellsPerBit+CatC;          % OrientationRight
Context1=1:nICellsPerRule;        % Color Rule
Context2=nICellsPerRule+Context1; % Orientation Rule
if nBits==2
  CatC=[]; CatD=[];
  EStimTargets=[CatA,CatB]; % incongruent stimulus
  ERule1Targets=CatA;
  ERule2Targets=CatB;
elseif nBits==4
  EStimTargets=[CatA,CatD]; % incongruent stimulus
  %EStimTargets=[CatB,CatD]; % congruent stimulus
  ERule1Targets=[CatA,CatB];
  ERule2Targets=[CatC,CatD];
else
  error('nBits must be set to 2 or 4');
end
Response1=1:nECellsPerResponse;         % Left
Response2=nECellsPerResponse+Response1; % Right
if nResponses<2
  Response2=[];
end
% Define ROIs for post-hoc analysis
Layer1ResponseROIs={[CatA,CatC],[CatB,CatD]};
Layer2ResponseROIs={Response1,Response2};

% Network connectivity
ruletask_ResponseLR_CB_topology;

% Connection strengths
tauE=2; tauNMDA=150; 
gEI=.2; gEI2=1.5;
gIE=.2; gIE2=1; % 1.2
gII=0; gII2=0;

gIIs=0;
gIsI=1.25; gIsI2=gIsI/2;
gEIs=1; gEIs2=gEI2;
gIsE=.35; gIsE2=gIE2; 
gIsIs=.45; gIsIs2=gIsIs;%0;
geeampa=0; geeampa2=1.8; %.14/(nECellsPerBit/2); 
geenmda=.2; geenmda2=1.8; %.5/(nECellsPerBit/2);
gEEampa=1.1; 
gEEnmda=.3;%gEEampa*1.25;
gEdE=.1; gEEd=gEdE/2.2;

% -- Ifast,Islow ASYMMETRIES -----------------
gEInmda=.1; gEInmda2=.2;
tauI=5; tauI2=5; tauIs=13; tauIs2=13;

% Thalamus connectivity
gETC=.25;
gTCE=0;
gTCRE=0.4; % .69 (Destexhe96), .4 (Ching2010)
gRETCa=0.1; % .069 (Destexhe96), .1 (Ching2010)
gRETCb=0.1; % .138 (Destexhe96), .001 (Ching2010)
gRERE=0.138; % .138 (Destexhe96), .1 (Ching2010)

% ############################################
% Input parameters
% ############################################
% -------------
% CASE 1 trial: Context --> (Context+Stimulus) --> off
% -------------
Context1Time=[500 3500]; S1time=[1000 3500]; % tspan=[0 3500]; or [0 7500]
% % Context1Time=[500 6500]; S1time=[1000 6500]; % tspan=[0 7500] (one trial)
% % Context1Time=[500 1500]; S1time=[1000 1500]; % tspan=[0 2000];
% -------------
% CASE 2 trial: Context --> (Context+Stimulus) --> more Context --> off:
% -------------
% Context1Time=[500 4000]; S1time=[1000 3000]; % tspan=[0 7500] (two trials)
% % Context1Time=[500 3500]; S1time=[1000 3000]; % tspan=[0 3500]; or [0 7500]
% % Context1Time=[500 6500]; S1time=[1000 5500]; % tspan=[0 7500]; or [0 15000]
% -------------
% CASE 3 trial: Stimulus --> Context --> off:
% Context1Time=[2000 3000]; S1time=[500 2000]; % tspan=[0 3500];
% % Context1Time=[3500 6500]; S1time=[500 3500]; % tspan=[0 7500]; or [0 15000]
% -------------

Context2Time=Context1Time+max(S1time(2),Context1Time(2)); 
S2time=S1time+max(S1time(2),Context1Time(2));
TOILIM=[Context1Time(1) S2time(2)]; % for population analysis
StimOn=S1time(1); RuleOn=Context1Time(1); TrialOff=Context1Time(2); % for task analysis

% feedforward input to middle layers
StimAmpE=2;
StimAmpI=1.2; % feedforward input to middle layers (PV+ INs)
StimAmpIs=0;
% feedback input to upper layers
RuleAmpI=0;
RuleAmpE=.2; % .0; % drive to E-cells to maintain rE | inhibition
RuleAmpIs=17;% 12; % feedback input to upper layers (CB+ INs)
kickstart=0;

% filtered poisson parameters
lambda_modulation=1000; % <-- this lambda is always used
lambda_baseline=300;    % <-- additionally, this lambda used if modulation_frequency>0
modulation_frequency=0;
ninputs=1; sharedfraction=0;
inp_tauD=2; inp_tauR=.5;

% noise sources
IC_noise=.1; 
ENOISE=0; % 0,.25
INOISE=0; % 0,.5
Isnoise=0;
EnoiseType9=0; % noise to deep layer E-cells

% -----------------------
% ADDED layer 2 LTS cells
% -----------------------

% gIsI2=1.25; gIsE2=1; gEIs2=1.5; gIsIs2=0;
% gIsI2=0*gIsI; gIsE2=0*gIE2; gEIs2=1*gEI2; gIsIs2=1*gIsIs;
% gIsI2=gIsI/2; gIsE2=1*gIE2; gEIs2=1*gEI2; gIsIs2=1*gIsIs;

gIsI2=gIsI;%1.25; %1.25; %1 %.75*gIsI;
gIsE2=gIE2;%1.25; %1.25; %1 %.75*gIE2;

modulation_frequency=25;
VAR2='gIsE2'; VAL2=0;%[0:.25:2];%VAL2=[10:40]; %26,40  %[5:60];
VAR1='iteration'; VAL1=1;%1:5;

% VAR2='modulation_frequency'; VAL2=[10:40]; %26,40  %[5:60];
% VAR1='iteration'; VAL1=1:5;

% -----------------------
% ADDED (TC,RE) cells
% -----------------------
if include_thalamus_flag==1 % asleep
  VAR2='gETC'; VAL2=.25;
  VAR1='gTCE'; VAL1=.0;
elseif include_thalamus_flag==2 % awake
  VAR2='gETC'; VAL2=1.2;
  VAR1='gTCE'; VAL1=.0;%.2; %0
  gRETCa=.069; gRETCb=.138; gRERE=0.138; % Destexhe96
%   gRETCa=.1; gRETCb=.001; gRERE=0.1; % Ching2010
%   %gRETCa=0; gRETCb=0; gRERE=.138;
end

auxvar2=''; auxval2=[];
auxvar1=''; auxval1=[];

%{
<------- without thalamus ---------->
VAR2='gIsI2'; VAL2=[0 .25 .5 .75 1 1.25 1.5];

VAR2='gIE2'; VAL2=[.25 .5 .75 1 1.25];
auxvar2='gIsE2'; auxval2=VAL2;

VAR1='iteration'; VAL1=[1 2];
auxvar1=''; auxval1=[];

VAR1='geenmda2'; VAL1=1.8;%1.2;%[0:.5:2];
auxvar1='geeampa2'; auxval1=VAL1; % []

VAR2='modulation_frequency'; VAL2=5:60;
VAR1='iteration'; VAL1=[1:5];

<------- with thalamus ---------->

% VAR2='gETC'; VAL2=.25;
% VAR1='gTCE'; VAL1=.0;

%}






% test: comparison (with and without competition):
% VAR2='gIE2'; VAL2=[0 .5 1 1.5 2 2.5 5];        %[0 .5 1 1.5];%[1.5 2 2.5];%.1;%5:2:19;%[1:5];
% VAR1='iteration'; VAL1=[1:10];      %[5 13];%2.5 [0 .5 1 1.5 2];%.5:.25:2.5; %[0:.25:3];

% test (once gIE2>0 increases sync winner effects):
% VAR2='tauI2'; VAL2=[tauI tauIs];        %[0 .5 1 1.5];%[1.5 2 2.5];%.1;%5:2:19;%[1:5];
% VAR1='iteration'; VAL1=[1:5];      %[5 13];%2.5 [0 .5 1 1.5 2];%.5:.25:2.5; %[0:.25:3];

% salva: increase recurrent excitation and lateral inhibition in deep layer
% --> to get attractor dynamics that go hi and low

% note: gEEnmda amplifies differences between L and R (given no competition between them)

%{

studied:
tauI [5-19]
StimAmpE [.5-2.5]
RuleAmpE [0-1]
StimAmpI [0-2]
RuleAmpI 0
StimAmpIs 0
RuleAmpIs [0-20]

explore: (StimAmpE,StimAmpI) | tauI2=5ms,13ms

%}

% VAR2='junk1'; VAL2=0;
% VAR1='junk2'; VAL1=0;
% iteration=1:10

if plotinputs && plotconnect
  ruletask_ResponseLR_CB_inputs;
  return
end

%{
Save plots (example):
  ftype='routing-by-rhythm-not-rates';
  fname=['../ruletask_ResponseLR_CatAD_context-stim_works4b-LOVELY-SWITCHING-with-delay_no-output-competition_' ftype];
  print(gcf,'-djpeg',[fname '.jpg'])
  saveas(gcf,fname,'fig')
%}

% --------------------------------------------------------------------

% VAR2='StimAmpE'; VAL2=[0:.25:3];%0;.1;1
% VAR1='StimAmpI'; VAL1=[0 .25 .5 .75 1];%0;.1;1

% % nE=2: ASYNC --> SYNC transition (without Working Memory): (DONE! must repeat until converges)
% S1time=[1000 2500],S2time=S1time,Context1Time=[2500 4000],TOILIM=[1000 4000]
% gEI=.2,gIE=.1,StimAmpE=1,RuleAmpE=1.9,StimAmpI=-2.4,RuleAmpI=.25
% SDIM1=1,SDIM2=2,Context1=1,NOTRULE=1, geeampa=0, geenmda=0
% 
% % nE=2: ASYNC --> SYNC transition (with Working Memory): (DONE!)
% S1time=[1000 1500],S2time=S1time,Context1Time=[2500 3000],TOILIM=[1000 4000]
% gEI=.1,gIE=.1,StimAmpE=1,RuleAmpE=1.9,StimAmpI=-2.1,RuleAmpI=.25
% SDIM1=1,SDIM2=2,Context1=1,NOTRULE=1, geeampa=.08, geenmda=.15, gEInmda=0
% 
% % nE=2: Context(Rule) -> Stimulus (with Working Memory): (tspan=[0 3000])
% Context1Time=[500 2500],S1time=[1000 1500]
% Context2Time=[3000 5000],S2time=[3500 4000],TOILIM=[500 3000]
% gEI=.5,gIE=.1,StimAmpE=1,RuleAmpE=.2,StimAmpI=0,RuleAmpI=.6
% SDIM1=1,SDIM2=2,Context1=1,Context2=2,NOTRULE=[],CatA=[1 2],CatB=[1 2]
% geeampa=.08, geenmda=.2, gEInmda=0, gEEampa=.3
% % no rule: gEI=.5,gIE=.1,StimAmpE=1,RuleAmpE=0,StimAmpI=0,RuleAmpI=0
% 
% % nE=4: Coh switch (wip): Context(Rule) -> Stimulus (with Working Memory):
% Context1Time=[500 2500],S1time=[1000 1500]
% Context2Time=[2500 4500],S2time=[3000 3500],TOILIM=[500 4500]
% StimAmpE=1,RuleAmpE=.1,StimAmpI=-2,RuleAmpI=.25
% %gEI=.2,gIE=.1,StimAmpE=.5,RuleAmpE=.15,StimAmpI=0,RuleAmpI=.6
% SDIM1=1:2,SDIM2=3:4,Context1=1,Context2=2,NOTRULE=[],CatA=[1 2],CatB=[3 4]
% gEI=.2,gIE=.2,geeampa=.08, geenmda=.5, gEInmda=0, gEEampa=.8
% 
% Context2Time=[3000 4500],S2time=[3500 4000],TOILIM=[500 4500]
% gEEampa=.85, geeampa=.14,RuleAmpE=.09
% 
% VAR2='StimAmpI'; VAL2=0;%-2;
% VAR1='StimAmpE'; VAL1=3;%1,3

%{

4-cell network:

% ASYNC
StimAmpE=1;%1;%2.5;
StimAmpI=.25;%.5;
VAR2='gEI'; VAL2=.2;%.04:.005:.08;%0:.02:.1;%[0:.1:1];%0;.1;1
VAR1='gIE'; VAL1=.0;%.1;%.04:.005:.08;%0:.02:.1;%[0:.1:1];%0;.1;1

% SYNC
StimAmpE=1.8;%1;%2.5;
StimAmpI=.25;%.5;
VAR2='gEI'; VAL2=.2;%.04:.005:.08;%0:.02:.1;%[0:.1:1];%0;.1;1
VAR1='gIE'; VAL1=.1;%.1;%.04:.005:.08;%0:.02:.1;%[0:.1:1];%0;.1;1

% ASYNC --> SYNC transition (without Working Memory):
S1time=[1000 2500],S2time=S1time,Context1Time=[2500 4000],TOILIM=[1000 4000]
gEI=.2,gIE=.1,StimAmpE=1,RuleAmpE=1.8,StimAmpI=-.05,RuleAmpI=.25
SDIM1=1,SDIM2=2,RULE=1,NOTRULE=1, geeampa=0, geenmda=0

% ASYNC --> SYNC transition (with Working Memory):
S1time=[1000 1500],S2time=S1time,Context1Time=[2500 3000],TOILIM=[1000 4000]
gEI=.1,gIE=.1,StimAmpE=1,RuleAmpE=1.9,StimAmpI=-.04,RuleAmpI=.25
SDIM1=1,SDIM2=2,RULE=1,NOTRULE=1, geeampa=.08, geenmda=.15, gEInmda=0

%}

%{
% assess results
cellfun(@length,spiketimes{1}) % input
cellfun(@length,spiketimes{2}) % output
%}

% connection probabilities
if noiseless_flag || force_all2all
  PrSYNei=1;
  PrSYNie=1;
  PrSYNii=1;
  PrSYNee=1; % within-layer connection probability
  PrSYNEE=1; % between-layer connection probability
else
  PrSYNei=.65;
  PrSYNie=.6;
  PrSYNii=.55;
  PrSYNee=.3; % within-layer connection probability
  PrSYNEE=.3; % between-layer connection probability
end
if noiseless_flag
  IC_noise=0;
  Inoise=0; 
else  
  %IC_noise=.1;
  Inoise=0;
end

% ---------------------------------------------------------
% Input parameters
% InputType=1;%6;
  % 1=tonic input
  % 6=Poisson (Jason's simple method)
  % 7=Poisson (Ben's method)
InputType=1; 
switch InputType
  case 1
    InputMechanism='InputGenerator2';
  case 2
    InputMechanism='tonic';
  case {6,7}
    InputMechanism='iMultiPoissonExp';
  otherwise
    InputMechanism='InputGenerator2';
end

% common input parameters
onset=250;     % ms, stimulus start time
offset=1000;%inf;   % ms, stimulus stop time
IstimFrac=0;  % note: Istim=IstimFrac*Estim % IstimFrac=.2;
inp_targets=nE/2+1;%5; % input-layer E-cells (indices) receiving the input. (eg, 1:nE)
ESTIM=0; % force Estim to this value below

% poisson input parameters
Ninputs=10;       % number of independent inputs
sharedfraction=0; % fraction of inputs shared across population
inp_lambda=5;     % Hz, Poisson spike rate

% exponential synaptic filter
inp_tauD=2;       % ms, exponential decay time
inp_tauR=.5;      % ms, exponential rise time

% ---------------------------------------------------------
% Network parameters
% connectivity
ConnType=1;
switch ConnType
  case 1 % all-to-all
    span=inf;
  case 2 % neighbors
    span=.5;
  case 3 % gaussian
    span=.5;
end
normalize=1;  % whether to normalize weights by the number of presynaptic sources
zerodiag=0;   % whether to remove self connections when connecting a population to itself

% synaptic time constants and reversal potentials
% tauE=2; 
% tauNMDA=150; 
tauNMDArise=1E-5;
% tauI=5;
Egaba=-80;

% synaptic weights
% gEI=.2;%.05;%8;
% gIE=.2;%.5;%(nI/nE)*gEI;
% gII=.5;%.5*gIE;
% geeampa=.5;%0;%.1;              % within-layer E->E gAMPA
% geenmda=geeampa*1.25;%0;%.5;    % within-layer E->E gNMDA
% gEEampa=.15;%1.5*geeampa;%.1;                     % between-layer E->E gAMPA
% gEEnmda=gEEampa*1.25;           % between-layer E->E gNMDA

% compartment connection parameters
% gEEd=2.5; % E->Ed
% gEdE=5.5; % Ed->E
% gEdE=.1; 
% gEEd=gEdE/2.2;

% % connection probabilities
% if noiseless_flag || force_all2all
%   PrSYNei=1;
%   PrSYNie=1;
%   PrSYNii=1;
%   PrSYNee=1; % within-layer connection probability
%   PrSYNEE=1; % between-layer connection probability
% else
%   PrSYNei=.65;
%   PrSYNie=.6;
%   PrSYNii=.55;
%   PrSYNee=.3; % within-layer connection probability
%   PrSYNEE=.3; % between-layer connection probability
% end
% if noiseless_flag
%   IC_noise=0;
%   Inoise=0; 
% else  
%   %IC_noise=.1;
%   Inoise=0;
% end
% ---------------------------------------------------------
% cell model selection ("*" means used in my published work)
  % Durstewitz 832:  (EsomaType=8, EdendType=3, IfastType=2)
  % Papoutsi 7156:   (EsomaType=2/7, EdendType=1/2, IfastType=5, IslowType=6)
  % Kramer XY34:     (EsomaType=X, EdendType=Y, IfastType=3, IslowType=4)
EsomaType=[2 9 2]; %2 '2h'
  % 1: HH model
  % 2*: Tallie ACd ci196 model (only naf, kdr, pas) -- for homogeneous network
  % 2h: same as 2 with h-current
  % 3: Wang/Tegner L2/3 E-cells (only E_iK, E_iNa, E_ileak)
  % 4: Tallie ACd ci196 model (all currents)
  % 5: Tallie ACd ci176h w/o wbNa (all currents except wbNa)
  % 7*: Papoutsi model (all currents, heterogeneous parameters) -- used for Tallie paper
  % 8: (Durstewitz-basedability neuron) (iL,iNa,iK,iCa,iNap,CaDyn,iCan) 'pfcsuper' in routing SESSION006e  
  % 9: (reduced Durstewitz-based neuron) (iL,iNa,iK) - removes iNap/iCan based buildup in excitability
EdendType=0;%1;
  % 0*: no dendrite
  % 1: Papoutsi PYda
  % 2: Papoutsi PYdb
  % 3: (Durstewitz-based neuron) Edend (iNa,iKs,iNap)
IfastType=5;%5; % <-- targets Esoma
  % 1: HH model
  % 2*: Wang-Buzsaki IN
  % 3: Kramer08 FS (basket cell) model
  % 4: Kramer08 LTS model
  % 5: Papoutsi FS
  % 6: Papoutsi CB (Durstewitz)
IslowType=7;%1; % <-- targets Edend
  % 0*: no Islow
  % 1: Papoutsi CB (Durstewitz)
  % 5: Papoutsi FS
  % 7: Kramer LTS modified for DLPFC (longer AP)
TC_Type=1;
  % 1: Austin's model with parameters from Destexhe 1996
RE_Type=1;
  % 1: Austin's model with parameters from Destexhe 1996

% define population indices
if nIs==0 && nIs2==0, IslowType=0; end
if nI==0 && nI2==0, IfastType=0; end
EsomaID=1:numLayers;
id=EsomaID;
if EdendType>0
  EdendID=id+numLayers;
  id=EdendID;
end
if IfastType>0 && (nI>0 || nI2>0)
  if nI2==0
    IfastID=id(end)+1;
  else
    IfastID=id+numLayers;
  end
  id=IfastID;
end
if ~isequal(IslowType,0) && (nIs>0 || nIs2>0)
  if nIs2==0
    IslowID=id(end)+1;
  else
    IslowID=id+numLayers;
  end
  id=IslowID;
end
if include_thalamus_flag
  TC_ID=id(end)+1;
  RE_ID=TC_ID+1;
  id=RE_ID;
end

% ---------------------------------------------------------
% set cell type-specific parameters (adjusted to match response properties)
switch EsomaType(1)
  case 1 % HH model
    Enoise=1.5;   % -> std(VE|stim=0)=8mV
    Estim=5;      % -> rE=32Hz | gIE=0
  case {2,'2h'} % Tallie ACd ci196 model (only naf, kdr, pas)
    Enoise=1;     % -> std(VE|stim=0)=9mV
    if InputType==2 % tonic
      Estim=.5;     % -> rE=32Hz | gIE=0, tonic
      %tauI=13; Estim=2;  % <-- beta2
      %tauI=5; Estim=2;   % <-- gamma
    else
      Estim=.3; % Poisson (InputType==6)
    end
  case 3 % Wang/Tegner L2/3 E-cells (only E_iK, E_iNa, E_ileak)
    Enoise=5;      % -> std(VE|stim=0)=8mV
    if InputType==2 % tonic
      Estim=15;      % -> rE=32Hz | gIE=0, tonic
    else
      Estim=2;       % Poisson input (InputType==6)
    end
    dt=min(dt,.005);
  case 4 % Tallie ACd ci196 model (all currents)
    Enoise=1;      % -> std(VE|stim=0)=5mV
    Estim=3.5;     % -> rE=32Hz | gIE=0
    %Estim=0;%Estim=.5;%Estim=2;
    %tauI=13; Estim=2;  % <-- beta2
    %tauI=5; Estim=2;   % <-- gamma
  case 5 % Tallie ACd ci176h w/o wbNa (all currents except wbNa)
    Enoise=1.5;         % -> std(VE|stim=0)=8mV
    Estim=1;
    %   Estim=1.75;       % -> rE=32Hz | gIE=0
    %   tauI=13; Estim=2;  % <-- 17Hz (beta1)
    %   tauI=5; Estim=2;   % <-- 33Hz (gamma)
    %   tauI=5; Estim=1; % <-- 22Hz (beta2)|gh=0,.01
  case 6
    Enoise=1.5;     
    Estim=1;
  case 7 % Papoutsi model (all currents, heterogeneous parameters)
    Enoise=1;
    Estim=3;
    % tauI=13; Estim=3;  % <-- 19Hz (beta2)
%     tauI=5; Estim=3;   % <--  33Hz (gamma)
  case 8 % (Durstewitz-based neuron)
    Enoise=1;
    Estim=1.75;   % -> rE=32Hz | gIE=0
    % gEEampa=1; gEI=.2; gIE=.5; gII=.5;
    % 
  otherwise
    Estim=0;
    Enoise=0;
end

switch EdendType
  case {1,2}  % Papoutsi PYda/PYdb
    Estim=2.5;
  case 3      % (Durstewitz-based neuron)
    Estim=1.75;   % -> rE=32Hz | gIE=0
end

if noiseless_flag, Enoise=0; end

% -------------------------------------------------------------------------
% probe noise-driven response:
if 0
  Estim=0;
end
% probe single cell response:
if 0
  nE=1; nI=1; gIE=0;
end
% -------------------------------------------------------------------------
% WHAT TO VARY (across multiple simulations)
setnum=0; % 0=vary nothing
switch setnum
  case 1 % characterize effect of tauI on response to tonic input
    var1='tauI'; % smaller set
    val1=10:15;
    var2='Estim'; % bigger set
    val2=.25:.5:25;
  case 2 % Q: how does inhibition strength affect sparse rhythm?
    var1='gIE';
    val1=0:10;
    var2='Estim';
    val2=0:20;
  case 3 % play
    Estim=2.25;%2.25;
    dt=.01;
    gEI=.5;%.2;
    gIE=1;%.5;
    gII=2;%.5;
    %gEEampa=0;%1;
    sharedfraction=0;
    var1='sharedfraction';
    val1=sharedfraction;% .8 [.4:.1:1];%[0 .25 .5 .75];%[0 1];%0:.25:1;
    var2='gEI';
    val2=gEI; % .2
  otherwise 
    if exist('VAR1','var') && ~isempty(VAR1) && exist('VAL1','var') && ~isempty(VAL1)
      var1=VAR1; 
      val1=VAL1; 
    else % vary nothing
      var1='tauE';
      val1=tauE;      
    end
    if exist('VAR2','var') && ~isempty(VAR2) && exist('VAL2','var') && ~isempty(VAL2)
      var2=VAR2; 
      val2=VAL2; 
    else % vary nothing
      var2='Estim';
      val2=Estim;  
    end
end

% #########################################################################
nval2=numel(val2);
nval1=numel(val1);
ndt=diff(tspan)/dt/dsfact+1;

% preallocate matrices for analysis results
edges=2:2:60;
muaisirates=zeros(length(edges),nval1*nval2);
cellisirates=zeros(length(edges),nval1*nval2);

% Start simulations
tstart=tic;
cnt=0;
for freq=1:nval2
  % set first variable parameter
  eval(sprintf('%s=%g;',var2,val2(freq)));
  if length(auxval2)==length(val2)
    eval(sprintf('%s=%g;',auxvar2,auxval2(freq)));
  end  
for stim=1:nval1
  % set second variable parameter
  eval(sprintf('%s=%g;',var1,val1(stim)));
  if length(auxval1)==length(val1)
    eval(sprintf('%s=%g;',auxvar1,auxval1(stim)));
  end  
  % update inputs
  ruletask_ResponseLR_CB_inputs;
  if exist('ESTIM','var') && ~isempty(ESTIM)
    Estim=ESTIM;
  end
  if exist('ENOISE','var') && ~isempty(ENOISE)
    Enoise=ENOISE;
  end  
  if exist('INOISE','var') && ~isempty(INOISE)
    Inoise=INOISE;
  end  
  cnt=cnt+1;
  fprintf('processing sim %g of %g\n',cnt,nval1*nval2);
  % ---------------------------------------------------------
  % update derived parameters that might depend on those varied
  Istim=IstimFrac*Estim;
  % ---------------------------------------------------------
  % Specify cell models
  % E soma
  inputparams={'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,'targets',inp_targets};
  EsomaTemplates={};
  for i=1:length(EsomaType)
    EsomaTemplate=[];
    EsomaTemplate.label = 'E';
    EsomaTemplate.multiplicity = nE;
    if ~isempty(InputMatrixE) && EdendType==0
      EsomaTemplate.dynamics = {'V''=InputMatrix(:,floor(1+(t/dt)))+(current)./Cm'};
    else
      EsomaTemplate.dynamics = {'V''=(current)./Cm'};
    end
    switch EsomaType(i)
      case 1 % HH
        EsomaTemplate.mechanisms = {'K','Na','leak',InputMechanism,'randn'};
        EsomaTemplate.parameters = {'Cm',1,'V_IC',-70,'noise',Enoise,...
          'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};
      case 2 % Tallie ACd ci196 (reduced)
        ECa=126.1; EK=-80; ENa=55; Cm=1; gpas=0.017; epas=-66; % passive biophysical parameters
        gnaf=50.25; gkdr=4; gAHP=.014; gM=.084; gcan=.0056; gkca=.2; % active biophysical parameters
        taurca=28.5714;
        %gnaf=gnaf*ones(nE,1);
        %gnaf=gnaf*unifrnd(0,1,[nE 1]);
        EsomaTemplate.mechanisms = {'kdr','naf',InputMechanism,'pas','randn'};
        EsomaTemplate.parameters = {...
          'Cm',Cm,'V_IC',-65,'noise',Enoise,...
          'gnaf',gnaf,'taurca',taurca,'gcan',gcan,'gM',gM,'gkca',gkca,'gkdr',gkdr,'gAHP',gAHP,...
          'gpas',gpas,'epas',epas,'ek',EK,'eca',ECa,'ena',ENa,...
          inputparams{:}};    
      case {'2h'} % Tallie ACd ci196 (reduced + h)
        ECa=126.1; EK=-80; ENa=55; Cm=1; gpas=0.017; epas=-66; % passive biophysical parameters
        gnaf=50.25; gkdr=4; gAHP=.014; gM=.084; gcan=.0056; gkca=.2; % active biophysical parameters
        taurca=28.5714; eh=-10; 
        if ~isequal(var1,'gh') && ~isequal(var2,'gh'), gh=.1; end
        %gnaf=gnaf*ones(nE,1);
        %gnaf=gnaf*unifrnd(0,1,[nE 1]);
        EsomaTemplate.mechanisms = {'kdr','naf','h',InputMechanism,'pas','randn'};
        EsomaTemplate.parameters = {...
          'Cm',Cm,'V_IC',-65,'noise',Enoise,...
          'gnaf',gnaf,'taurca',taurca,'gh',gh,'eh',eh,'gcan',gcan,'gM',gM,'gkca',gkca,'gkdr',gkdr,'gAHP',gAHP,...
          'gpas',gpas,'epas',epas,'ek',EK,'eca',ECa,'ena',ENa,...
          inputparams{:}};    
      case 3 % Wang/Tegner ??
        Ecm=.8; % Enoise=4;
        %EsomaTemplate.mechanisms = {'E_iNa','E_iK','E_iCa','E_iNap','E_CaDyn','E_iCan','E_ileak','randn',InputMechanism};
        EsomaTemplate.mechanisms = {'E_iNa','E_iK','E_ileak','randn',InputMechanism};
        EsomaTemplate.parameters = {'Cm',Ecm,'noise',Enoise,'V_IC',-70,'IC_noise',IC_noise,...
                                    inputparams{:}};
                                    %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};
      case 4 % Tallie ACd ci196
        ECa=126.1; EK=-80; ENa=55; Cm=1; gpas=0.017; epas=-66; % passive biophysical parameters
        gnaf=50.25; gkdr=4; gAHP=.014; gM=.084; gcan=.0056; gkca=.2; % active biophysical parameters
        taurca=28.5714;
        EsomaTemplate.mechanisms = {'cadyn','kdr','AHP','M','naf','kca','can',InputMechanism,'pas','randn'};
        EsomaTemplate.parameters = {...
          'Cm',Cm,'V_IC',-65,'noise',Enoise,...
          'gnaf',gnaf,'taurca',taurca,'gcan',gcan,'gM',gM,'gkca',gkca,'gkdr',gkdr,'gAHP',gAHP,...
          'gpas',gpas,'epas',epas,'ek',EK,'eca',ECa,'ena',ENa,...
          inputparams{:}};
          %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};    
      case 5 % Tallie ACd ci176h w/o wbNa
        if ~isequal(var1,'gh') && ~isequal(var2,'gh'), gh=.01; end
        ECa=126.1; EK=-80; ENa=55; Cm=1;
        gnaf = 65; wbgNa=0; gkdr = 4;
        gAHP = .01; gnap = .0005; gks = 0; taurca = 100; gcat = .001; gcan = 0;
        gkca = .6; taumin = 0; eh = -10; gpas = 0.04; epas = -75; %noise=.15;
        EsomaTemplate.mechanisms = {'cadyn','cat','kdr','h','AHP','iks','kca','naf','nap','pas','randn',InputMechanism};
        EsomaTemplate.parameters = {'c',Cm,'V_IC',-65,'noise',Enoise,...
          'gnaf',gnaf,'gnap',gnap,'taurca',taurca,'gcat',gcat,'gcan',gcan,...
          'gkca',gkca,'taumin',taumin,'gkdr',gkdr,'gks',gks,'gAHP',gAHP,'wbgNa',wbgNa,...
          'gh',gh,'gpas',gpas,'epas',epas,'eh',eh,'ek',EK,'eca',ECa,'ena',ENa,...
          inputparams{:}};
          %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};      
      case 6 % Tallie ACd ci176h (only naf,kdr,h)
        if ~isequal(var1,'gh') && ~isequal(var2,'gh'), gh=.01; end
        ECa=126.1; EK=-80; ENa=55; Cm=1;
        gnaf = 65; wbgNa=0; gkdr = 4;
        gAHP = .01; gnap = .0005; gks = 0; taurca = 100; gcat = .001; gcan = 0;
        gkca = .6; taumin = 0; eh = -10; gpas = 0.04; epas = -75; %noise=.15;
        EsomaTemplate.mechanisms = {'kdr','h','naf','pas','randn',InputMechanism};
        EsomaTemplate.parameters = {'c',Cm,'V_IC',-65,'noise',Enoise,...
          'gnaf',gnaf,'gnap',gnap,'taurca',taurca,'gcat',gcat,'gcan',gcan,...
          'gkca',gkca,'taumin',taumin,'gkdr',gkdr,'gks',gks,'gAHP',gAHP,'wbgNa',wbgNa,...
          'gh',gh,'gpas',gpas,'epas',epas,'eh',eh,'ek',EK,'eca',ECa,'ena',ENa,...
          inputparams{:}};
          %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};      
      case 7 % Papoutsi (all currents, heterogeneous parameters)
        ECa=126.1; EK=-80; ENa=55; Cm=1; 
        EsomaTemplate.mechanisms = {'AHPdist','cadyndist','candist','catdist','hdist','iksdist','kdrdist','kcadist','nafdist','nap','pasdist','randndist',InputMechanism};
        EsomaTemplate.parameters = {...
          'Cm',Cm,'V_IC',-65,'ek',EK,'eca',ECa,'ena',ENa,...
          inputparams{:}};
          %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};    
      case 8 % (Durstewitz-based neuron)
        dE=21; % soma diameter [um]. area=pi*(d/2)^2 dendritic area=4*soma area.
        SE=pi*(dE/2)^2;
        cm    =.8;
        gm    =.54;  Eleak =-80;
        gNa   =80;   ENa   =55;
        gK    =36;   EK    =-80;
        gKs   =25;   EKs   = EK;
        gNap  = .6;  ENap  = ENa;  napshift=5;
        gCan  =.07;  ECan  =-20;   tauca=80; alphacaf=.002;  % gCan=.1
        gCa   =(29.4/SE)*100; ECa=150;
        EsomaTemplate.mechanisms = {'iL','iNa','iK','iCa','iNap','CaDyn','iCan','randn',InputMechanism};
        EsomaTemplate.parameters = {...
          'Cm',Cm,'g_l',gm,'E_l',Eleak,'V_IC',-70,'noise',Enoise,'IC_noise',IC_noise,...
          'gNa',gNa,'ENa',ENa,'gNap',gNap,'napshift',napshift,'ENap',ENap,'gKf',gK,'EKf',EK,'gCaf',gCa,'ECaf',ECa,... 
          'gCan',gCan,'ECan',ECan,'tauca',tauca,'alphacaf',alphacaf,...
          inputparams{:}};
          %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};    
      case 9 % (reduced Durstewitz-based neuron)
        if exist('EnoiseType9','var')
          this_enoise=EnoiseType9;
        else
          this_enoise=Enoise;
        end
        dE=21; % soma diameter [um]. area=pi*(d/2)^2 dendritic area=4*soma area.
        SE=pi*(dE/2)^2;
        cm    =.8;
        gm    =.54;  Eleak =-80;
        gNa   =80;   ENa   =55;
        gK    =36;   EK    =-80;
        gKs   =25;   EKs   = EK;
        gNap  = .6;  ENap  = ENa;  napshift=5;
        gCan  =.07;  ECan  =-20;   tauca=80; alphacaf=.002;  % gCan=.1
        gCa   =(29.4/SE)*100; ECa=150;
        EsomaTemplate.mechanisms = {'iL','iNa','iK','randn',InputMechanism};
        EsomaTemplate.parameters = {...
          'Cm',Cm,'g_l',gm,'E_l',Eleak,'V_IC',-70,'noise',this_enoise,'IC_noise',IC_noise,...
          'gNa',gNa,'ENa',ENa,'gNap',gNap,'napshift',napshift,'ENap',ENap,'gKf',gK,'EKf',EK,'gCaf',gCa,'ECaf',ECa,... 
          'gCan',gCan,'ECan',ECan,'tauca',tauca,'alphacaf',alphacaf,...
          inputparams{:}};
          %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};    
      case 10
        % PY-cell                 L, dgNaf, dgKdr, dgKS, A (iAPoirazi), h, CaN, CaT,
        %                         sAHP, fAHP, dADP. also CaL and CaR, & buffer.      
    end
    EsomaTemplates{i}=EsomaTemplate;
  end
  
  % E dendrite
  EdendTemplate=[];
  EdendTemplate.label = 'Ed';
  EdendTemplate.multiplicity = nE;
  if ~isempty(InputMatrixE)
    EdendTemplate.dynamics = {'V''=InputMatrix(:,floor(1+(t/dt)))+(current)./Cm'};
  else
    EdendTemplate.dynamics = {'V''=(current)./Cm'};
  end
  switch EdendType
    case 1 % Papoutsi PYda
      % ~/models/Papoutsi/PYda1_spec-reduced.mat
      % ~/models/Papoutsi/PYda2_spec-reduced.mat
      gcan = [0.0167];
      gcat = [0.0167];
      gh = [0.0039];
      gks = [0.333];
      gkca = [0.76];
      gkdr = [0.006];
      gnaf = [1.3889];
      gnap = [0.015];
      gpas = [0.047];
      EdendTemplate.mechanisms = {'PYda1_cadyn','PYda1_can','PYda1_cat','PYda1_kca','PYda1_kdr','PYda1_naf','PYda1_nap','PYda1_iks','PYda1_h','PYda1_pas','randn',InputMechanism};
      EdendTemplate.parameters = {'Cm',1.2,'noise',Enoise,'V_IC',-70,'IC_noise',IC_noise,...
                                  inputparams{:},...
                                  'gcan',gcan,'gcat',gcat,'gh',gh,...
                                  'gks',gks,'gkca',gkca,'gkdr',gkdr,...
                                  'gnaf',gnaf,'gnap',gnap,'gpas',gpas};
    case 2 % Papoutsi PYdb
      % ~/models/Papoutsi/PYdb_spec-reduced.mat
      gcan = [0.0167];
      gh = [0.0025];
      gks = [0.1667];
      gkdr = [1.5];
      gnaf = [0.5];
      gnap = [0.005];
      gpas = [0.047];
      EdendTemplate.mechanisms = {'PYdb_cadyn','PYdb_can','PYdb_kdr','PYdb_naf','PYdb_nap','PYdb_iks','PYdb_h','PYdb_pas','randn',InputMechanism};
      EdendTemplate.parameters = {'Cm',1.2,'noise',Enoise,'V_IC',-70,'IC_noise',IC_noise,...
                                  %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
                                  inputparams{:},...
                                  'gcan',gcan,'gh',gh,...
                                  'gks',gks,'gkdr',gkdr,...
                                  'gnaf',gnaf,'gnap',gnap,'gpas',gpas};                                  
    case 3 % (Durstewitz-based neuron)
      cm=2;
      gm    =.54;  Eleak =-80;
      gNa   =80;   ENa   =55;
      gK    =36;   EK    =-80;
      gKs   =25;   EKs   = EK;
      gNap=1.5*.6; ENap  = ENa;  napshift=5;
      EdendTemplate.mechanisms = {'iNa','iKs','iL','iNap','randn',InputMechanism};
      EdendTemplate.parameters = {...
        'Cm',cm,'g_l',gm,'E_l',Eleak,'V_IC',-70,'noise',Enoise,'IC_noise',IC_noise,...
        'gNap',gNap,'napshift',napshift,'ENap',ENap,'gKs',gKs,'EKs',EKs,'gNa',gNa,'ENa',ENa,...
        inputparams{:}};
        %'inputtype',InputType,'onset',onset,'offset',offset,'amp',Estim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};
  end
  % I soma (fast)
  IfastTemplate=[];
  IfastTemplate.label = 'Ifast';
  IfastTemplate.multiplicity = nI;
  if ~isempty(InputMatrixIf)
    IfastTemplate.dynamics = {'V''=InputMatrix(:,floor(1+(t/dt)))+(current)./Cm'};
  else
    IfastTemplate.dynamics = {'V''=(current)./Cm'};
  end
  switch IfastType
    case 1 % HH
      IfastTemplate.mechanisms = {'K','Na','leak',InputMechanism,'randn'};
      IfastTemplate.parameters = {'Cm',1,'V_IC',-70,'noise',Inoise,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};
    case 2 % Wang-Buzsaki interneuron
      Icm=1; % Inoise=2;
      IfastTemplate.mechanisms = {'I_wbNa','I_wbK','I_ileak',InputMechanism,'randn'};
      IfastTemplate.parameters = {'Cm',Icm,'noise',Inoise,'V_IC',-70,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'IC_noise',IC_noise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda};
    case 3 % Kramer08 FS (Basket cell)
      % ~/models/dnsim/Kramer08/ForDave/kramer/B_*
      Cm=.9; ENa=50; E_EKDR=-95; I_EKDR=-100; Eh=-35; IB_Eh=-25; ECa=125; 
      IfastTemplate.mechanisms = {'B_iNaF','B_iKDR','B_leak',InputMechanism,'randn'};
      IfastTemplate.parameters = {...
        'V_IC',-65,'IC_noise',IC_noise,'Cm',Cm,'E_l',-65,'g_l',3,'noise',Inoise,...
        'gNaF',200,'E_NaF',ENa,'NaF_V0',38,'NaF_V1',58.3,'NaF_d1',6.7,'NaF_V2',37,'NaF_d2',15,'NaF_c0',.15,'NaF_c1',1.15,...
        'gKDR',20,'E_KDR',I_EKDR,'KDR_V1',27,'KDR_d1',11.5,'KDR_V2',10,'KDR_d2',10,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        };
    case 4 % Kramer08 LTS
      % ~/models/dnsim/Kramer08/ForDave/kramer/LTS_*
      gAR_L=50; Cm=.9; ENa=50; E_EKDR=-95; I_EKDR=-100; Eh=-35; IB_Eh=-25; ECa=125; IC_noise=0;
      IfastTemplate.mechanisms = {'LTS_iNaF','LTS_iKDR','LTS_iAR','LTS_leak',InputMechanism,'randn'};
      IfastTemplate.parameters = {...
        'V_IC',-65,'IC_noise',IC_noise,'Cm',Cm,'E_l',-65,'g_l',6,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'noise',Inoise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        'gNaF',200,'E_NaF',ENa,'NaF_V0',38,'NaF_V1',58.3,'NaF_d1',6.7,'NaF_V2',37,'NaF_d2',15,'NaF_c0',.15,'NaF_c1',1.15,...
        'gKDR',10,'E_KDR',I_EKDR,'KDR_V1',27,'KDR_d1',11.5,'KDR_V2',10,'KDR_d2',10,...
        'gAR',gAR_L,'E_AR',Eh,'AR_V12',-87.5,'AR_k',-5.5,'c_ARaM',1,'c_ARbM',1,'AR_L',1,'AR_R',1,...
        };
    case 5 % Papoutsi FS (soma)
      % ~/models/Papoutsi/FS-spec.mat
      % PV-cell: L, dgNaf, dgKdr, iAPoirazi, dgKS, CaN, h (ihPoirazi or iAR), & buffer
      gks = [0.0201];
      gkdr = [5];
      gnaf = [62.5];
      gpas = [0.019];
      Cm=1.2;
      IfastTemplate.mechanisms = {'FSs_kdr','FSs_naf','FSs_iks','FSs_pas',InputMechanism,'randn'};
      IfastTemplate.parameters = {'V_IC',-65,'IC_noise',IC_noise,'Cm',Cm,'noise',Inoise,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        'btauoffset',200,'btauscale',3200,...
        'gks',gks,'gkdr',gkdr,'gnaf',gnaf,'gpas',gpas};
    case 6 % Papoutsi CB
      % CB-cell (Poirazi 2014): L, dgNaf, dgKdr, iAPoirazi,       CaT, h (ihPoirazi or iAR), & buffer
      cm=2; gm=.3; Eleak=-70;
      gDR=25; EDR=-80; gA=1; EA=-80; gNa=100; ENa=55;
      gAR=25; E_AR=-35; % ECaf=120; gCaf=0.25; gCan=0.025; ECan=-20; 
      IfastTemplate.mechanisms = {'iL','iNa','dgKdr','iA','iAR',InputMechanism,'randn'};
      IfastTemplate.parameters = {...
        'V_IC',-65,'IC_noise',IC_noise,'Cm',cm,'E_l',Eleak,'g_l',gm,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'noise',Inoise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        'gDR',gDR,'EDR',EDR,'gA',gA,'EA',EA,'gAR',gAR,'E_AR',E_AR,'gNa',gNa,'ENa',ENa};
        % ,'ECaf',ECaf,'gCaf',gCaf,'gCan',gCan,'ECan',ECan      
        % TODO: compare iNa and iAR b/w ~/code/scc_code... and /base/mechanisms
  end
  % I soma (slow)
  IslowTemplate=[];
  IslowTemplate.label = 'Islow';
  IslowTemplate.multiplicity = nIs;
  if ~isempty(InputMatrixIs)
    IslowTemplate.dynamics = {'V''=InputMatrix(:,floor(1+(t/dt)))+(current)./Cm'};
  else
    IslowTemplate.dynamics = {'V''=(current)./Cm'};
  end
  switch IslowType
    case {1,6} % Papoutsi CB
      % CB-cell (Poirazi 2014): L, dgNaf, dgKdr, iAPoirazi,       CaT, h (ihPoirazi or iAR), & buffer
      cm=2; gm=.3; Eleak=-70;
      gDR=25; EDR=-80; gA=1; EA=-80; gNa=100; ENa=55;
      gAR=25; E_AR=-35; % ECaf=120; gCaf=0.25; gCan=0.025; ECan=-20; 
      IslowTemplate.mechanisms = {'iL','iNa','dgKdr','iA','iAR',InputMechanism,'randn'};
      IslowTemplate.parameters = {...
        'V_IC',-65,'IC_noise',IC_noise,'Cm',cm,'E_l',Eleak,'g_l',gm,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'noise',Isnoise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        'gDR',gDR,'EDR',EDR,'gA',gA,'EA',EA,'gAR',gAR,'E_AR',E_AR,'gNa',gNa,'ENa',ENa};
    case 5 % Papoutsi FS (soma)
      % ~/models/Papoutsi/FS-spec.mat
      % PV-cell: L, dgNaf, dgKdr, iAPoirazi, dgKS, CaN, h (ihPoirazi or iAR), & buffer
      gks = [0.0201];
      gkdr = [5];
      gnaf = [62.5];
      gpas = [0.019];
      Cm=1.2;
      IslowTemplate.mechanisms = {'FSs_kdr','FSs_naf','FSs_iks','FSs_pas',InputMechanism,'randn'};
      IslowTemplate.parameters = {'V_IC',-65,'IC_noise',IC_noise,'Cm',Cm,'noise',Isnoise,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        'btauoffset',200,'btauscale',3200,...
        'gks',gks,'gkdr',gkdr,'gnaf',gnaf,'gpas',gpas};
    case 7 % based on Kramer08 LTS, modified for DLPFC CB+ IN (gKDR 10->5)
      % ~/models/dnsim/Kramer08/ForDave/kramer/LTS_*
      gAR_L=50; Cm=.9; ENa=50; E_EKDR=-95; I_EKDR=-100; Eh=-35; IB_Eh=-25; ECa=125; IC_noise=0;
      gKDR=5; % decreased from 10 to achieve longer AP
      % to achieve spiking: set Istim>=21.871
        % -------------------------------------------------
      % MODIFICATION NOTES: how to make it a "PFC CB+" cell
        % -------------------------------------------------
      % reduced gKDR 10->5 to achieve longer AP observed in DLPFC CB+ INs (Zaitsev)
      % RMP: -64mV (Kawaguchi; same as Papoutsi CB RS)
      % threshold: -52.5mV @ Istim>21.87 (a little lower than Papoutsi CB RS-51mV)
      IslowTemplate.mechanisms = {'LTS_iNaF','LTS_iKDR','LTS_iAR','LTS_leak',InputMechanism,'randn'};
      IslowTemplate.parameters = {...
        'V_IC',-65,'IC_noise',IC_noise,'Cm',Cm,'E_l',-65,'g_l',6,...
        'inputtype',InputType,'onset',onset,'offset',offset,'amp',Istim,'noise',Isnoise,'sharedfraction',sharedfraction,'Ninputs',Ninputs,'tauD',inp_tauD,'tauR',inp_tauR,'lambda',inp_lambda,...
        'gNaF',200,'E_NaF',ENa,'NaF_V0',38,'NaF_V1',58.3,'NaF_d1',6.7,'NaF_V2',37,'NaF_d2',15,'NaF_c0',.15,'NaF_c1',1.15,...
        'gKDR',gKDR,'E_KDR',I_EKDR,'KDR_V1',27,'KDR_d1',11.5,'KDR_V2',10,'KDR_d2',10,...
        'gAR',gAR_L,'E_AR',Eh,'AR_V12',-87.5,'AR_k',-5.5,'c_ARaM',1,'c_ARbM',1,'AR_L',1,'AR_R',1,...
        };      
  end  
  
  if include_thalamus_flag
    
    % basic: awake = depolarized, asleep = hyperpolarized
    % in thalamus, hyperpolarization causes weird oscillations (due to h and T-currents)
    
    % to make an "awake" thalamic model (source: Austin Soplata discussion):
    % change the RMP (using potassium leak and h-current)
    %   iH_TC.gH (default sleep: 0.005; awake: try .02-.1), >.01, up to .06 or .08)
    %   iKLeak_TC.gKLeak (default sleep: 0.0172; awake: try >.0172)
    %   note: the effects of gH and gKLeak depend on each other
    % RMP increases with gH, keep RMP <= -50mV or -55mV, and > -80mV
      % avg sleep RMP: -72mV (down to -80mV)
      % avg awake RMP: more depolarized (-55mV to -65mV or -70mV)
      %   - occassional bursting begins for RMP ~ -70mV, the increases below)
    % awake thalamus has non-bursting TC (i.e., in relay mode)
    
    % Thalamocortical cell
    TC_Template=[];
    TC_Template.label = 'TC';
    TC_Template.multiplicity = nTC;
    if ~isempty(InputMatrixTC)
      TC_Template.dynamics = {'V''=InputMatrix(:,floor(1+(t/dt)))+(current)./Cm'};
    else
      TC_Template.dynamics = {'V''=(current)./Cm'};
    end
    switch TC_Type
      case 1 % Austin's model with parameters from Destexhe 1996
        % ref: Destexhe, A., Bal, T., McCormick, D. A., & Sejnowski, T. J. (1996). Ionic mechanisms underlying synchronized oscillations and propagating waves in a model of ferret thalamic slices. Journal of Neurophysiology, 76(3), 2049–2070.
        % email: https://mail.google.com/mail/u/0/#search/thalamic+model+from%3Aaustin/14ca06949c38ec2d
        gLeak=.03;%.01; 
        gK=10;
        % --------------------------------
        % How to "awaken" the thalamus:
        % --------------------------------
        % Experimental targets:
        % 1. awake RMP: -55mV to -65mV
        % 2. awake TC should not burst
        % Parameter changes:
        % 1. try iH_TC.gH = .02 or .08
        % 2. try iKLeak_TC.gKLeak > .0172
        % --------------------------------
        if include_thalamus_flag==1 % asleep
          fprintf('PFC coupled to "asleep" thalamus.\n');
          % asleep:
          gH=0.005;       % mS/cm^2
          gKLeak=0.0172;  % mS/cm^2 
        elseif include_thalamus_flag==2 % awake
          fprintf('PFC coupled to "awake" thalamus.\n');
          % awake: (needs more work; experimental targets not met)
          gH=0.08;        % .07
          gKLeak=0.025;   % .025
        end
        % --------------------------------
        TC_Template.mechanisms = {'Ca_TC','iH_TC','iK_TC','iKLeak_TC','iLeak_TC','iNa_TC','iT_TC'};
        TC_Template.parameters = {'Cm',1,'V_IC',-65,'gLeak',gLeak,'gK',gK,'gH',gH,'gKLeak',gKLeak};        
    end
    % TRN/RE inhibitory cell
    RE_Template=[];
    RE_Template.label = 'RE';
    RE_Template.multiplicity = nRE;
    if ~isempty(InputMatrixRE)
      RE_Template.dynamics = {'V''=InputMatrix(:,floor(1+(t/dt)))+(current)./Cm'};
    else
      RE_Template.dynamics = {'V''=(current)./Cm'};
    end
    switch RE_Type
      case 1 % Austin's model with parameters from Destexhe 1996
        % ref: Destexhe, A., Bal, T., McCormick, D. A., & Sejnowski, T. J. (1996). Ionic mechanisms underlying synchronized oscillations and propagating waves in a model of ferret thalamic slices. Journal of Neurophysiology, 76(3), 2049–2070.
        % email: https://mail.google.com/mail/u/0/#search/thalamic+model+from%3Aaustin/14ca06949c38ec2d
        RE_Template.mechanisms = {'iK_RE','iLeak_RE','iNa_RE','iT_RE'};
        RE_Template.parameters = {'Cm',1,'V_IC',-85};
    end
  end
  
  % create multi-layer specification
  spec=[];
  for i=1:numLayers
    E=EsomaID(i);
    if length(EsomaTemplates)>=i
      EsomaTemplate=EsomaTemplates{i};
    end
    % specify nodes in this layer
    spec.nodes(E) = EsomaTemplate;
    spec.nodes(E).label = sprintf('%s%g',EsomaTemplate.label,i);
    if i>1
      spec.nodes(E).multiplicity=nE2;
      this_geeampa=geeampa2;
      this_geenmda=geenmda2;
    else
      this_geeampa=geeampa;
      this_geenmda=geenmda;              
    end
    if ~isempty(InputMatrixE) && EdendType==0
      spec.nodes(E).parameters={spec.nodes(E).parameters{:},'InputMatrix',InputMatrixE{i}};
    end
    % ---------------------------------------------------------------------
    % optional nodes:
    % ---------------------------------------------------------------------
    if EdendType>0
      Ed=EdendID(i);
      spec.nodes(Ed) = EdendTemplate;
      spec.nodes(Ed).label = sprintf('%s%g',EdendTemplate.label,i);
      if i>1
        % update multiplicity
        spec.nodes(Ed).multiplicity=nE2;
        % remove inputs from non-input layers
        spec.nodes(Ed).mechanisms = setdiff(spec.nodes(Ed).mechanisms,InputMechanism);
      end
      if ~isempty(InputMatrixE)
        spec.nodes(Ed).parameters={spec.nodes(Ed).parameters{:},'InputMatrix',InputMatrixE{i}};
      end
      % specify inter-compartmental connections within this layer
      spec.connections(E,Ed).label = [spec.nodes(E).label '-' spec.nodes(Ed).label];
      spec.connections(E,Ed).mechanisms = {'iCOM','AMPA','NMDA'};
      spec.connections(E,Ed).parameters = {'gcore',gEEd,'tauDx',tauE,'g_SYN',this_geeampa,'g_NMDA',this_geenmda,'NtauD',tauNMDA,'NtauR',tauNMDArise,'normalize',normalize,'PrSYN',PrSYNee,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',zerodiag,'mask',EEmask{i}};
      spec.connections(Ed,E).label = [spec.nodes(Ed).label '-' spec.nodes(E).label];
      spec.connections(Ed,E).mechanisms = {'iCOM'};
      spec.connections(Ed,E).parameters = {'gcore',gEdE};
    else
      spec.connections(E,E).label = [spec.nodes(E).label '-' spec.nodes(E).label];
      spec.connections(E,E).mechanisms = {'AMPA','NMDA'};
      spec.connections(E,E).parameters = {'tauDx',tauE,'g_SYN',this_geeampa,'g_NMDA',this_geenmda,'NtauD',tauNMDA,'NtauR',tauNMDArise,'normalize',normalize,'PrSYN',PrSYNee,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',zerodiag,'mask',EEmask{i}};
    end
    % ---------------------------------------------------------------------
    if IfastType>0 && (i==1 || nI2>0)
      Ifast=IfastID(i);
      spec.nodes(Ifast) = IfastTemplate;
      spec.nodes(Ifast).label = sprintf('%s%g',IfastTemplate.label,i);
      if i>1
        spec.nodes(Ifast).multiplicity = nI2;
        taui=tauI2;
        gei=gEI2;
        gie=gIE2;
        gii=gII2;
        geinmda=gEInmda2;
        %gie=gIE*(nI/nI2); % nE/nE2 nI/nI2 nIs/nIs2
        %gei=gEI*(nE/nE2);
      else
        taui=tauI;
        gie=gIE;
        gei=gEI;
        gii=gII;
        geinmda=gEInmda;
      end
      if ~isempty(InputMatrixIf)
        spec.nodes(Ifast).parameters={spec.nodes(Ifast).parameters{:},'InputMatrix',InputMatrixIf{i}};
      end
      % connections
      spec.connections(E,Ifast).label = [spec.nodes(E).label '-' spec.nodes(Ifast).label];
      spec.connections(E,Ifast).mechanisms = {'AMPA','NMDA'};
      spec.connections(E,Ifast).parameters = {'tauDx',tauE,'g_SYN',gei,'g_NMDA',geinmda,'NtauD',tauNMDA,'NtauR',tauNMDArise,'normalize',normalize,'PrSYN',PrSYNei,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',EIfmask{i}};
      spec.connections(Ifast,E).label = [spec.nodes(Ifast).label '-' spec.nodes(E).label];
      spec.connections(Ifast,E).mechanisms = {'GABAa'};
      spec.connections(Ifast,E).parameters = {'tauDx',taui,'g_SYN',gie,'E_SYN',Egaba,'normalize',normalize,'PrSYN',PrSYNie,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',IfEmask{i}};
      spec.connections(Ifast,Ifast).label = [spec.nodes(Ifast).label '-' spec.nodes(Ifast).label];
      spec.connections(Ifast,Ifast).mechanisms = {'GABAa'};
      spec.connections(Ifast,Ifast).parameters = {'tauDx',taui,'g_SYN',gii,'E_SYN',Egaba,'normalize',normalize,'PrSYN',PrSYNii,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',zerodiag,'mask',IfIfmask{i}};
    end
    % ---------------------------------------------------------------------
    if ~isequal(IslowType,0) && (i==1 || nIs2>0)
      Islow=IslowID(i);
      spec.nodes(Islow) = IslowTemplate;
      spec.nodes(Islow).label = sprintf('%s%g',IslowTemplate.label,i);      
      if i>1
        spec.nodes(Islow).multiplicity = nIs2;
        tauis=tauIs2;
        gise=gIsE2;
        geis=gEIs2;
        gisis=gIsIs2;
        gisi=gIsI2;
      else
        tauis=tauIs;
        gise=gIsE;
        geis=gEIs;
        gisis=gIsIs;
        gisi=gIsI;
      end
      if ~isempty(InputMatrixIs)
        spec.nodes(Islow).parameters={spec.nodes(Islow).parameters{:},'InputMatrix',InputMatrixIs{i}};
      end
      % connections
      spec.connections(E,Islow).label = [spec.nodes(E).label '-' spec.nodes(Islow).label];
      spec.connections(E,Islow).mechanisms = {'AMPA'};
      spec.connections(E,Islow).parameters = {'tauDx',tauE,'g_SYN',geis,'normalize',normalize,'PrSYN',PrSYNei,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',EIsmask{i}};
      if EdendType>0
        Etarget=Ed;
      else
        Etarget=E;
      end
      spec.connections(Islow,Etarget).label = [spec.nodes(Islow).label '-' spec.nodes(Etarget).label];
      spec.connections(Islow,Etarget).mechanisms = {'GABAa'};
      spec.connections(Islow,Etarget).parameters = {'tauDx',tauis,'g_SYN',gise,'E_SYN',Egaba,'normalize',normalize,'PrSYN',PrSYNie,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',IsEmask{i}};
      spec.connections(Islow,Islow).label = [spec.nodes(Islow).label '-' spec.nodes(Islow).label];
      spec.connections(Islow,Islow).mechanisms = {'GABAa'};
      spec.connections(Islow,Islow).parameters = {'tauDx',tauis,'g_SYN',gisis,'E_SYN',Egaba,'normalize',normalize,'PrSYN',PrSYNii,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',zerodiag,'mask',IsIsmask{i}};      
      if IfastType>0
        spec.connections(Islow,Ifast).label = [spec.nodes(Islow).label '-' spec.nodes(Ifast).label];
        spec.connections(Islow,Ifast).mechanisms = {'GABAa'};
        spec.connections(Islow,Ifast).parameters = {'tauDx',tauis,'g_SYN',gisi,'E_SYN',Egaba,'normalize',normalize,'PrSYN',PrSYNii,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',IsIfmask{i}};      
        spec.connections(Ifast,Islow).label = [spec.nodes(Islow).label '-' spec.nodes(Islow).label];
        spec.connections(Ifast,Islow).mechanisms = {'GABAa'};
        spec.connections(Ifast,Islow).parameters = {'tauDx',taui,'g_SYN',gIIs,'E_SYN',Egaba,'normalize',normalize,'PrSYN',PrSYNii,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',IfIsmask{i}};      
      end
    end
    % ---------------------------------------------------------------------
    % remove inputs from non-input layers
    if IfastType>0
      spec.nodes(Ifast).mechanisms = setdiff(spec.nodes(Ifast).mechanisms,InputMechanism);
    end
    if ~isequal(IslowType,0)
      spec.nodes(Islow).mechanisms = setdiff(spec.nodes(Islow).mechanisms,InputMechanism);
    end
    if EdendType>0 % remove from soma if dendrite is present
      spec.nodes(E).mechanisms = setdiff(spec.nodes(E).mechanisms,InputMechanism);
      if i>1 % remove from dendrite in downstream layers
        spec.nodes(Ed).mechanisms = setdiff(spec.nodes(Ed).mechanisms,InputMechanism);
      end
    elseif i>1 % remove from soma in downstream layers if dendrite not present
      spec.nodes(E).mechanisms = setdiff(spec.nodes(E).mechanisms,InputMechanism);
    end
    % ---------------------------------------------------------------------
    % specify connections from previous layer to this layer
    if i>1
      Elast=EsomaID(i-1); 
      if EdendType>0 % connect to dendrite
        Etarget=Ed;
      else % connect to soma
        Etarget=E;
      end
      spec.connections(Elast,Etarget).label = [spec.nodes(Elast).label '-' spec.nodes(Etarget).label];
      spec.connections(Elast,Etarget).mechanisms = {'AMPA','NMDA'};
      spec.connections(Elast,Etarget).parameters = {'tauDx',tauE,'g_SYN',gEEampa,'g_NMDA',gEEnmda,'NtauD',tauNMDA,'NtauR',tauNMDArise,'normalize',normalize,'PrSYN',PrSYNEE,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',EEmaskFF{i-1}};
    end
  end
  if include_thalamus_flag
    spec.nodes(TC_ID) = TC_Template;
    if ~isempty(InputMatrixTC)
      spec.nodes(TC_ID).parameters={spec.nodes(TC_ID).parameters{:},'InputMatrix',InputMatrixTC};
    end
    spec.nodes(RE_ID) = RE_Template;    
    if ~isempty(InputMatrixRE)
      spec.nodes(RE_ID).parameters={spec.nodes(RE_ID).parameters{:},'InputMatrix',InputMatrixRE};
    end
    % cortico-thalamic projection
    src=EsomaID(end);
    dst=TC_ID;
    spec.connections(src,dst).label = [spec.nodes(src).label '-' spec.nodes(dst).label];
    spec.connections(src,dst).mechanisms = {'AMPA'};
    spec.connections(src,dst).parameters = {'tauDx',tauE,'g_SYN',gETC,'normalize',normalize,'PrSYN',1,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',ETCmask};
    % thalamocortical projection
    src=TC_ID;
    dst=EsomaID(1);
    spec.connections(src,dst).label = [spec.nodes(src).label '-' spec.nodes(dst).label];
    spec.connections(src,dst).mechanisms = {'AMPA'};
    spec.connections(src,dst).parameters = {'tauDx',tauE,'g_SYN',gTCE,'normalize',normalize,'PrSYN',1,'IC_noise',IC_noise,'ConnType',ConnType,'span',span,'zerodiag',0,'mask',TCEmask};
    % TC<->RE
    spec.connections(TC_ID,RE_ID).label = [spec.nodes(TC_ID).label '-' spec.nodes(RE_ID).label];%'TC-RE';
    spec.connections(TC_ID,RE_ID).mechanisms = {'iAMPA'};
    spec.connections(TC_ID,RE_ID).parameters = {'gAMPA',gTCRE,'netcon',TCREmask};
    spec.connections(RE_ID,TC_ID).label = [spec.nodes(RE_ID).label '-' spec.nodes(TC_ID).label];%'RE-TC';
    spec.connections(RE_ID,TC_ID).mechanisms = {'iGABAA','iGABAB'};
    spec.connections(RE_ID,TC_ID).parameters = {'gGABAA_base',gRETCa,'gGABAB',gRETCb,'netcon',RETCmask};
    spec.connections(RE_ID,RE_ID).label = [spec.nodes(RE_ID).label '-' spec.nodes(RE_ID).label];%'RE-RE';
    spec.connections(RE_ID,RE_ID).mechanisms = {'iGABAA'};
    spec.connections(RE_ID,RE_ID).parameters = {'gGABAA_base',gRERE,'netcon',REREmask};
    % spec.connections(RE_ID,RE_ID).mechanisms = {'iAMPA_PY_faux','iGABAA'};    
  end

  % process specification and simulate model
  if buildonly==1
    tmp=buildmodel(spec);
    break;
  end
  [data,model] = runsim(spec,'timelimits',tspan,'dt',dt,'dsfact',dsfact,'SOLVER',solver,'coder',usecoder,'debug',1,'verbose',verbose);
  if skip_plotv==0 && (plot_series_flag || (nval1==1 && nval2==1))
    plotv(data,spec,'varlabel','V','xlim',xlims/1000);
    set(gcf,'position',[77 500 1000 450]);
  end
  % plotpow(data,spec);
  % dnsim(spec);

  % preallocate matrices for analysis results
  if cnt==1
    Epops={spec.nodes(EsomaID).label};
    Ipops={spec.nodes(IfastID).label};
    npops=length(spec.nodes);
    poplabels={spec.nodes.label};
    muaOscFreqs = zeros(npops,nval1,nval2); % 2pops
    muaAreaPowers = zeros(npops,nval1,nval2);
    suaOscFreqs = nan(npops,nE,nval1,nval2); % 2pops
      % to capture STO freqs, set: gIE=0, Estim=.9, Enoise=0, InputType=1
    suaAreaPowers = nan(npops,nE,nval1,nval2);
    suaOscFreqsAvg = nan(npops,nval1,nval2); % 2pops
    suaAreaPowersAvg = nan(npops,nval1,nval2);    
    rateOscFreqs = zeros(npops,nval1,nval2);
    rateAreaPowers = zeros(npops,nval1,nval2);
    AvgSpikeRates = zeros(npops,nval1,nval2); % E;I
    NumSpikes = zeros(npops,nval1,nval2);
    MUAs = zeros(ndt,npops,nval1,nval2);
    LFPs = zeros(ndt,npops,nval1,nval2);
    spkcoh = zeros(npops,npops,nval1,nval2); %  % E1-E1, E2-E2, E1-E2
  end
  
  % #######################################################################
  % ANALYZE ALL POPULATIONS
  % -----------------------------------------------------------------------
  % Parameters
  Fmin=2;
  powthreshprc=95; PxxSmooth=5;
  Fmax=150; Fwin=5; NOVERLAP=[]; % spectral parameters
  t=data(1).epochs.time;
  selE=1:nE; % E subset in MUA and LFP
  selI=1:nI; % I subset in MUA and LFP
  
  if InputType==1
    window_size=t(end)/10;
  else
    window_size=t(end)/(20*nval2);
  end
  dW=window_size/10; spikethreshold=0; smoothfact=1;
  
  % calculate spike rates
  [h,rates,tmins,spiketimes,spikeinds]=plotspk(data,spec,'plot_flag',0,'window_size',window_size,'dW',dW,'spikethreshold',spikethreshold); % firing rate(t) and FRH
  if cnt==1
    allrates = zeros(length(tmins),npops,nval1,nval2);
    allspiketimes = cell(npops,nval1,nval2);
  end
  allspiketimes(:,stim,freq)=spiketimes;
  % -----------------------------------------------------------------------
  % Per population analysis
  % -----------------------------------------------------------------------
  for a=1:npops
    % cell selection
    popname=spec.nodes(a).label;
    if ismember(popname,Epops)
      if a==1
        cellsel=1:nE;
      else
        cellsel=1:nE2;
      end
    elseif ismember(popname,Ipops)
      if isequal(popname,'Ifast1') %a==1
        cellsel=1:nI;
      else
        cellsel=1:nI2;
      end
    else
      cellsel=1:spec.nodes(a).multiplicity;
    end
    % extract data
    Vind=find(strcmp([popname '_V'],{data(a).sensor_info.label}));
    Vdat=squeeze(data(a).epochs.data(Vind,:,cellsel));
    if length(cellsel)==1, Vdat=Vdat'; end
    Vavg=double(nanmean(data(a).epochs.data(Vind,:,cellsel),3))';
    ravg=mean(rates{a}(cellsel,:),1);
    ravg=smooth(ravg,smoothfact);

    MUAs(:,a,stim,freq)=Vavg;
    allrates(:,a,stim,freq)=ravg;
    t1=nearest(tmins,TOILIM(1)/1000);
    t2=nearest(tmins,TOILIM(2)/1000);
    AvgSpikeRates(a,stim,freq)=mean(allrates(t1:t2,a,stim,freq));
    %AvgSpikeRates(a,stim,freq)=mean(cellfun(@length,spiketimes{a}))/(t(end)-t(1));
    NumSpikes(a,stim,freq)=sum(cellfun(@length,spiketimes{a}))/length(spiketimes{a});
    
    % MUA spectrum
    t1=nearest(t,TOILIM(1)/1000);
    t2=nearest(t,TOILIM(2)/1000);    
    X=Vavg(t1:t2);
    %X=Vavg;
    Fs = fix(1/(t(2)-t(1)));
    NFFT=2^(nextpow2(length(X)-1)-2);
    WINDOW=2^(nextpow2(NFFT-1)-3);
    FreqRange=[max(Fmin,2/t(end)) Fmax]; % frequencies to consider
    X=smooth(X,ceil(.1/dt));
    [Pxx,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    tmpPxx=smooth(Pxx,PxxSmooth);
    ht=prctile(log10(tmpPxx(sel)),powthreshprc);
    [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
    if ~isempty(PPind)  
      PPind=PPind(1);
      OscFreq = f(sel(PPind));
      flo=OscFreq-Fwin/2;
      fhi=OscFreq+Fwin/2;
      sel2=find(flo<=f & f<=fhi);
      AreaPower = sum(Pxx(sel2))*(f(2)-f(1));
      muaOscFreqs(a,stim,freq)=OscFreq;
      muaAreaPowers(a,stim,freq)=AreaPower;
    else
      muaOscFreqs(a,stim,freq)=nan;
      muaAreaPowers(a,stim,freq)=nan;        
    end
    if cnt==1 && a==1
      MUApows = zeros(length(Pxx),npops,nval1,nval2);
    end
    MUApows(:,a,stim,freq)=Pxx;
    muaf=f;

    % Average SUA spectrum
    t1=nearest(t,TOILIM(1)/1000);
    t2=nearest(t,TOILIM(2)/1000);    
    Fs = fix(1/(t(2)-t(1)));
    FreqRange=[max(Fmin,2/t(end)) Fmax]; % frequencies to consider
    Pxx=0;
    for i=1:size(Vdat,2)
      X=Vdat(t1:t2,i);
      NFFT=2^(nextpow2(length(X)-1)-1);%2); % <-- use higher resolution to capture STO freq variation
      WINDOW=2^(nextpow2(NFFT-1)-3);
      %X=smooth(X,ceil(.1/dt));
      [tmpPxx,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);
      if all(isnan(tmpPxx(:)))
        tmpPxx=zeros(size(tmpPxx));
      end
      tmpPxx=double(tmpPxx);
      %tmpPxx=smooth(Pxx,PxxSmooth);
      sel = find(FreqRange(1)<=f & f<=FreqRange(end));
      ht=prctile(log10(tmpPxx(sel)),powthreshprc);
      [PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
      if ~isempty(PPind)  
        PPind=PPind(1);
        OscFreq = f(sel(PPind));
        flo=OscFreq-Fwin/2;
        fhi=OscFreq+Fwin/2;
        sel2=find(flo<=f & f<=fhi);
        AreaPower = sum(tmpPxx(sel2))*(f(2)-f(1));
        suaOscFreqs(a,i,stim,freq)=OscFreq;
        suaAreaPowers(a,i,stim,freq)=AreaPower;
      else
        suaOscFreqs(a,i,stim,freq)=nan;
        suaAreaPowers(a,i,stim,freq)=nan;        
      end
      if i==1
        Pxx=tmpPxx;
      else
        Pxx=Pxx+tmpPxx/size(Vdat,2);
      end
    end
    if cnt==1 && a==1
      avgSUApows = zeros(length(Pxx),npops,nval1,nval2);
    end
    avgSUApows(:,a,stim,freq)=Pxx;
    suaOscFreqsAvg(a,stim,freq)=nanmean(suaOscFreqs(a,:,stim,freq));
    suaAreaPowersAvg(a,stim,freq)=nanmean(suaAreaPowers(a,:,stim,freq));
    suaf=f;
    
    % Spike rate spectrum
    X=ravg;
    Fs = fix(1/(tmins(2)-tmins(1)));
    NFFT=2^(nextpow2(length(tmins)-1)-1);
    WINDOW=2^(nextpow2(NFFT-1)-3);
    FreqRange=[max(Fmin,2/tmins(end)) Fmax]; % frequencies to consider
    X=smooth(X,ceil(1/(1/Fs)));
    [Pxx,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    ht=prctile(log10(Pxx(sel)),90);
    [PeakPower,PPind]=findpeaks(log10(Pxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
    if ~isempty(PPind)
      PPind=PPind(PeakPower==max(PeakPower));
      OscFreq = f(sel(PPind));
      flo=OscFreq-Fwin/2;
      fhi=OscFreq+Fwin/2;
      sel2=find(flo<=f & f<=fhi);
      AreaPower = sum(Pxx(sel2))*(f(2)-f(1));
      rateOscFreqs(a,stim,freq)=OscFreq;
      rateAreaPowers(a,stim,freq)=AreaPower;
    else
      rateOscFreqs(a,stim,freq)=nan;
      rateAreaPowers(a,stim,freq)=nan;    
    end
    if cnt==1 && a==1
      ratepows = zeros(length(Pxx),npops,nval1,nval2);
    end
    ratepows(:,a,stim,freq)=Pxx;
    ratef=f;  
  end
  
  % -----------------------------------------------------------------------
  % Pairwise population analysis
  % -----------------------------------------------------------------------
  
  % Spike coherence
  % method: convolve spike train with EPSP then calc spike-spike coherence
  % note: decrease 'tau' to achieve stricter coherence estimate

%   tau = 10;
%   psp = exp(-tmins/(tau/1000));
%   psp = [zeros(1,length(psp)) psp];  
  
  t1=nearest(t,TOILIM(1)/1000);
  t2=nearest(t,TOILIM(2)/1000);    
  tt=t; % tmins
  taud=5;%5;%10
  taur=.9*taud;
  psp = (exp(-tt/(taud/1000))-exp(-tt/(taur/1000)));
  % psp = gaussian
  psp = psp/max(psp);
  psp = [zeros(1,length(psp)) psp];
% 
%   taud=5;%10;
%   taur=.9*taud;  
%   psp = (exp(-tmins/(taud/1000))-exp(-tmins/(taur/1000)));
%   psp = psp/max(psp);
%   psp = [zeros(1,length(psp)) psp];
  
  for a=1:npops
    E1=find(strcmp(poplabels{a},{spec.nodes.label}));
    n1=spec.nodes(E1).multiplicity;
    for b=a:npops
      E2=find(strcmp(poplabels{b},{spec.nodes.label}));
      n2=spec.nodes(E2).multiplicity;
      coh=nan(n1,n2);
      if b==a, lfp=zeros(length(tt),1); end
      for i=1:n1
        spks=spiketimes{E1}{i};
        spks=spks(spks>=TOILIM(1)/1000&spks<=TOILIM(2)/1000);
        x1=histc(spks,tt); if size(x1,1)==1, x1=x1'; end
        x1=conv(x1,psp,'same');
        if b==a, lfp=lfp+x1; end
        for j=1:n2
          spks=spiketimes{E2}{j};
          spks=spks(spks>=TOILIM(1)/1000&spks<=TOILIM(2)/1000);
          x2=histc(spks,tt); if size(x2,1)==1, x2=x2'; end
          x2=conv(x2,psp,'same');
          coh(i,j)=sum(x1.*x2)./sqrt(sum(x1.^2).*sum(x2.^2));
        end
      end
      if b==a, LFPs(:,a,stim,freq)=lfp; end
      avgcoh=nanmean(coh(:));
      spkcoh(a,b,stim,freq)=avgcoh;
      fprintf('%s-%s spkcoh: %g\n',poplabels{a},poplabels{b},avgcoh);
      if plot_coh_flag %plot_series_flag
        if a==1 && b==1, figure('position',[730 15 950 850]); end
        subplot(npops,npops,b+(a-1)*npops); imagesc(coh); caxis([0 1]);
        ylabel(poplabels{a}); xlabel(poplabels{b}); if a==1, title('spkcoh'); end
      end
      
      if a==1 && b==1
        % look for assemblies, defined wrt pairwise spike coherence
        prcthresh=80;
        cutoff=prctile(coh(:),prcthresh); % .6
        bigcoh=zeros(size(coh));
        bigcoh(coh>cutoff)=coh(coh>cutoff);
        symi=symrcm(bigcoh);
        if 0
          figure; 
          subplot(2,2,1); imagesc(bigcoh); axis square; xlabel('spkcoh'); ylabel('spkcoh');
          title(sprintf('cutoff=%3.3g (%g%%-tile)',cutoff,prcthresh));
          subplot(2,2,3); imagesc(bigcoh(symi,symi)); axis square; xlabel('spkcoh'); ylabel('spkcoh');
          title(sprintf('cutoff=%3.3g (%g%%-tile), symrcm-sorted',cutoff,prcthresh));
          popname=spec.nodes(a).label;
          Vind=find(strcmp([popname '_V'],{data(a).sensor_info.label}));
          Vdat=squeeze(data(a).epochs.data(Vind,:,1:size(coh,1)));
          subplot(2,2,2);
          imagesc(t,1:size(Vdat,2),Vdat'); axis xy; colormap(1-gray); %colorbar
          title(sprintf('%s V(t)',popname)); ylabel(popname); xlabel('t');
          subplot(2,2,4);
          imagesc(t,1:size(Vdat,2),Vdat(:,symi)'); axis xy; colormap(1-gray); %colorbar
          title(sprintf('%s V(t)',popname)); ylabel([popname ' (symrcm-sorted)']); xlabel('t');
          
          figure; plot(tm0can,symi); xlabel('tm0can'); ylabel('symi');
        end
      end
      
    end
  end

  % MUA synchrony
  % ...
  
  % #######################################################################  
  % VEind->Vind, selE->cellsel, Eind->a, VE->Vavg, rE->ravg

  % what to plot in series? (1st two populations)
  Eind=EsomaID(1); 
  Iind=IfastID(1);
  VE=MUAs(:,Eind,stim,freq);
  VE2=MUAs(:,EsomaID(end),stim,freq);
  VI=MUAs(:,Iind,stim,freq);
  rE=allrates(:,Eind,stim,freq);
  rE2=allrates(:,EsomaID(end),stim,freq);
  rI=allrates(:,Iind,stim,freq);
  Ename=poplabels{Eind};
  Ename2=poplabels{EsomaID(end)};
  Iname=poplabels{Iind};
  ne=length(spiketimes{Eind});
  
  % MUA-based spectral peaks
  tmp=smooth(VE,3);
  [pkval,pkind]=findpeaks(tmp,'MinPeakHeight',prctile(tmp,90));    
  % (ISI to next cycle)
  tcycle=t(pkind);
  nextISI=zeros(ne,3); % (mean,stdev,n)
  allisis=[];
  for i=1:ne
    spkt=spiketimes{Eind}{i};
    isis=zeros(length(spkt),1);
    for j=1:length(spkt)
      this=nearest(tcycle,spkt(j));
      if this>=numel(tcycle), break; end
      next=this+1;
      isis(j)=tcycle(next)-tcycle(this);
    end
    isis(isis==0)=[];
    nextISI(i,:)=[mean(isis) std(isis) numel(isis)];
    if ~isempty(isis), allisis=cat(1,allisis,isis); end
  end
  isi=diff(t(pkind));
  [n,bin]=histc(round(1./isi),edges); 
  if ~isempty(n), muaisirates(:,cnt)=n; end

  [n,bin]=histc(round(1./allisis),edges); 
  if ~isempty(n), cellisirates(:,cnt)=n; end
  
  % -----------------------------------------------------------
  if QC_flag || plot_series_flag || (nval1==1 && nval2==1)      
    popname=spec.nodes(1).label;
    Vind=find(strcmp([popname '_V'],{data(Eind).sensor_info.label}));
    Vdat=squeeze(data(Eind).epochs.data(Vind,:,selE));
    lfp=LFPs(:,Eind,stim,freq);
    
    if 0
      figure('position',[200 130 1050 840]);
      subplot(2,2,2); plot(1./nextISI(:,1),nextISI(:,2),'o'); xlabel('1/mean cellnetISI [Hz]'); ylabel('stdev ISI');
      subplot(2,2,3); plot(1:nE,1./nextISI(:,1),'o'); ylabel('1/mean cellnetISI [Hz]'); xlabel(Ename);
      subplot(2,2,1); plot(1:nE,nextISI(:,2),'o'); ylabel('stdev ISI'); xlabel(Ename);
      if ~isempty(allisis)
        is=1./allisis; tmpbins=min(is):max(is);
        [h1,b1]=hist(is,tmpbins);
        subplot(2,2,4); bar(b1,h1/length(allisis)); xlabel('1/(cellnetISIs) [Hz]'); ylabel('density');
        try xlim([min(edges) max(edges)]); end
      end      
      figure; 
      subplot(3,1,1); plot(t,tmp,'-',t(pkind),pkval,'*'); ylabel('MUA'); xlabel('time [s]');
      isi=diff(t(pkind)); ylim([-100 50]);
      subplot(3,1,2); 
      [n,bin]=histc(round(1./isi),edges); 
      bar(edges,n/length(isi),'histc'); xlabel('1/muaISI [Hz]'); ylabel('density');
      try xlim([min(edges) max(edges)]); end; title('rhythms based on MUA peaks');
      muaisirates(:,cnt)=n;
      subplot(3,1,3);
      [n,bin]=histc(round(1./allisis),edges); 
      bar(edges,n/length(allisis),'histc'); xlabel('1/cellnetISIs [Hz]'); ylabel('density');    
      try xlim([min(edges) max(edges)]); end; title('rhythms based on spikes to MUA peaks');
      cellisirates(:,cnt)=n;
    end
    
    [sortedISIs,ISIind]=sort(1./nextISI(:,1));
    sorti=symi; % ISIind

    if 0 % ~plot_raster_only && rulesims_flag
      figure('position',[1120 300 850 660]);%[100 150 1700 800])
      subplot(3,2,4);
      plot(t,VI+50,'r',t,VE2-50,'g'); 
      hold on; plot(t,VE,'b','linewidth',3); 
      %plot(t,VI+50,'r',t,VE2-50,'g',t,VE,'b','linewidth',2); 
      legend(Iname,'E(out)','E(in)','Location','SouthWest'); ylabel('MUA'); xlabel('time (s)');
      if isempty(xlims), xlim([0 min(2,t(end))]); else xlim(xlims/1000); end

      subplot(3,2,1); 
      I=InputMatrixE{1}; imagesc(I(sorti,:)); title('Iext->E(in)'); colorbar;
      ylabel([Ename ' index (coh-sorted)']); axis xy; xlabel('time (s)');

      subplot(3,2,3);
      imagesc(t,1:size(Vdat,2),Vdat(:,sorti)'); 
      colormap(1-gray); axis xy; colorbar
      tmp=VE-min(VE); tmp=.2*size(Vdat,2)*tmp/prctile(tmp,99);
      hold on; plot(t,tmp+.5,'b-','linewidth',.5);
      tmp=VE2-min(VE2); tmp=.2*size(Vdat,2)*tmp/prctile(tmp,99.9);
      hold on; plot(t,tmp+.5,'g-','linewidth',1);    
      if isempty(xlims), xlim([0 min(2,t(end))]); else xlim(xlims/1000); end
      ylabel([Ename ' index (coh-sorted)']); title('E(in) Voltage'); xlabel('time (s)');

      subplot(3,2,2);
      tmp=VE2-min(VE2); tmp=10*(tmp/max(tmp));
  %     plot(t,tmp,'g'); xlabel('time (s)'); hold on
      ypos=1;
      if nIs>0
        spks=spiketimes{IslowID(1)};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','k');
          end
        end
      end
      ypos=ypos+i;
      if nI>0
        spks=spiketimes{IfastID(1)};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','r');
          end
        end
      end
      ypos=ypos+i;
      spks=spiketimes{Eind};
      for i=1:length(spks)
        for j=1:length(spks{i})
          line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','b');
        end
      end
      spks=spiketimes{EsomaID(end)};
      ypos=ypos+i;
      for i=1:length(spks)
        for j=1:length(spks{i})
          line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','g');
        end
      end    
      ylim([0 i+1+ypos]); legend('E(out)','E(in)','Location','SouthWest');

      subplot(3,2,6);
      plot(tmins,rI,'r-',tmins,rE2,'g-'); hold on
      plot(tmins,rE,'b-','linewidth',3); 
      %plot(tmins,rI,'r-',tmins,rE2,'g-',tmins,rE,'b-','linewidth',2); 
      legend(Iname,'E(out)','E(in)','Location','SouthWest');%legend(Ename,Iname,Ename2);
      xlabel('time (s)'); ylabel('pop spike rate (Hz)'); 
      xlim([0 max(tmins)]);
      rein=AvgSpikeRates(Eind,stim,freq);
      riin=AvgSpikeRates(Iind,stim,freq);
      rout=AvgSpikeRates(EsomaID(end),stim,freq);
      hline(rein,'b'); hline(riin,'r'); hline(rout,'g');
      text(mean(TOILIM)/1000,rout,sprintf('%gHz',round(rout)));
      text(mean(TOILIM)/1000,rein,sprintf('%gHz',round(rein)));
      text(mean(TOILIM)/1000,riin,sprintf('%gHz',round(riin)));
      vline(TOILIM(1)/1000,'k'); vline(TOILIM(2)/1000,'k');

      subplot(3,2,5); 
      f=muaf;
      Pxx=MUApows(:,Eind,stim,freq);
      OscFreq=muaOscFreqs(1,stim,freq);
      AreaPower=muaAreaPowers(1,stim,freq);
      %plot(f,log10(Pxx),'.-','linewidth',3); xlim([0 50]);% xlim([0 Fmax]); %axis tight    
      plot(f,Pxx,'.-','linewidth',4); xlim([0 50]);% xlim([0 Fmax]); %axis tight
      hold on;
      f=suaf;
      Pxx=avgSUApows(:,1,stim,freq);
      %plot(f,log10(Pxx),'b-'); legend([Ename '-MUA'],[Ename '<SUA>']);    
      plot(f,smooth(Pxx,5),'b-'); 
      legend('E(in)-MUA','E(in)<SUA>');%legend([Ename '-MUA'],[Ename '<SUA>']);    
      vline(OscFreq,'k'); vline(OscFreq-Fwin/2,'r'); vline(OscFreq+Fwin/2,'r');
      %title(sprintf('(f0=%3.3gHz, area=%3.3g)|stim=%g,EsomaType(1)=%g,InputType=%g',OscFreq,AreaPower,Estim,EsomaType(1),InputType));
      title(sprintf('(f0=%3.3gHz, area=%3.3g)',OscFreq,AreaPower));
      ylabel('Power'); xlabel('freq (Hz)');
      ym1=prctile(MUApows(muaf<100,Eind,stim,freq),98);
      ym2=prctile(avgSUApows(suaf<100,Eind,stim,freq),98);
      ylim([0 max(ym1,ym2)])
    end
    
    if plot_spectrogram_flag
      figure('position',[750 290 770 670]);%[1030 550 500 420]);%[600 400 550 600]);
      subplot(3,1,1);
      plot(t,lfp); ylabel('<spikes*psp>'); legend('E LFP');
      title('LFP estimate');%title(sprintf('%s=%g, %s=%g',var1,val1(stim),var2,val2(freq)))
      subplot(3,1,2);
      Fs = fix(1/(t(2)-t(1)));
      NFFT=2^(nextpow2(length(t)-1)-2);
      WINDOW=2^(nextpow2(NFFT-1)-3);
      [yo,fo,to] = spectrogram(lfp,WINDOW,[],NFFT,Fs);
      fsel=find(fo>1 & fo<100);
      yy = (abs(yo(fsel,:))+eps).^2;
      imagesc(to,fo(fsel),yy); axis xy; axis tight; title('spectrogram');      
      ylabel('freq (Hz)'); ylim([0 80]); %colorbar
      caxis([0 prctile(yy(:),99)]);
      subplot(3,1,3)
      yz = (yy-repmat(mean(yy,2),[1 size(yy,2)]))./(repmat(std(yy,0,2),[1 size(yy,2)]));
      imagesc(to,fo(fsel),yz); axis xy; axis tight; title('spectrogram (z-score)');
      xlabel('time (s)'); ylabel('freq (Hz)'); ylim([0 80]); %colorbar;%caxis([-1 1]); %colorbar
      caxis(prctile(yz(:),[2 98]));
    end
    
    if QC_flag && (nval1+nval2)>2, pause, end
%     if Estim==0 %&& gIE==0
%       fprintf('tauI=%g,tauI=%g,gIE=%g,Estim=%g,Enoise=%g: std(V%s)=%-3.2gmV, <r%s>=%3.2gHz, <r%s>=%3.2gHz\n',tauI,tauI,gIE,Estim,Enoise,Ename,std(VE(1000:4000)),Ename,mean(rE),Iname,mean(rI));
%     else
      %fprintf('tauI=%g,tauI=%g,Estim=%g: fnet=%3.2gHz, <r%s>=%3.2gHz, <r%s>=%3.2gHz\n',tauI,tauI,Estim,OscFreq,Ename,mean(rE),Iname,mean(rI));
      fprintf('fnet=%3.2gHz, <r%s>=%3.2gHz, <r%s>=%3.2gHz\n',OscFreq,Ename,mean(rE),Iname,mean(rI));
%     end        
    
  end
  
  if rulesims_flag
    nroi=length(Layer1ResponseROIs);
    nroi2=length(Layer2ResponseROIs);
    
    % define all ROIs
    allrois={Layer1ResponseROIs{:} Layer2ResponseROIs{:}};
    allroid=num2cell(EsomaID(1)*ones(1,nroi));
    allroid=cat(2,allroid,num2cell(EsomaID(end)*ones(1,nroi)));
    if nI2>0
      allrois{end+1}=1:nI2;
      allroid{end+1}=IfastID(end);
    end
    if nIs2>0
      allrois{end+1}=1:nIs2;
      allroid{end+1}=IslowID(end);
    end
    % define all TOIs
    clear alltois
    alltois{1} = [max(StimOn,RuleOn)/1000 TrialOff/1000]; % CONTEXT + STIMULUS (Trial 1)
    StimOn=S1time(1); RuleOn=Context1Time(1); TrialOff=Context1Time(2); % for task analysis
    if tspan(2)>=Context2Time(2)
      alltois{2} = [max(S2time(1),Context2Time(1))/1000 Context2Time(2)/1000];
    end
    if cnt==1 % preallocate result matrices
      % store: (R,p), f0
      ruleRdp=zeros(nval1,nval2);
      rulepdp=nan(nval1,nval2);
      ruleRdr=zeros(nval1,nval2);
      rulepdr=nan(nval1,nval2);
      rulef0=zeros(nval1,nval2);
      ruleRout=zeros(nroi2,nval1,nval2);
      ruleRin=zeros(nroi,nval1,nval2);
      ruleCin=zeros(nroi,nval1,nval2);
    end
    % define center frequencies for freq-avg power
    pkthresh=90;
    smoothfraction=.1; % set to eps instead of 0
    f0predef=[];%40;
    if 1
      Pxx = MUApows(:,Eind,stim,freq);    
      f = muaf;
    else
      Pxx = avgSUApows(:,Eind,stim,freq);
      f = suaf;
    end
    sel = find(FreqRange(1)<=f & f<=FreqRange(end));
    tmpPxx=smooth(Pxx,PxxSmooth);
    ht=prctile(tmpPxx(sel),pkthresh);
    [PeakPower,PPind]=findpeaks(tmpPxx(sel),'MinPeakHeight',ht,'NPeaks',3);
    if ~isempty(PPind)  
      f0s = ceil(f(sel(PPind)));
    else
      f0s = [];
    end
    f0s = [f0s' f0predef];
    f0legend=cellfun(@(x)['f=' num2str(x) 'Hz'],num2cell(f0s),'uni',0);
    nf0=length(f0s);
    
    % spectrogram parameters
    %f0=muaOscFreqs(1,stim,freq); % center freq for spectral integration
    %fpad=5; % : avg over f0+/-fpad
    Fs = fix(1/(t(2)-t(1)));
    NFFT=2^(nextpow2(length(t)-1)-2);
    WINDOW=2^(nextpow2(NFFT-1)-3);
    E1=EsomaID(1); E2=EsomaID(end);
    % 
    taud2=10;%5;%10
    taur2=.9*taud;
    psp2 = (exp(-tt/(taud2/1000))-exp(-tt/(taur2/1000)));
    % psp = gaussian
    psp2 = psp2/max(psp2);
    psp2 = [zeros(1,length(psp2)) psp2];    
    % compute ROI averages (LFP spectrograms, ..)
    Rroi=zeros(length(tmins),nroi);
    covroi=zeros(length(t),nroi);
    % TFRs=
    for roi=1:nroi
      I=Layer1ResponseROIs{roi};
      % compute ROI LFP
      lfp=zeros(length(tt),1);
      thiscov=zeros(length(t),1);
      nelm=length(I);
      denom=0;
      for i=1:nelm
        spks=spiketimes{E1}{I(i)};
        X=histc(spks,tt); if size(X,1)==1, X=X'; end
        X=conv(X,psp,'same');
        lfp=lfp+X/nelm;
        X=conv(X,psp2,'same');
        for j=i+1:nelm
          spks=spiketimes{E1}{I(j)};
          Y=histc(spks,tt); if size(Y,1)==1, Y=Y'; end
          Y=conv(Y,psp2,'same');
          tmpcov=conv(X,Y,'same');
          thiscov=thiscov+tmpcov;
          denom=denom+1;
          %[acor,lag] = xcorr(X,Y);
        end
      end
      covroi(:,roi)=smooth(thiscov,100)/denom;
      % compute ROI spectrogram
      [yo,fo,to] = spectrogram(lfp,WINDOW,[],NFFT,Fs);
      if roi==1
        Ptfr = zeros(length(to),nroi,nf0);
      end
      if plot_spectrogram_flag
        fsel=find(fo>(0) & fo<(100)); % freq range to avg
        yy = (abs(yo(fsel,:))+eps).^2;
        if roi==1
          figure('position',[1045 220 550 680]); 
          clims=prctile(yy(:),[0 97]);%[min(yy(:)) max(yy(:))];
        end
        subplot(nroi+1,1,nroi-roi+1);
        imagesc(to,fo(fsel),yy); axis xy; axis tight; 
        title(sprintf('E(in) roi %g',roi)); caxis(clims); colorbar
        ylabel('freq (Hz)'); xlabel('time (sec)');
        for ff=1:nf0
          hline(f0s(ff),'k');
        end      
        if roi==nroi
          subplot(nroi+1,1,nroi+1);
          Pxx2 = MUApows(:,Eind,stim,freq); 
          plot(muaf,Pxx2,'b-','linewidth',2); xlim([0 80]);
          hold on
          Pxx2 = avgSUApows(:,Eind,stim,freq);
          plot(suaf,smooth(Pxx2,5),'b--');
          Pxx2 = MUApows(:,EsomaID(end),stim,freq);
          plot(muaf,Pxx2,'g-','linewidth',2); xlim([0 80]);          
          Pxx2 = avgSUApows(:,EsomaID(end),stim,freq);
          plot(suaf,smooth(Pxx2,5),'g--');
          legend('E(in)-MUA','E(in)<SUA>','E(out)-MUA','E(out)-SUA');
          for ff=1:nf0
            vline(f0s(ff),'k');
          end      
          ylabel('Power'); xlabel('freq (Hz)'); 
          ym1=prctile(MUApows(muaf<100,Eind,stim,freq),98);
          ym2=prctile(avgSUApows(suaf<100,Eind,stim,freq),98);
          ylim([0 max(ym1,ym2)])          
          %ylim([0 10])
        end
      end      
      for ff=1:nf0
        f0=f0s(ff);
        fpad=.2*f0;
        fsel=find(fo>(f0-fpad) & fo<(f0+fpad)); % freq range to avg
        P = (abs(yo(fsel,:))+eps).^2; % compute power P(f,t)
        Pavg = mean(P,1); % freq-average power
        Ptfr(:,roi,ff)=smooth(Pavg,ceil(smoothfraction*numel(to)));
      end
      % compute ROI spike rates
      Rroi(:,roi)=smooth(mean(rates{E1}(I,:),1),ceil(smoothfraction*size(Rroi,1)));
    end
    
    % compute output ROI rates
    I=Layer2ResponseROIs{1};
    rL=mean(rates{E2}(I,:),1);
    if numel(Layer2ResponseROIs)>1
      I=Layer2ResponseROIs{2};
      rR=mean(rates{E2}(I,:),1);
    else
      rR=zeros(size(rL));
    end
    rL=smooth(rL,ceil(smoothfraction*numel(rL)));
    rR=smooth(rR,ceil(smoothfraction*numel(rL)));
    % compute differentials
    if nroi>1
      dP = 100*(squeeze(Ptfr(:,1,:))-squeeze(Ptfr(:,2,:)))/max(Ptfr(:)); % plot(to,dP)
      drAD = 100*(Rroi(:,1)-Rroi(:,2))/max(Rroi(:)); % plot(tmins,drAD)
      dcov = 100*(covroi(:,1)-covroi(:,2))/max(covroi(:));
    else
      dP = 100*squeeze(Ptfr(:,1,:))/max(Ptfr(:)); % plot(to,dP)
      drAD = 100*Rroi(:,1)/max(Rroi(:)); % plot(tmins,drAD)
      dcov = 100*covroi(:,1)/max(covroi(:));
    end
    drLR = 100*(rL-rR)/max(max(rL),max(rR)); % plot(tmins,drLR)
    [jnk,dPord]=sort(range(dP,1),2,'descend');
    dp = dP(:,dPord(1));
    S1t=Context1Time(1)/1000; S2t=Context2Time(1)/1000;
    
    dp = smooth(dp,ceil(smoothfraction*numel(dp)));
    drAD = smooth(drAD,ceil(smoothfraction*numel(drAD)));
    drLR = smooth(drLR,ceil(smoothfraction*numel(drLR)));
    
    a=S1time(1)/1000; b=S1time(2)/1000; c=S2time(1)/1000; d=S2time(2)/1000;
    tselFR = find((tmins>=a&tmins<=b)|(tmins>=c&tmins<=d));
    tselPA = find((to>=a&to<=b)|(to>=c&to<=d));
    
    dplot=zeros(size(dp));
    for i=1:length(to)
      ind=nearest(tmins,to(i));
      dplot(i)=drLR(ind);
    end
    if 0
      selP=tselPA;
      selR=tselFR;
      Tstr='stim-epochs';
    else
      selP=1:length(dp);
      selR=1:length(rL);
      Tstr=sprintf('%g-%gs',xlims); %'full-trials';
    end
    [Rdp,pdp]=corr(dp(selP),dplot(selP));
    [Rdr,pdr]=corr(drAD(selR),drLR(selR));

    ruleRdp(stim,freq)=Rdp;
    rulepdp(stim,freq)=pdp;
    ruleRdr(stim,freq)=Rdr;
    rulepdr(stim,freq)=pdr;
    rulef0(stim,freq)=f0s(dPord(1));
    
    if ~plot_raster_only || summary_plots
      % Plot differential output:
      figure('position',[620 80 840 830],'visible',visible_status); % [500 110 840 825]
      xlims=[0 max(tmins)]; 
      tmp=max([max(abs(dP(:))) max(abs(dP(:)))]); ylimsP=[-tmp tmp];
      if ylimsP(1)==ylimsP(2), ylimsP(1)=ylimsP(1)-.5; ylimsP(2)=ylimsP(2)+.5; end
      subplot(3,2,1); % P(A)-P(D)
      plot(to,dP(:,dPord),'linewidth',2); ylim(ylimsP); xlim(xlims);
      xlabel('time (sec)'); ylabel('P(A)-P(D) (%max)'); hline(0,'k');
      vline(S1t,'k'); vline(S2t,'k');
      hold on; plot(to,dp,'b-','linewidth',4);
      legend({f0legend{dPord},'smoothed'},'Location','SouthWest');
      ht=.9*ylimsP(2); line([a b],[ht ht]); line([c d],[ht ht]);
      subplot(3,2,3); % r(A)-r(D)
      plot(tmins,drAD,'linewidth',2); xlabel('time (sec)'); ylabel('r(A)-r(D) (%max)'); hline(0,'k');
      tmp=max([max(abs(drAD(:))) max(abs(drLR(:)))]); 
      if isnan(tmp), tmp=0; end; ylimsR=[-tmp tmp];
      if ylimsR(1)==ylimsR(2), ylimsR(1)=ylimsR(1)-.5; ylimsR(2)=ylimsR(2)+.5; end      
      xlim(xlims); ylim(ylimsR); vline(S1t,'k'); vline(S2t,'k'); 
      ht=.9*ylimsR(2); line([a b],[ht ht]); line([c d],[ht ht]);
      subplot(3,2,5); % r(L)-r(R)
      plot(tmins,drLR,'linewidth',2); xlabel('time (sec)'); ylabel('r(L)-r(R) (%max)'); hline(0,'k');
      xlim(xlims); ylim(ylimsR); vline(S1t,'k'); vline(S2t,'k');  
      ht=.9*ylimsR(2); line([a b],[ht ht]); line([c d],[ht ht]);
      subplot(3,2,2); % r(L)-r(R) vs. P(A)-P(D). lsline. text(R^2)
      scatter(dp(selP),dplot(selP)); xlabel(sprintf('P(A)-P(D) | t=(%s), f=%gHz',Tstr,f0s(dPord(1)))); ylabel('r(L)-r(R)');
      xlim(ylimsP); ylim(ylimsR);
      lsline; title(sprintf('R=%3.3g (p=%3.3g)',Rdp,pdp)); hline(0,'k');  vline(0,'k');      
      subplot(3,2,4); % r(L)-r(R) vs. r(A)-r(D). lsline. text(R^2)
      scatter(drAD(selR),drLR(selR)); xlim(ylimsR); ylim(ylimsR); xlabel(sprintf('r(A)-r(D) | t=(%s)',Tstr)); ylabel('r(L)-r(R)');
      lsline; title(sprintf('R=%3.3g (p=%3.3g)',Rdr,pdr)); hline(0,'k');  vline(0,'k');
      xlim(ylimsR); ylim(ylimsR);
      subplot(3,2,6); % spike covariance
      %plot(t,covroi); legend(cellfun(@(x)['roi' num2str(x)],num2cell(1:nroi),'uni',0));
      denv=abs(hilbert(dcov));
      plot(t,dcov,'b-',t,denv,'r-'); %plot(t,dcov); 
      legend('roi1-roi2','envelope'); xlim(xlims);
      vline(S1t,'k'); vline(S2t,'k'); hline(0,'k');
      xlabel('time (sec)'); ylabel('spike covariance (A-D)');
    end
    
    if ~plot_raster_only && plotinputs
      % Plot inputs (context cues and stimuli) and outputs (L,R):
      C1=R1; C2=R2; S=S1+S2;
      figure('position',[70 80 560 890]);
      subplot(4,1,1); % C1
      ts=(0:size(C1,2)-1)*dt/1000;
      plot(ts,C1); xlabel('time'); ylabel('Context 1');
      xlim(xlims); vline(S1t,'k'); vline(S2t,'k'); 
      subplot(4,1,2); % C2
      plot(ts,C2); xlabel('time'); ylabel('Context 2');
      xlim(xlims);vline(S1t,'k'); vline(S2t,'k'); 
      subplot(4,1,3); % S
      plot(ts,S); xlabel('time'); ylabel('Stimulus');
      xlim(xlims);vline(S1t,'k'); vline(S2t,'k'); 
      subplot(4,1,4); % r(L)
      plot(tmins,drLR); xlabel('time'); ylabel('r(L)-r(R)'); hline(0,'k');
      xlim(xlims); vline(S1t,'k'); vline(S2t,'k'); 
  %     subplot(5,1,4); % r(L)
  %     plot(tmins,rL); xlabel('time'); ylabel('FR(left)');
  %     subplot(5,1,5); % r(R)    
  %     plot(tmins,rR); xlabel('time'); ylabel('FR(right)');
    end
    
    % Stats
    fprintf('\nROI#: L1[#spiking](<#spikes>,spkcoh)-->L2[#spiking](<#spikes>)\n');
    E1=EsomaID(1);
    E2=EsomaID(2);
    % -------------------------------------
    % TOI 1: constrain to spikes BEFORE (RULE+STIM)
    % -------------------------------------
%     spk1=spiketimes{E1}; % layer 1 E-cells
%     spk2=spiketimes{E2}; % layer 2 E-cells
%     % restrict to times of interest
%     mintime=max(StimOn,RuleOn);
%     for i=1:length(spk1), spk1{i}=spk1{i}(spk1{i}<=mintime/1000); end
%     for i=1:length(spk2), spk2{i}=spk2{i}(spk2{i}<=mintime/1000); end
%     % spike counts
%     nspk1=cellfun(@length,spk1); % layer 1: spike counts (all E-cells)
%     nspk2=cellfun(@length,spk2); % layer 2: spike counts (all E-cells)
%     for i=1:length(Layer1ResponseROIs) % loop over response ROIs
%       a=Layer1ResponseROIs{i}; % this input-layer response ROI
%       b=Layer2ResponseROIs{i}; % this output-layer response ROI
%       aa=a(nspk1(a)~=0); % spiking cells in this L1 ROI
%       bb=b(nspk2(b)~=0); % spiking cells in this L2 ROI
%       % calc mean spike count for input L1 ROI
%       ra=mean(nspk1(aa)); na=length(aa); % # spiking cells
%       % calc spike coherence for input L1 ROI
%       coh=nan(na,na);
%       for j=1:na
%         x1=histc(spk1{aa(j)},tt); if size(x1,1)==1, x1=x1'; end
%         x1=conv(x1,psp,'same');
%         for k=1:na
%           x2=histc(spk1{aa(k)},tt); if size(x2,1)==1, x2=x2'; end
%           x2=conv(x2,psp,'same');
%           coh(j,k)=sum(x1.*x2)./sqrt(sum(x1.^2).*sum(x2.^2));
%         end
%       end
%       coha=nanmean(coh(:));
%       % calc mean spike count for output L2 ROI
%       rb=mean(nspk2(bb)); nb=length(bb); % # spiking cells
%       % display results: 
%       if na==0
%         fprintf('ROI%g(null): L1[%g/%g](%gspk,%3.3gcoh)-->L2[%g/%g](%gspk)',i,na,length(a),ra,coha,nb,length(b),rb);
%       else
%         fprintf('ROI%g(%g-%g): L1[%g/%g](%gspk,%3.3gcoh)-->L2[%g/%g](%gspk)',i,aa(1),aa(end),na,length(a),ra,coha,nb,length(b),rb);
%       end
%       if ismember(i,Context1), fprintf(' [rule]\n'); else fprintf('\n'); end
%     end
%     fprintf('  (over %g-%gms)\n\n',tspan(1),mintime);    
    
    % -------------------------------------
    % NOW AGAIN: constrain to spikes DURING (RULE+STIM)
    % -------------------------------------
    % calc coherence and mean rates
    spk1=spiketimes{E1}; % layer 1 E-cells
    spk2=spiketimes{E2}; % layer 2 E-cells
    % restrict to times of interest
    mintime=max(StimOn,RuleOn);
    maxtime=TrialOff/1000;
    tdur=maxtime-mintime/1000;
    for i=1:length(spk1), spk1{i}=spk1{i}(spk1{i}>=mintime/1000&spk1{i}<=maxtime); end
    for i=1:length(spk2), spk2{i}=spk2{i}(spk2{i}>=mintime/1000&spk2{i}<=maxtime); end
    % spike counts
    nspk1=cellfun(@length,spk1); % layer 1: spike counts (all E-cells)
    nspk2=cellfun(@length,spk2); % layer 2: spike counts (all E-cells)
    rout=zeros(1,length(Layer2ResponseROIs));
    rin=zeros(1,length(Layer1ResponseROIs));
    Cin=zeros(1,length(Layer1ResponseROIs));
    for i=1:length(Layer1ResponseROIs) % loop over response ROIs
      a=Layer1ResponseROIs{i}; % this input-layer response ROI
      b=Layer2ResponseROIs{i}; % this output-layer response ROI
      aa=a(nspk1(a)~=0); % spiking cells in this L1 ROI
      bb=b(nspk2(b)~=0); % spiking cells in this L2 ROI
      % calc mean spike count for input L1 ROI
      ra=mean(nspk1(aa)); na=length(aa); % # spiking cells
      % calc spike coherence for input L1 ROI
      coh=nan(na,na);
      for j=1:na
        x1=histc(spk1{aa(j)},tt); if size(x1,1)==1, x1=x1'; end
        x1=conv(x1,psp,'same');
        for k=1:na
          x2=histc(spk1{aa(k)},tt); if size(x2,1)==1, x2=x2'; end
          x2=conv(x2,psp,'same');
          coh(j,k)=sum(x1.*x2)./sqrt(sum(x1.^2).*sum(x2.^2));
        end
      end
      coha=nanmean(coh(:));
      % calc mean spike count for output L2 ROI
      rb=mean(nspk2(bb)); nb=length(bb); % # spiking cells
      % display results: 
      if na==0
        fprintf('ROI%g(null): L1[%g/%g](%gspk,%3.3gcoh)-->L2[%g/%g](%gspk)',i,na,length(a),ra,coha,nb,length(b),rb);
      else
        fprintf('ROI%g(%g-%g): L1[%g/%g](%gspk,%3.3gcoh)-->L2[%g/%g](%gspk)',i,aa(1),aa(end),na,length(a),ra,coha,nb,length(b),rb);
      end
      if ismember(i,Context1), fprintf(' [rule]\n'); else fprintf('\n'); end
      rin(i)=ra/tdur;
      rout(i)=rb/tdur;
      Cin(i)=coha;
    end
    fprintf('  (over %g-%gms)\n\n',mintime,tspan(2));
    
    ruleRout(:,stim,freq)=rout;
    ruleRin(:,stim,freq)=rin;
    ruleCin(:,stim,freq)=Cin;
    
    FreqRange=[max(Fmin,2/t(end)) Fmax]; % frequencies to consider
    % -------------------------------------------------------
%     % ROIs
%     Layer1ResponseROIs (EsomaID(1)), Layer2ResponseROIs (EsomaID(end))
%     if nI2>0, {1:nI2} IfastID(end); end
%     if nIs2>0, {1:nIs2} IfastID(end); end
%     % TOIs
%     stim: [max(StimOn,RuleOn) TrialOff/1000]    
    nrois=length(allrois);
    ntois=length(alltois);
    if cnt==1
      ROIrate=zeros(nrois,ntois,nval1,nval2);
      ROImuaf0=zeros(nrois,ntois,nval1,nval2);
      ROImuaPA=zeros(nrois,ntois,nval1,nval2);
      ROImuaPmax=zeros(nrois,ntois,nval1,nval2);
      ROIsuaf0=zeros(nrois,ntois,nval1,nval2);
      ROIsuaPA=zeros(nrois,ntois,nval1,nval2);
      ROIsuaPmax=zeros(nrois,ntois,nval1,nval2);
      ROIlfpf0=zeros(nrois,ntois,nval1,nval2);
      ROIlfpPA=zeros(nrois,ntois,nval1,nval2);
      ROIlfpPmax=zeros(nrois,ntois,nval1,nval2);
    end
    for i=1:ntois
      tmin=alltois{i}(1);
      tmax=alltois{i}(2);
      t1=nearest(t,tmin); % sec, TOI start time
      t2=nearest(t,tmax); % sec, TOI stop time
      tdur=tmax-tmin;
      for j=1:nrois
        id=allroid{j};  % population index
        cellsel=allrois{j}; % cell indices
        % ROI data
        popname=spec.nodes(id).label;
        Vind=find(strcmp([popname '_V'],{data(id).sensor_info.label}));
        ROIdata=squeeze(data(id).epochs.data(Vind,t1:t2,cellsel));    
        if length(cellsel)==1, ROIdata=ROIdata'; end
        MUA=double(nanmean(ROIdata,2))';
        % Spike rate
        spk=spiketimes{id}(cellsel);
        for k=1:length(spk), spk{k}=spk{k}(spk{k}>=tmin&spk{k}<=tmax); end
        nspk=cellfun(@length,spk);
        ROIrate(j,i,stim,freq)=mean(nspk(nspk~=0))/tdur;
        % SUA spectrum
        NFFT=2^(nextpow2(size(ROIdata,1)-1)-1);%2); % <-- use higher resolution to capture STO freq variation
        WINDOW=2^(nextpow2(NFFT-1)-3);
        Pxx=0;
        for k=1:size(ROIdata,2)
          X=ROIdata(:,k);
          [tmpPxx,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);                  
          if all(isnan(tmpPxx(:)))
            tmpPxx=zeros(size(tmpPxx));
          end
          Pxx=Pxx+tmpPxx/size(ROIdata,2);
        end
        fsel = find(FreqRange(1)<=f & f<=FreqRange(end));
        tmpPxx=smooth(Pxx,PxxSmooth);
        ht=prctile((tmpPxx(fsel)),powthreshprc);
        [PeakPower,PPind]=findpeaks((tmpPxx(fsel)),'MinPeakHeight',ht,'NPeaks',3);
        if ~isempty(PPind)
          PPind=PPind(PeakPower==max(PeakPower));
          OscFreq = f(fsel(PPind));
          flo=OscFreq-Fwin/2;
          fhi=OscFreq+Fwin/2;
          sel2=find(flo<=f & f<=fhi);
          ROIsuaf0(j,i,stim,freq)=OscFreq;
          ROIsuaPA(j,i,stim,freq)=sum(tmpPxx(sel2))*(f(2)-f(1));%sum(Pxx(sel2))*(f(2)-f(1));
          ROIsuaPmax(j,i,stim,freq)=max(tmpPxx(fsel));%max(Pxx(fsel));
        end
        suaPxx=tmpPxx;%Pxx;
        suafaxis=f;
        
        % MUA spectrum
        X=MUA;
        NFFT=2^(nextpow2(length(X)-1)-2);
        WINDOW=2^(nextpow2(NFFT-1)-3);
        X=smooth(X,ceil(.1/dt));
        [Pxx,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);
        fsel = find(FreqRange(1)<=f & f<=FreqRange(end));
        tmpPxx=smooth(Pxx,PxxSmooth);
        ht=prctile((tmpPxx(fsel)),powthreshprc);
        [PeakPower,PPind]=findpeaks((tmpPxx(fsel)),'MinPeakHeight',ht,'NPeaks',3);
        if ~isempty(PPind)
          PPind=PPind(PeakPower==max(PeakPower));
          OscFreq = f(fsel(PPind));
          flo=OscFreq-Fwin/2;
          fhi=OscFreq+Fwin/2;
          sel2=find(flo<=f & f<=fhi);
          ROImuaf0(j,i,stim,freq)=OscFreq;
          ROImuaPA(j,i,stim,freq)=sum(tmpPxx(sel2))*(f(2)-f(1));%sum(Pxx(sel2))*(f(2)-f(1));
          ROImuaPmax(j,i,stim,freq)=max(tmpPxx(fsel));%max(Pxx(fsel));
        end
        muaPxx=tmpPxx;%Pxx;
        muafaxis=f;

        % LFP and LFP spectrum
        taud2=10; taur2=.9*taud;
        psp2 = (exp(-t(t1:t2)/(taud2/1000))-exp(-t(t1:t2)/(taur2/1000)));
        psp2 = psp2/max(psp2);
        psp2 = [zeros(1,length(psp2)) psp2];    
        lfp=0;
        for k=1:size(ROIdata,2)
          X=histc(spk{k},t); if size(X,1)==1, X=X'; end
          X=conv(X,psp2,'same');
          lfp=lfp+X/size(ROIdata,2);
        end
        %[yo,fo,to] = spectrogram(lfp,WINDOW,[],NFFT,Fs);
        X=smooth(lfp,ceil(.1/dt))*range(MUA(:));
        NFFT=2^(nextpow2(length(X)-1)-2);
        WINDOW=2^(nextpow2(NFFT-1)-3);
        [Pxx,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);
        fsel = find(FreqRange(1)<=f & f<=FreqRange(end));
        tmpPxx=smooth(Pxx,PxxSmooth);
        ht=prctile((tmpPxx(fsel)),powthreshprc);
        [PeakPower,PPind]=findpeaks((tmpPxx(fsel)),'MinPeakHeight',ht,'NPeaks',3);
        if ~isempty(PPind)
          PPind=PPind(PeakPower==max(PeakPower));
          OscFreq = f(fsel(PPind));
          flo=OscFreq-Fwin/2;
          fhi=OscFreq+Fwin/2;
          sel2=find(flo<=f & f<=fhi);
          ROIlfpf0(j,i,stim,freq)=OscFreq;
          ROIlfpPA(j,i,stim,freq)=sum(tmpPxx(sel2))*(f(2)-f(1));%sum(Pxx(sel2))*(f(2)-f(1));
          ROIlfpPmax(j,i,stim,freq)=max(tmpPxx(fsel));%max(Pxx(fsel));
        end
        lfpPxx=tmpPxx;%Pxx;
        lfpfaxis=f;
        
        if ~plot_raster_only
          % plot (ROI,TOI) spectra: suaPxx, muaPxx, lfpPxx
          flims=[0 100];
          if j==1
            figure('position',[470 70 1060 850]); 
            nr=3; nc=2;
          end
          subplot(nr,nc,j);
          plot(muafaxis,muaPxx,'b-','linewidth',2); xlim(flims);
          hold on
          plot(suafaxis,suaPxx,'g-','linewidth',2);
%           if strncmp('E',popname,1)
%             plot(lfpfaxis,lfpPxx,'r-','linewidth',2);
%             legend('MUA','<SUA>','LFP->');
%           else
            legend('MUA','<SUA>');
%           end
          xlabel('freq (Hz)'); ylabel('power');
          r=ROIrate(j,i,stim,freq);
          f0=ROImuaf0(j,i,stim,freq);
          pa=ROImuaPA(j,i,stim,freq);
          pmax=ROImuaPmax(j,i,stim,freq);
          ylim([0 100]); vline(f0,'k'); 
          text(f0,pmax,sprintf('f=%-3.3gHz (PA=%-3.3g)',f0,pa));
          text(5,pmax,sprintf('FR=%-3.3gHz',r));
          title(sprintf('%s(%g-%g) [%g-%gs]',popname,cellsel(1),cellsel(end),tmin,tmax));
        end
      end
      if ~plot_raster_only
        tmp1=ROIsuaPmax(:,:,stim,freq);
        tmp2=ROImuaPmax(:,:,stim,freq);
        ylims=[0 max(max(tmp1(:)),max(tmp2(:)))];
        hAllAxes = findobj(gcf,'type','axes');
        hLeg = findobj(hAllAxes,'tag','legend');
        hAxes = setdiff(hAllAxes,hLeg);      
        set(hAxes(1:4),'ylim',ylims);
      end
    end
    if ~plot_raster_only
      figure; 
      for j=1:nrois
        id=allroid{j};  % population index
        cellsel=allrois{j}; % cell indices
        popname=spec.nodes(id).label;
        ravg=mean(rates{id}(cellsel,:),1);
        subplot(nr,nc,j);
        plot(tmins,ravg,'b'); axis tight
        xlabel('time (s)'); ylabel('FR (Hz)');
        title(sprintf('%s(%g-%g)',popname,cellsel(1),cellsel(end)));
      end
    end
    
    if 1%plot_raster_only
      figure('position',[220 300 810 590],'visible',visible_status);
      xmin=min(StimOn,RuleOn)/1000;
      ypos=1; yticks=[]; yticklabels={};
      xlim([xmin t(end)]);
      if nIs>0
        spks=spiketimes{IslowID(1)};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','k');
          end
          if i==nIsCellsPerRule
            hline(i+ypos+.5,'color','k','linestyle','--');
          end
        end
        yticks(end+1)=ypos+i/2; yticklabels{end+1}='Islow';
        hline(i+ypos+.5,'color','k','linewidth',3);
        hline(1+.5,'color','k','linewidth',3);
      end
      ypos=ypos+i;
      if nI>0
        spks=spiketimes{IfastID(1)};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','r');
          end
          if i==nICellsPerRule
            hline(i+ypos+.5,'color','r','linestyle','--');
          end
        end
        yticks(end+1)=ypos+i/2; yticklabels{end+1}='Ifast';
        hline(i+ypos+.5,'color','k','linewidth',3);
      end
      ypos=ypos+i;
      spks=spiketimes{Eind};
      for i=1:length(spks)
        for j=1:length(spks{i})
          line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','b');
        end
        if i==nECellsPerBit
          hline(i+ypos+.5,'color','b','linestyle','--');
        end
      end
      hline(i+ypos+.5,'color','k','linewidth',3);
      yticks(end+1)=ypos+i/2; yticklabels{end+1}='E(in)';
      % add metrics
      for j=1:length(rin)
        xloc=xmin+.1;
        yloc=ypos+nECellsPerBit/2+(j-1)*nECellsPerBit;
        yoff=.25*nECellsPerBit;
        text(xloc,yloc+yoff,sprintf('<r>=%gHz',round(rin(j))));
        text(xloc,yloc     ,sprintf('ro/ri=%3.3g',rout(j)/rin(j)));
        text(xloc,yloc-yoff,sprintf('coh=%3.3g',Cin(j)));
      end
      spks=spiketimes{EsomaID(end)};
      ypos=ypos+i;
      for i=1:length(spks)
        for j=1:length(spks{i})
          line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','g');
        end
        if i==nECellsPerResponse
          hline(i+ypos+.5,'color','g','linestyle','--');
        end
      end
      hline(i+ypos+.5,'color','k','linewidth',3);
      yticks(end+1)=ypos+i/2; yticklabels{end+1}='E(out)';
      % add metrics
      for j=1:length(rout)
        xloc=xmin+.1;
        yoff=.5*nECellsPerResponse;
        yloc=ypos+nECellsPerResponse/2+(j-1)*nECellsPerResponse;
        text(xloc,yloc+yoff,sprintf('<r>=%gHz',round(rout(j))));
      end      
      if nI2>0
        ypos=ypos+i;
        spks=spiketimes{IfastID(2)};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','r');
          end
        end
        yticks(end+1)=ypos+i/2; yticklabels{end+1}='Ifast2';
        hline(i+ypos+.5,'color','k','linewidth',3);
      end
      if nIs2>0
        ypos=ypos+i;
        spks=spiketimes{IslowID(2)};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','k');
          end
        end
        yticks(end+1)=ypos+i/2; yticklabels{end+1}='Islow2';
        hline(i+ypos+.5,'color','k','linewidth',3);
      end     
      
      if include_thalamus_flag
        hline(i+ypos+.5,'color','k','linewidth',6);
        ypos=ypos+.25;
        ypos=ypos+i;
        spks=spiketimes{TC_ID};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','b');
          end
          if i==nTCPerResponse
            hline(i+ypos+.5,'color','b','linestyle','--');
          end
        end
        yticks(end+1)=ypos+i/2; yticklabels{end+1}='TC';
        hline(i+ypos+.5,'color','k','linewidth',3);
        ypos=ypos+i;
        spks=spiketimes{RE_ID};
        for i=1:length(spks)
          for j=1:length(spks{i})
            line([spks{i}(j) spks{i}(j)],[i+ypos-.5 i+ypos+.5],'color','r');
          end
        end
        yticks(end+1)=ypos+i/2; yticklabels{end+1}='RE';
        hline(i+ypos+.5,'color','k','linewidth',3);
      end
      ylim([0 i+1+ypos]);
      set(gca,'ytick',yticks,'yticklabel',yticklabels);
      xlabel('time (sec)');
      title(sprintf('rastergram (%s=%g, %s=%g)',var1,val1(stim),var2,val2(freq)))
      vline(maxtime,'k');
    end
    
  end
  
  % plot currents; analyze IPSPs and h-current
  if plotcurrents_flag
    pset.p=model.model.parameters;
    timelimits=tspan;
    av=model.model.auxvars;
    for a=1:size(av,1)
      eval(sprintf('%s=%s;',av{a,1},av{a,2}));
    end
    fn=model.model.functions;
    for a=1:size(fn,1)
      eval(sprintf('%s=%s;',fn{a,1},fn{a,2}));
    end
    
    SDIM1=CatA(1); SDIM2=CatB(2);
    
    nt=length(t);
    E1=find(strcmp('E1',{spec.nodes.label}));
    If1=find(strcmp('Ifast1',{spec.nodes.label}));
    Is1=find(strcmp('Islow1',{spec.nodes.label}));    
    selIf = find(sum(EIfmask{1},2)>0,1,'first');
    selIs = find(sum(EIsmask{1},2)>0,1,'first');
    var=strcmp('E1_V',{data(E1).sensor_info.label});
    E_V = squeeze(data(E1).epochs.data(var,:,1:nE)); % selIf));
    var=strcmp('Ifast1_V',{data(If1).sensor_info.label});
    If_V = squeeze(data(If1).epochs.data(var,:,1:nI)); % selIf));
    var=strcmp('Islow1_V',{data(Is1).sensor_info.label});
    Is_V = squeeze(data(Is1).epochs.data(var,:,1:nIs)); % selIs));
    
    % plot E-cells
    % V(t): (E1(1-2)/If1(1); E1(7-8)/Is1(2))
    % I(t): IeeAMPA, IeeNMDA, IeeAMPA+IeeNMDA. IieGABA
    % g(t): sAMPAee, sNMDAee, sGABAie
    
    % Synaptic conductances of excitatory cells
    var=strcmp('E1_E1_AMPA_s',{data(E1).sensor_info.label});
    E_s_AMPA = squeeze(data(E1).epochs.data(var,:,1:nE));
    var=strcmp('E1_E1_NMDA_sNMDApre',{data(E1).sensor_info.label});
    E_s_NMDA = squeeze(data(E1).epochs.data(var,:,1:nE));
    var=strcmp('Ifast1_E1_GABAa_s',{data(E1).sensor_info.label});
    IfE_s_GABA = squeeze(data(E1).epochs.data(var,:,1:nI));
    var=strcmp('Islow1_E1_GABAa_s',{data(E1).sensor_info.label});
    IsE_s_GABA = squeeze(data(E1).epochs.data(var,:,1:nIs));
    % Synaptic currents to excitatory cells
    E_iAMPA=zeros(nt,nE);
    E_iNMDA=zeros(nt,nE);
    If_E_iGABA=zeros(nt,nE);
    Is_E_iGABA=zeros(nt,nE);
    for a=1:nt
      E_iAMPA(a,:) = E1_E1_AMPA_ISYN(E_V(a,:)',E_s_AMPA(a,:)');
      E_iNMDA(a,:) = E1_E1_NMDA_INMDA(E_V(a,:)',E_s_NMDA(a,:)');
      If_E_iGABA(a,:) = Ifast1_E1_GABAa_ISYN(E_V(a,:)',IfE_s_GABA(a,:)');
      Is_E_iGABA(a,:) = Islow1_E1_GABAa_ISYN(E_V(a,:)',IsE_s_GABA(a,:)');
    end
    % generate plots
    figure('position',[70 50 1800 850]);
    tlims=[1 1.5];%1.5];%[1 2];
    mintime=max(StimOn,RuleOn);
    subplot(4,2,1); plot(t,E_V(:,SDIM1),'b',t,E_V(:,SDIM1+1),'b--',t,If_V(:,selIf),'r'); ylabel('V(t)'); 
    xlim(tlims); ylims=ylim; title('voltage');
    legend(sprintf('E1(cell%g)',SDIM1),sprintf('E1(cell%g)',SDIM1+1),sprintf('If1(cell%g)',selIf),'Location','NorthWest');
    subplot(4,2,2); plot(t,E_V(:,SDIM2),'b',t,E_V(:,SDIM2+1),'b--',t,Is_V(:,selIs),'r'); ylabel('V(t)'); 
    xlim(tlims); ylim(ylims); title('voltage'); 
    legend(sprintf('E1(cell%g)',SDIM2),sprintf('E1(cell%g)',SDIM2+1),sprintf('Is1(cell%g)',selIs),'Location','NorthWest');
    subplot(4,2,3);
    plot(t,E_iAMPA(:,SDIM1)+E_iNMDA(:,SDIM1),'k',t,If_E_iGABA(:,SDIM1),'r'); 
    legend(sprintf('E1(cell%g)<-iAMPA+iNMDA',SDIM1),'E1<-If1.iGABAa','Location','NorthWest'); ylabel('I(t)'); 
    xlim(tlims); ylims=ylim; title(sprintf('E1(cell%g) inputs (E/I)',SDIM1));
    subplot(4,2,4);
    plot(t,E_iAMPA(:,SDIM2)+E_iNMDA(:,SDIM2),'k',t,Is_E_iGABA(:,SDIM2),'r'); 
    legend(sprintf('E1(cell%g)<-iAMPA+iNMDA',SDIM2),'E1<-Is1.iGABAa','Location','NorthWest'); ylabel('I(t)'); 
    xlim(tlims); ylim(ylims); title(sprintf('E1(cell%g) inputs (E/I)',SDIM2));
    subplot(4,2,5); 
    plot(t,E_iAMPA(:,SDIM1),'b',t,E_iNMDA(:,SDIM1),'g');%,t,If_AMPA+If_NMDA,'k--'); 
    legend('E1<-iAMPA','E1<-iNMDA','Location','NorthWest'); ylabel('I(t)'); xlim(tlims);
    title(sprintf('E1(cell%g) inputs (AMPA,NMDA)',SDIM1)); ylims=ylim;    
    subplot(4,2,6);
    plot(t,E_iAMPA(:,SDIM2),'b',t,E_iNMDA(:,SDIM2),'g');%,t,If_AMPA+If_NMDA,'k--'); 
    legend('E1<-iAMPA','E1<-iNMDA','Location','NorthWest'); ylabel('I(t)'); xlim(tlims);
    title(sprintf('E1(cell%g) inputs (AMPA,NMDA)',SDIM2)); ylim(ylims);   
    subplot(4,2,7); 
    plot(t,E_s_AMPA(:,SDIM1),'b',t,E_s_NMDA(:,SDIM1),'g',t,IfE_s_GABA(:,selIf),'r'); 
    legend('E1<-sAMPA','E1<-sNMDA','E1<-If1.sGABAa','Location','NorthWest'); ylabel('g(t)'); xlim(tlims);
    title(sprintf('E1(cell%g) synaptic conductance',SDIM1));
    ylim([0 1]);
    subplot(4,2,8);
    plot(t,E_s_AMPA(:,SDIM2),'b',t,E_s_NMDA(:,SDIM2),'g',t,IsE_s_GABA(:,selIs),'r'); 
    legend('E1<-sAMPA','E1<-sNMDA','E1<-Is1.sGABAa','Location','NorthWest'); ylabel('g(t)'); xlim(tlims);
    title(sprintf('E1(cell%g) synaptic conductance',SDIM2));
    ylim([0 1]);
    
    if 0 %IsIfsame
      % Synaptic conductances of interneurons
      var=strcmp('E1_Ifast1_AMPA_s',{data(If1).sensor_info.label});
      If_s_AMPA = squeeze(data(If1).epochs.data(var,:,1:nE)); % selIf));
      var=strcmp('E1_Ifast1_NMDA_sNMDApre',{data(If1).sensor_info.label});
      If_s_NMDA = squeeze(data(If1).epochs.data(var,:,1:nE)); % selIf));
      var=strcmp('Ifast1_Ifast1_GABAa_s',{data(If1).sensor_info.label});
      If_s_GABA = squeeze(data(If1).epochs.data(var,:,1:nI)); % selIf));
      var=strcmp('E1_Islow1_AMPA_s',{data(Is1).sensor_info.label});
      Is_s_AMPA = squeeze(data(Is1).epochs.data(var,:,1:nE)); % selIs));
      var=strcmp('Islow1_Islow1_GABAa_s',{data(Is1).sensor_info.label});
      Is_s_GABA = squeeze(data(Is1).epochs.data(var,:,1:nIs)); % selIs));
      % Synaptic currents to interneurons
      E_If_iAMPA=zeros(nt,nI);
      E_If_iNMDA=zeros(nt,nI);
      If_If_iGABA=zeros(nt,nI);
      E_Is_iAMPA=zeros(nt,nIs);
      Is_Is_iGABA=zeros(nt,nIs);
      for a=1:nt
        E_If_iAMPA(a,:) = E1_Ifast1_AMPA_ISYN(If_V(a,:)',If_s_AMPA(a,:)');
        E_If_iNMDA(a,:) = E1_Ifast1_NMDA_INMDA(If_V(a,:)',If_s_NMDA(a,:)');
        If_If_iGABA(a,:) = Ifast1_Ifast1_GABAa_ISYN(If_V(a,:)',If_s_GABA(a,:)');
        E_Is_iAMPA(a,:) = E1_Islow1_AMPA_ISYN(Is_V(a,:)',Is_s_AMPA(a,:)');
        Is_Is_iGABA(a,:) = Islow1_Islow1_GABAa_ISYN(Is_V(a,:)',Is_s_GABA(a,:)');
      end
      % select cells to plot
      If_AMPA = E_If_iAMPA(:,selIf);
      If_NMDA = E_If_iNMDA(:,selIf);
      If_GABA = If_If_iGABA(:,selIf);
      Is_AMPA = E_Is_iAMPA(:,selIs);
      Is_GABA = Is_Is_iGABA(:,selIs);
      % generate plots
      figure('position',[70 50 1800 850]);
      %tlims=[1 1.3];%1.5];%[1 2];
      mintime=max(StimOn,RuleOn);
      subplot(4,2,1); plot(t,E_V(:,SDIM1),'b',t,If_V(:,selIf),'r'); ylabel('V(t)'); 
      xlim(tlims); ylims=ylim; title('voltage');
      legend(sprintf('E1(cell%g)',SDIM1),sprintf('If1(cell%g)',selIf)); % legend('E','If');
      subplot(4,2,2); plot(t,E_V(:,SDIM2),'b',t,Is_V(:,selIs),'r'); ylabel('V(t)'); 
      xlim(tlims); ylim(ylims); title('voltage'); 
      legend(sprintf('E1(cell%g)',SDIM2),sprintf('Is1(cell%g)',selIs)); % legend('E','Is');
      subplot(4,2,3);
      plot(t,If_AMPA+If_NMDA,'k',t,If_GABA,'r'); 
      legend('If<-iAMPA+iNMDA','If<-iGABAa'); ylabel('I(t)'); 
      xlim(tlims); ylims=ylim; title(sprintf('If1(cell%g) inputs (E/I)',selIf));
      subplot(4,2,4);
      plot(t,Is_AMPA,'b',t,Is_GABA,'r'); 
      legend('Is<-iAMPA','Is<-iGABAa'); ylabel('I(t)'); xlim(tlims); ylim(ylims);
      title(sprintf('Is1(cell%g) inputs (E/I)',selIs)); 
      subplot(4,2,5); 
      plot(t,If_AMPA,'b',t,If_NMDA,'g');%,t,If_AMPA+If_NMDA,'k--'); 
      legend('If<-iAMPA','If<-iNMDA'); ylabel('I(t)'); xlim(tlims);
      title(sprintf('If1(cell%g) inputs (AMPA,NMDA)',selIf)); ylims=ylim;    
      subplot(4,2,6);
      plot(t,Is_AMPA,'b'); legend('Is<-AMPA'); ylabel('I(t)'); 
      xlim(tlims); ylim(ylims); title(sprintf('Is1(cell%g) inputs (AMPA,NMDA)',selIs)); 
      subplot(4,2,7); 
      plot(t,If_s_AMPA(:,selIf),'b',t,If_s_NMDA(:,selIf),'g',t,If_s_GABA(:,selIf),'r'); 
      legend('If<-sAMPA','If<-sNMDA','If<-sGABAa'); ylabel('g(t)'); xlim(tlims);
      title(sprintf('If1(cell%g) synaptic conductance',selIf));
      ylim([0 1]);
      subplot(4,2,8);
      plot(t,Is_s_AMPA(:,selIs),'b',t,Is_s_GABA(:,selIs),'r'); 
      legend('Is<-sAMPA','Is<-sGABAa'); ylabel('g(t)'); xlim(tlims);
      title(sprintf('Is1(cell%g) synaptic conductance',selIs)); 
      ylim([0 1]);
    end
    
  end  
  
  % analyze persistent activity
  if 0 % offset<tspan(2)
    %model=buildmodel(spec,'verbose',0);
    pset.p=model.model.parameters;
    timelimits=tspan;
    av=model.model.auxvars;
    for a=1:size(av,1)
      eval(sprintf('%s=%s;',av{a,1},av{a,2}));
    end
    fn=model.model.functions;
    for a=1:size(fn,1)
      eval(sprintf('%s=%s;',fn{a,1},fn{a,2}));
    end
    l={data(1).sensor_info.label};
    t=data(1).epochs.time;
    nt=length(t);
    E_V=squeeze(data(1).epochs.data(strcmp('E_V',l),:,:));
    E_E_AMPA_s=squeeze(data(1).epochs.data(strcmp('E_E_AMPA_s',l),:,1:nE));
    E_E_NMDA_sNMDApre=squeeze(data(1).epochs.data(strcmp('E_E_NMDA_sNMDApre',l),:,1:nE));
    I_E_GABAa_s=squeeze(data(1).epochs.data(strcmp('I_E_GABAa_s',l),:,1:nI));
    EE_NMDA=zeros(nt,nE);
    EE_AMPA=zeros(nt,nE);
    IE_GABA=zeros(nt,nE);
    for a=1:nt
      EE_NMDA(a,:)=E_E_NMDA_INMDA(E_V(a,:)',E_E_NMDA_sNMDApre(a,:)');
      EE_AMPA(a,:)=E_E_AMPA_ISYN(E_V(a,:)',E_E_AMPA_s(a,:)');
      IE_GABA(a,:)=I_E_GABAa_ISYN(E_V(a,:)',I_E_GABAa_s(a,:)');
    end
    [min(IE_GABA(:)) max(IE_GABA(:))]
    [min(EE_AMPA(:)) max(EE_AMPA(:))]
    [min(EE_NMDA(:)) max(EE_NMDA(:))]

%     isyn=EE_NMDA+EE_AMPA+IE_GABA;
%     figure; plot(t,isyn(:,1));
%     
%     figure; xlims=[.75 1.5];
%     subplot(2,1,1); plot(t,E_V); xlim(xlims)
%     subplot(2,1,2); plot(t,EE_NMDA,'b',t,max(EE_NMDA(:))+IE_GABA,'r'); xlim(xlims)

    figure; xlims=[1 1.1];%[.75 1.5];
    subplot(3,1,1); plot(t,E_V); xlim(xlims)
    subplot(3,1,2); plot(t,-EE_NMDA,'b',t,-(max(EE_NMDA(:))+IE_GABA),'r'); xlim(xlims)
    subplot(3,1,3); plot(t,max(0,-EE_NMDA-EE_AMPA-IE_GABA)); xlim(xlims); ylim([-1 20])
%     subplot(3,1,3); plot(t,-EE_NMDA-EE_AMPA-IE_GABA); xlim(xlims); ylim([-20 20])
    hold on
    plot(t,20+5*E_E_NMDA_sNMDApre,'b'); hline(20,'k');
    plot(t,25+5*I_E_GABAa_s,'r'); ylim([-1 30]); hline(25,'k');
%     E_E_NMDA_BMg
%     E_E_NMDA_NT
%     E_E_NMDA_INMDA
%     E_E_AMPA_ISYN
%     I_E_GABAa_ISYN    
%     E_V
%     E_E_AMPA_s
%     E_E_NMDA_sNMDApre
%     I_E_GABAa_s
  end
  
end % end loop over stim amps
end % end loop over stim freqs
toc(tstart)

if strcmp('off',visible_status)
  set(findobj('type','figure'),'visible','on')
end

if ~summary_plots || buildonly, return; end
    
nsim=nval1*nval2;

if nsim>1 && rulesims_flag
  
  % -----------------------------------------
  
  PAmu=mean(ROImuaPA,3);
  PAsd=std(ROImuaPA,[],3);
  f0mu=mean(ROImuaf0,3);
  f0sd=std(ROImuaf0,[],3);
  FRmu=mean(ROIrate,3);
  FRsd=std(ROIrate,[],3);
  
  showINs=1;
  x=val2; markersize=15; outwidth=3;
  nr=3;     % (f0,PA,FR)
  nc=ntois; % (trial1,trial2)
  rois=1:6; % (A,D,L,R)
  figure('position',[425 83 1066 860]);
  X=repmat(x,[length(rois) 1])';
  for toi=1:ntois
    subplot(nr,nc,toi) % f0
    mu=squeeze(f0mu(:,toi,1,:));
    e=squeeze(f0sd(:,toi,1,:))/nval1;
    i=rois(1); errorbar(x,mu(i,:),e(i,:),e(i,:),'b*--','markersize',markersize); hold on
    i=rois(2); errorbar(x,mu(i,:),e(i,:),e(i,:),'r*--','markersize',markersize);
    i=rois(3); errorbar(x,mu(i,:),e(i,:),e(i,:),'bo-','linewidth',outwidth,'markersize',markersize);
    i=rois(4); errorbar(x,mu(i,:),e(i,:),e(i,:),'ro-','linewidth',outwidth,'markersize',markersize);
    if showINs
      i=rois(5); errorbar(x,mu(i,:),e(i,:),e(i,:),'k-','markersize',markersize);
      i=rois(6); errorbar(x,mu(i,:),e(i,:),e(i,:),'g-','markersize',markersize);
    end
    if toi==1
      if showINs
        legend('A','D','L','R','If2','Is2'); 
      else
        legend('A','D','L','R'); 
      end
    end
    xlabel(var2); ylabel('f0'); title(sprintf('toi %g',toi));
    xlim([min(x) max(x)]); ylim([min(ROImuaf0(:)) max(ROImuaf0(:))]);
    subplot(nr,nc,toi+nc) % PA
    mu=squeeze(PAmu(:,toi,1,:));
    e=squeeze(PAsd(:,toi,1,:))/nval1;
    i=rois(1); errorbar(x,mu(i,:),e(i,:),e(i,:),'b*--','markersize',markersize); hold on
    i=rois(2); errorbar(x,mu(i,:),e(i,:),e(i,:),'r*--','markersize',markersize);
    i=rois(3); errorbar(x,mu(i,:),e(i,:),e(i,:),'bo-','linewidth',outwidth,'markersize',markersize);
    i=rois(4); errorbar(x,mu(i,:),e(i,:),e(i,:),'ro-','linewidth',outwidth,'markersize',markersize);
    if showINs
      i=rois(5); errorbar(x,mu(i,:),e(i,:),e(i,:),'k-','markersize',markersize);
      i=rois(6); errorbar(x,mu(i,:),e(i,:),e(i,:),'g-','markersize',markersize);
    end
    xlabel(var2); ylabel('PA'); %title(sprintf('toi %g',toi));
    xlim([min(x) max(x)]); ylim([min(ROImuaPA(:)) max(ROImuaPA(:))]);
    subplot(nr,nc,toi+nc*2) % FR
    mu=squeeze(FRmu(:,toi,1,:));
    e=squeeze(FRsd(:,toi,1,:))/nval1;
    i=rois(1); errorbar(x,mu(i,:),e(i,:),e(i,:),'b*--','markersize',markersize); hold on
    i=rois(2); errorbar(x,mu(i,:),e(i,:),e(i,:),'r*--','markersize',markersize);
    i=rois(3); errorbar(x,mu(i,:),e(i,:),e(i,:),'bo-','linewidth',outwidth,'markersize',markersize);
    i=rois(4); errorbar(x,mu(i,:),e(i,:),e(i,:),'ro-','linewidth',outwidth,'markersize',markersize);
    if showINs
      i=rois(5); errorbar(x,mu(i,:),e(i,:),e(i,:),'k-','markersize',markersize);
      i=rois(6); errorbar(x,mu(i,:),e(i,:),e(i,:),'g-','markersize',markersize);
    end
    xlabel(var2); ylabel('FR'); %title(sprintf('toi %g',toi));
    xlim([min(x) max(x)]); ylim([min(ROIrate(:)) max(ROIrate(:))]);
  end  
    
  
%   ROIrate=zeros(nrois,ntois,nval1,nval2);
%   ROImuaf0=zeros(nrois,ntois,nval1,nval2);
%   ROImuaPA=zeros(nrois,ntois,nval1,nval2);
%   ROImuaPmax=zeros(nrois,ntois,nval1,nval2);
%   ROIsuaf0=zeros(nrois,ntois,nval1,nval2);
%   ROIsuaPA=zeros(nrois,ntois,nval1,nval2);
%   ROIsuaPmax=zeros(nrois,ntois,nval1,nval2);
%   ROIlfpf0=zeros(nrois,ntois,nval1,nval2);
%   ROIlfpPA=zeros(nrois,ntois,nval1,nval2);
%   ROIlfpPmax=zeros(nrois,ntois,nval1,nval2);  
  
  % plot (ROI,TOI) stats: (alltois,allrois,allroid)
  % ...
  
%       tmin=alltois{i}(1);
%       tmax=alltois{i}(2);
%       t1=nearest(t,tmin); % sec, TOI start time
%       t2=nearest(t,tmax); % sec, TOI stop time
%       tdur=tmax-tmin;
%       for j=1:nrois
%         id=allroid{j};  % population index
%         cellsel=allrois{j}; % cell indices
%         % ROI data
%         popname=spec.nodes(id).label;
%           title(sprintf('%s(%g-%g) [%g-%gs]',popname,cellsel(1),cellsel(end),tmin,tmax));

  if 0
    scaleup=10;
    figure
    subplot(2,1,1);
    roi=3; % L
    toi=1; % trial 1
    f01=squeeze(ROImuaf0(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
    pa1=squeeze(ROImuaPA(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
    scatter(val2,f01,pa1*scaleup,'b','o'); hold on
    roi=4; % R
    toi=1; % trial 1
    f02=squeeze(ROImuaf0(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
    pa2=squeeze(ROImuaPA(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
    scatter(val2,f02,pa2*scaleup,'g','o');
    legend('L','R'); xlabel(var2); ylabel('f0 (size=PA)');
    plot(val2,f01,'b');
    plot(val2,f02,'g');
    title('trial 1');
    if ntois>1
      subplot(2,1,2);
      roi=3; % L
      toi=2; % trial 1
      f01=squeeze(ROImuaf0(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
      pa1=squeeze(ROImuaPA(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
      scatter(val2,f01,pa1*scaleup,'b','o'); hold on
      roi=4; % R
      toi=2; % trial 1
      f02=squeeze(ROImuaf0(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
      pa2=squeeze(ROImuaPA(roi,toi,1,:)); % L, trial 1, (~), {gIE2}
      scatter(val2,f02,pa2*scaleup,'g','o');
      legend('L','R'); xlabel(var2); ylabel('f0 (size=PA)');
      plot(val2,f01,'b');
      plot(val2,f02,'g');
      title('trial 2');    
    end
  end  
  
  if 0
    x=val2;
    nr=2;     % (PA,FR)
    nc=ntois; % (trial1,trial2)
    rois=1:4; % (A,D,L,R)
    figure
    for toi=1:ntois
      subplot(nr,nc,toi) % PA
      y=squeeze(ROImuaPA(rois(1),toi,:,:));
      plot(x,y,'b--'); hold on
      y=squeeze(ROImuaPA(rois(2),toi,:,:));
      plot(x,y,'r--');
      y=squeeze(ROImuaPA(rois(3),toi,:,:));
      plot(x,y,'b-','linewidth',2);
      y=squeeze(ROImuaPA(rois(4),toi,:,:));
      plot(x,y,'r-','linewidth',2);
      xlabel(var2); ylabel('PA'); title(sprintf('toi %g',toi));
      xlim([min(x) max(x)]); ylim([min(ROImuaPA(:)) max(ROImuaPA(:))]);
      subplot(nr,nc,toi+nc) % FR
      y=squeeze(ROIrate(rois(1),toi,:,:));
      plot(x,y,'b--'); hold on
      y=squeeze(ROIrate(rois(2),toi,:,:));
      plot(x,y,'r--');
      y=squeeze(ROIrate(rois(3),toi,:,:));
      plot(x,y,'b-','linewidth',2);
      y=squeeze(ROIrate(rois(4),toi,:,:));
      plot(x,y,'r-','linewidth',2);
      xlabel(var2); ylabel('FR'); title(sprintf('toi %g',toi));
      xlim([min(x) max(x)]); ylim([min(ROIrate(:)) max(ROIrate(:))]);
    end
  end
  % -----------------------------------------
  
  
  showbounds=1; LB=9; UB=17;
  figure('position',[680 40 560 930]);
  subplot(3,1,1);
  muRdp=mean(ruleRdp,1);
  se=std(ruleRdp,[],1)/sqrt(nval1);
  errorbar(val2,muRdp,se,'bo-'); hold on; 
  muRdr=mean(ruleRdr,1);
  se=std(ruleRdr,[],1)/sqrt(nval1);
  errorbar(val2,muRdr,se,'ro-'); axis tight; hline(0,'k');
  legend('dPA','dFR'); xlabel(var2); ylabel('output correlation');
  if showbounds, vline(LB,'color','k','linestyle','--'); vline(UB,'color','k','linestyle','--'); end
  subplot(3,1,2);
  plot(val2,muRdp-muRdr,'bo-'); xlabel(var2); ylabel('dPA-dFR'); axis tight
  hline(0,'color','k','linestyle','--');
%   mu=mean(rulepdp,1);
%   se=std(rulepdp,[],1)/sqrt(nval1);
%   errorbar(val2,mu,se,'bo-'); hold on
%   mu=mean(rulepdr,1);
%   se=std(rulepdr,[],1)/sqrt(nval1);
%   errorbar(val2,mu,se,'ro-');
%   legend('dPA','dFR'); xlabel(var2); ylabel('p-values (for output prediction)');
%   hline(0.05,'color','k','linestyle','--');
%   if showbounds, vline(LB,'color','k','linestyle','--'); vline(UB,'color','k','linestyle','--'); end
  subplot(3,1,3);
  mu=mean(rulef0,1);
  se=std(rulef0,[],1)/sqrt(nval1);
  errorbar(val2,mu,se,'bo-'); axis tight
  xlabel(var2); ylabel('LFP freq (for readout)');
  if showbounds
    vline(LB,'color','k','linestyle','--'); vline(UB,'color','k','linestyle','--');
    hline(28,'color','k','linestyle','--');
    hline(17,'color','k','linestyle','--');
  end

  if nval1>1 && nval2>1
    figure('position',[680 560 680 420]);
    roi1=1; roi2=2;
    x=val1; y=val2;
    subplot(2,2,1);
    z = ruleRdp-ruleRdr;
    contour(x,y,z',15); colorbar; colormap(jet);
    xlabel(var1); ylabel(var2); title('Corr(Pow)-Corr(FR)');
    cmax=max(abs(caxis)); caxis([-cmax cmax]); axis square;
    subplot(2,2,2);
    z = squeeze(ruleRin(roi1,:,:)-ruleRin(roi2,:,:));
    contour(x,y,z',15); colorbar; colormap(jet);
    xlabel(var1); ylabel(var2); title('<FR(A)-FR(D)> (Hz)');
    cmax=max(abs(caxis)); caxis([-cmax cmax]); axis square;
    subplot(2,2,3);
    z = rulef0;
    contour(x,y,z',15); colorbar; colormap(jet);
    xlabel(var1); ylabel(var2); title('LFP freq (Hz)');
    caxis([0 max(caxis)]); axis square;
    subplot(2,2,4);
    io1=(ruleRout(roi1,:,:)./ruleRin(roi1,:,:));
    io2=(ruleRout(roi2,:,:)./ruleRin(roi2,:,:));
    z=squeeze(io1-io2);
    contour(x,y,z',15); colorbar; colormap(jet);
    xlabel(var1); ylabel(var2); title('gain(A->L) - gain(D->R)');
    cmax=max(abs(caxis)); caxis([-cmax cmax]); axis square;
  end
  
  dR=ruleRdp-ruleRdr;
  dRmax=max(dR,[],2);
  dRmaxInd=zeros(1,nval1);
  fRmax=zeros(1,nval1);
  IRmax=zeros(1,nval1);
  f0=rulef0;
  for i=1:nval1
    ind=find(dRmax(i)==dR(i,:));
    if ~isempty(ind)
      dRmaxInd(i)=ind;
      IRmax(i)=y(ind);
      fRmax(i)=f0(i,ind);
    end
  end
  if 0
    figure
    subplot(3,1,1)
    plot(x,fRmax,'o-'); xlabel(var1); ylabel('f0 at dRmax');
    subplot(3,1,2);
    plot(x,IRmax,'o-'); xlabel(var1); ylabel('Ie at dRmax');
    subplot(3,1,3);
    plot(fRmax,dRmax,'o-'); xlabel('f0 at dRmax'); ylabel('dRmax');
  end
  
  figure('position',[680 70 560 900]);
  dRlims=[min(dR(:)) max(dR(:))];
  f0lims=[min(f0(:)) max(f0(:))];
  for i=1:nval1
    val=val1(i);%13;
    ind=find(x==val);
    subplot(nval1,3,1+3*(i-1));
    plot(y,dR(ind,:),'bo-'); axis([min(y) max(y) dRlims]); hline(0,'k');
    xlabel(var2); ylabel('dR'); title(sprintf('%s=%g',var1,val));
    subplot(nval1,3,2+3*(i-1));
    f=f0(ind,:); %f(4)=20;
    plot(y,f,'bo-'); axis([min(y) max(y) f0lims]); hline(mean(f0(:)),'k');
    xlabel(var2); ylabel('f0'); title(sprintf('%s=%g',var1,val));
    subplot(nval1,3,3+3*(i-1));  
    plot(f,dR(ind,:),'bo-'); axis([f0lims dRlims]); hline(0,'k');
    xlabel('f0'); ylabel('dR'); title(sprintf('%s=%g',var1,val));
  end
  
%   figure
% %   z = squeeze(ruleRin(roi1,:,:));
% %   contour(x,y,z',15); colorbar; colormap(jet);
% %   xlabel(var1); ylabel(var2); title('FR(A)');
%   imagesc(x,y,dR'); colorbar;
%   xlabel(var1); ylabel(var2); title('Corr(Pow)-Corr(FR)'); axis xy
  
  
  dCorr=muRdp-muRdr;
  
  rout=squeeze(mean(ruleRout,2));
  rin=squeeze(mean(ruleRin,2));
  Cin=squeeze(mean(ruleCin,2));
  roi1=1; roi2=2;
  dRout=squeeze(mean(ruleRout(roi1,:,:)-ruleRout(roi2,:,:),2));
  dRin=squeeze(mean(ruleRin(roi1,:,:)-ruleRin(roi2,:,:),2));
  dCin=squeeze(mean(ruleCin(roi1,:,:)-ruleCin(roi2,:,:),2));
  dTransfer=squeeze(mean((ruleRout(roi1,:,:)./ruleRin(roi1,:,:))-(ruleRout(roi2,:,:)./ruleRin(roi2,:,:)),2));
  
  if 0
    figure
    subplot(2,2,1);
    plot(dRin,dCorr,'o-'); hline(0,'k'); vline(0,'k');
    title(sprintf('sync-specific routing vs differential rates (avg over %s)',var1));
    xlabel('<FR(A)-FR(D)>'); ylabel('<Corr(dPA)-Corr(dFR)>')
    subplot(2,2,2);
    plot(dCin,dCorr,'o-'); hline(0,'k'); vline(0,'k');
    title('sync-specific routing vs differential coherence (avg)');
    xlabel('<coh(A)-coh(D)>'); ylabel('<Corr(dPA)-Corr(dFR)>');
    subplot(2,2,3);
    plot(dRin,dTransfer,'o-'); hline(0,'k'); vline(0,'k');
    title('differential gain vs differential rates (avg)');
    xlabel('<FR(A)-FR(D)>'); ylabel('<gain(A)-gain(D)>');
    subplot(2,2,4);
    plot(dCin,dTransfer,'o-'); hline(0,'k'); vline(0,'k');
    title('differential gain vs differential coherence (avg)');
    xlabel('<coh(A)-coh(D)>'); ylabel('<gain(A)-gain(D)>');
  end  
  
  ro1=squeeze(ruleRout(1,:,:));
  ro2=squeeze(ruleRout(2,:,:));
  dro=ro1-ro2;
  mudro=nanmean(dro,1);
  sedro=std(dro,[],1)/sqrt(nval1);
  muro1=nanmean(ro1,1);
  sero1=std(ro1,[],1)/sqrt(nval1);
  muro2=nanmean(ro2,1);
  sero2=std(ro2,[],1)/sqrt(nval1);  
  ri1=squeeze(ruleRin(1,:,:));
  ri2=squeeze(ruleRin(2,:,:));
  dri=ri1-ri2;
  mudri=nanmean(dri,1);
  sedri=std(dri,[],1)/sqrt(nval1);
  muri1=nanmean(ri1,1);
  seri1=std(ri1,[],1)/sqrt(nval1);
  muri2=nanmean(ri2,1);
  seri2=std(ri2,[],1)/sqrt(nval1);  
  
  reldro=(ro1-ro2)./(ro1+ro2);
  mureldro=mean(reldro,1);
  sereldro=std(reldro,[],1)/sqrt(nval1);
  reldri=(ri1-ri2)./(ri1+ri2);
  mureldri=mean(reldri,1);
  sereldri=std(reldri,[],1)/sqrt(nval1);
  
  if 1
    figure('position',[400 100 1030 790]);
    subplot(2,2,1);
    plot(dri,dro,'o-'); axis tight; hold on;
    str=cellfun(@num2str,num2cell(val2),'uni',0); legend(str{:},'location','southeast');
    plot(mudri,mudro,'ko-','linewidth',5); title('differential input vs output');
    xlabel('<A-D> (Hz)'); ylabel('<L-R> (Hz)');
    hline(0,'color','k','linestyle','--'); vline(0,'color','k','linestyle','--');
    subplot(2,2,2);
    imagesc(1:nval2,1:nval1,dR); %colormap(1-gray); % colorbar; 
    xlabel(var2); ylabel(var1); title('Corr(Pow)-Corr(FR)'); axis xy    
    set(gca,'xtick',1:nval2,'xticklabel',val2); 
    set(gca,'ytick',1:nval1,'yticklabel',val1)    
    subplot(2,2,3);
    plot(val2,reldro); axis tight; hold on
    if nval1>1 && nval2>1
      errorbar(val2,mureldro,sereldro,'ko-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    end
    xlabel(var2); ylabel('<(L-R)/(L+R)>'); title('relative differential output');          
    str=cellfun(@num2str,num2cell(val1),'uni',0); legend(str{:},'location','northeast');
%     plot(val2,dro); axis tight; hold on
%     errorbar(val2,mudro,sedro,'go-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
%     xlabel(var2); ylabel('<L-R> (Hz)'); title('differential output');    
    subplot(2,2,4);
    mu=mean(dR,1);
    se=std(dR,[],1);
    plot(val2,dR); axis tight; hold on
    if nval1>1 && nval2>1
      errorbar(val2,mu,se,'ko-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    end
    xlabel(var2); ylabel('Corr(Pow)-Corr(FR)');
    str=cellfun(@num2str,num2cell(val1),'uni',0); legend(str{:},'location','northeast');
  end  
  
  if nval1>1 && nval2>1
  %   figure('position',[680 230 760 700]);
  %   subplot(3,2,2);
  %   errorbar(val2,muro1,sero1,'bo-'); axis tight; hold on
  %   errorbar(val2,muro2,sero2,'ro-'); axis tight    
  %   %plot(val2,muro1,'bo-',val2,muro2,'ro-'); axis tight
  %   xlabel(var2); ylabel('FR (Hz)'); legend('Left','Right'); title('output firing rates');
  %   subplot(3,2,4);
  %   errorbar(val2,mudro,sedro,'go-'); axis tight; hline(0,'color','k','linestyle','--');
  %   xlabel(var2); ylabel('<r(L)-r(R)> (Hz)'); title('differential output');
  %   subplot(3,2,1);
  %   errorbar(val2,muri1,seri1,'bo-'); axis tight; hold on
  %   errorbar(val2,muri2,seri2,'ro-'); axis tight    
  %   xlabel(var2); ylabel('FR (Hz)'); legend('CatA','CatD'); title('input firing rates');
  %   subplot(3,2,3);
  %   errorbar(val2,mudri,sedri,'go-'); axis tight; hline(0,'color','k','linestyle','--');
  %   xlabel(var2); ylabel('<r(A)-r(D)> (Hz)'); title('differential input');
  %   subplot(3,2,5);
  %   plot(val2,mudri./mudro,'go-'); title('differential input/output');
  %   xlabel(var2); ylabel('<r(A)-r(D)>/<r(L)-r(R)>');
  %   hline(0,'color','k','linestyle','--'); vline(0,'color','k','linestyle','--');
  %   subplot(3,2,6);
  %   plot(mudri,mudro,'go-'); title('differential input vs output');
  %   xlabel('<r(A)-r(D)> (Hz)'); ylabel('<r(L)-r(R)> (Hz)');
  %   hline(0,'color','k','linestyle','--'); vline(0,'color','k','linestyle','--');
    figure('position',[500 10 900 950]);
    subplot(4,2,2);
    errorbar(val2,muro1,sero1,'bo-','linewidth',3); axis tight; hold on
    errorbar(val2,muro2,sero2,'ro-','linewidth',3); axis tight    
    plot(val2,ro1,'o-',val2,ro2,'^-'); axis tight; hold on
    %plot(val2,muro1,'bo-',val2,muro2,'ro-'); axis tight
    xlabel(var2); ylabel('FR (Hz)'); legend('Left','Right','location','southeast'); title('output firing rates');
    subplot(4,2,4);
    plot(val2,dro); axis tight; hold on
    errorbar(val2,mudro,sedro,'go-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    xlabel(var2); ylabel('<L-R> (Hz)'); title('differential output');
    subplot(4,2,6)
    plot(val2,reldro); axis tight; hold on
    errorbar(val2,mureldro,sereldro,'ko-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    xlabel(var2); ylabel('<(L-R)/(L+R)>'); title('relative differential output');      
    subplot(4,2,1);
    errorbar(val2,muri1,seri1,'bo-','linewidth',3); axis tight; hold on
    errorbar(val2,muri2,seri2,'ro-','linewidth',3); axis tight    
    plot(val2,ri1,'o-',val2,ri2,'^-'); axis tight; hold on
    xlabel(var2); ylabel('FR (Hz)'); legend('CatA','CatD','location','southeast'); title('input firing rates');
    subplot(4,2,3);
    plot(val2,dri); axis tight; hold on
    errorbar(val2,mudri,sedri,'go-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    xlabel(var2); ylabel('<A-D> (Hz)'); title('differential input');
    subplot(4,2,5)
    plot(val2,reldri); axis tight; hold on
    errorbar(val2,mureldri,sereldri,'go-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    xlabel(var2); ylabel('<(A-D)/(A+D)>'); title('relative differential input');      
    subplot(4,2,7);
    plot(val2,dri./dro); axis tight; hold on
    plot(val2,mudri./mudro,'go-','linewidth',3); title('differential input/output');
    xlabel(var2); ylabel('<A-D>/<L-R>');
    hline(0,'color','k','linestyle','--'); vline(0,'color','k','linestyle','--');
    subplot(4,2,8);
    plot(val2,reldri./reldro); axis tight; hold on
    plot(val2,mureldri./mureldro,'go-','linewidth',3); title('relative differential input/output');
    xlabel(var2); ylabel('<(A-D)/(A+D)>/<(L-R)/(L+R)>');
    hline(0,'color','k','linestyle','--'); vline(0,'color','k','linestyle','--');
  %   plot(dri,dro,'o-'); axis tight; hold on;
  %   str=cellfun(@num2str,num2cell(val2),'uni',0); legend(str{:},'location','southeast');
  % %   plot(dri(:),dro(:),'o'); axis tight; hold on;
  %   plot(mudri,mudro,'ko-','linewidth',5); title('differential input vs output');
  %   xlabel('<A-D> (Hz)'); ylabel('<r(L)-R> (Hz)');
  %   hline(0,'color','k','linestyle','--'); vline(0,'color','k','linestyle','--');
  end
    
  if 0
    figure % plot (L-R)/(L+R) vs (L+R)/2  (relative increase vs mean output)
    plot(mudro,reldro); axis tight; hold on
    errorbar(mudro,mureldro,sereldro,'ko-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    xlabel('r(out)=(L+R)/2'); ylabel('<(L-R)/(L+R)>'); title('relative differential output');          
    str=cellfun(@num2str,num2cell(val1),'uni',0); legend(str{:},'location','northeast');
  end
    
  if 0
    figure
    plot(val2,reldro); axis tight; hold on
    errorbar(val2,mureldro,sereldro,'go-','linewidth',3); axis tight; hline(0,'color','k','linestyle','--');
    xlabel(var2); ylabel('<(r(L)-r(R))/(r(L)+r(R))> (Hz)'); title('relative differential output');    
  end
  
  if 0
    dCorr=ruleRdp-ruleRdr;
    dRout=squeeze(ruleRout(roi1,:,:)-ruleRout(roi2,:,:));
    dRin=squeeze(ruleRin(roi1,:,:)-ruleRin(roi2,:,:));
    dCin=squeeze(ruleCin(roi1,:,:)-ruleCin(roi2,:,:));  
    dTransfer=squeeze((ruleRout(roi1,:,:)./ruleRin(roi1,:,:))-(ruleRout(roi2,:,:)./ruleRin(roi2,:,:)));
    dRout=dRout(:);
    dRin=dRin(:);
    dCin=dCin(:);
    dTransfer=dTransfer(:);
    dCorr=dCorr(:);
    figure
    subplot(2,2,1);
    plot(dRin,dCorr,'o-'); 
    title('sync-specific routing vs differential rates');
    xlabel('FR(A)-FR(D)'); ylabel('Corr(dPA)-Corr(dFR)')
    subplot(2,2,2);
    plot(dCin,dCorr,'o-');
    title('sync-specific routing vs differential coherence');
    xlabel('coh(A)-coh(D)'); ylabel('Corr(dPA)-Corr(dFR)');
    subplot(2,2,3);
    plot(dRin,dTransfer,'o-');
    title('differential gain vs differential rates');
    xlabel('FR(A)-FR(D)'); ylabel('gain(A)-gain(D)');
    subplot(2,2,4);
    plot(dCin,dTransfer,'o-');
    title('differential gain vs differential coherence');
    xlabel('coh(A)-coh(D)'); ylabel('gain(A)-gain(D)');
  end
  
end

if nsim>1
  if 0
    figure
    cnt=0;
    for i=1:nval2
      for j=1:nval1
        cnt=cnt+1;
        subplot(nsim,1,cnt);
        bar(edges,muaisirates(:,cnt),'histc'); 
        xlabel('1/muaISIs [Hz]');
        title(sprintf('%s=%g, %s=%g',var1,val1(j),var2,val2(i)));
      end
    end
    figure
    cnt=0;
    for i=1:nval2
      for j=1:nval1
        cnt=cnt+1;
        subplot(nsim,1,cnt);
        bar(edges,cellisirates(:,cnt),'histc'); 
        xlabel('1/cellnetISIs [Hz]');
        title(sprintf('%s=%g, %s=%g',var1,val1(j),var2,val2(i)));
      end
    end
  end
  
  ylims=[min(spkcoh(:))-.1 max(spkcoh(:))+.1];
  figure('position',[70 150 1850 800])
  subplot(2,1,1);
  cnt=0; xticks=[];
  for i=1:nval1
    coh11 = squeeze(spkcoh(1,1,i,:)); % E1-E1, E2-E2, E1-E2
    coh22 = squeeze(spkcoh(2,2,i,:));
    coh12 = squeeze(spkcoh(1,2,i,:));
    xs=cnt+(1:nval2);
    plot(xs,coh11,'bo-',xs,coh22,'ro-',xs,coh12,'go-'); 
    ylim(ylims); hold on
    xticks=[xticks cnt+[1 nval2]];
    cnt=cnt+nval2; vline(cnt,'k');    
    if i==1
      text(1,mean(ylims),sprintf('%s=%g',var1,val1(i)));
    else
      text(cnt-nval2/2,mean(ylims),sprintf('%g',val1(i)));
    end
  end    
  set(gca,'xtick',xticks,'xticklabel',[val2(1) val2(end)]);
  legend('E1-E1','E2-E2','E1-E2','Location','NorthEast');
  xlabel(sprintf('%s',var2)); ylabel('spkcoh');
  subplot(2,1,2);
  cnt=0; xticks=[];
  for i=1:nval2
    coh11 = squeeze(spkcoh(1,1,:,i)); % E1-E1, E2-E2, E1-E2
    coh22 = squeeze(spkcoh(2,2,:,i));
    coh12 = squeeze(spkcoh(1,2,:,i));
    xs=cnt+(1:nval1);
    plot(xs,coh11,'bo-',xs,coh22,'ro-',xs,coh12,'go-'); 
    ylim(ylims); hold on
    xticks=[xticks cnt+[1 nval1]];
    cnt=cnt+nval1; vline(cnt,'k');    
    if i==1
      text(1,mean(ylims),sprintf('%s=%g',var2,val2(i)));
    else
      text(cnt-nval1/2,mean(ylims),sprintf('%g',val2(i)));
    end
  end    
  set(gca,'xtick',xticks,'xticklabel',[val1(1) val1(end)]);
  legend('E1-E1','E2-E2','E1-E2','Location','NorthEast');
  xlabel(sprintf('%s',var1)); ylabel('spkcoh');
  
%   coh11 = spkcoh(1,1,:,:); coh11=coh11(:); % E1-E1, E2-E2, E1-E2
%   coh22 = spkcoh(2,2,:,:); coh22=coh22(:);
%   coh12 = spkcoh(1,2,:,:); coh12=coh12(:);
%   nl=length(coh11);
%   plot(1:nl,coh11,'bo-',1:nl,coh22,'ro-',1:nl,coh12,'go-'); 
%   legend('E1-E1','E2-E2','E1-E2');
%   xlabel('sim'); ylabel('spkcoh');
%   ylim([min(spkcoh(:))-.1 max(spkcoh(:))+.1]);
end

if nval1>1 && nval2==1
  sel=find(muaOscFreqs(1,:)>5);
  ff=muaOscFreqs(1,sel);
  rr=AvgSpikeRates(1,sel);
  ss=val1(sel);
  try
    [pks,locs]=findpeaks(diff(muaAreaPowers(1,sel)));
    DenseFreq = ff(locs(max(pks)==pks));
    RateShift = rr(locs(max(pks)==pks));
    StimShift = ss(locs(max(pks)==pks));
  catch
    DenseFreq = nan;
    RateShift = nan;
    StimShift = nan;
  end
  figure('position',[100 100 1650 400]);
  subplot(1,4,1); plot(val1,muaOscFreqs(1,:),'b*-',val1,muaOscFreqs(2,:),'ro--'); xlabel(var1); ylabel('MUA freq [Hz]'); legend('E','I');
  axis tight; 
  if ~isempty(DenseFreq),hline(DenseFreq,'b'); vline(StimShift,'b'); end
  if InputType==1
    title(sprintf('EsomaType=%g,InputType=%g,tauI=%g,tauI=%g',EsomaType(1),InputType,tauI,tauI));
  else
    title(sprintf('EsomaType=%g,InputType=%g,tauI=%g,tauI=%g',EsomaType(1),InputType,tauI,tauI));
  end
  subplot(1,4,3); plot(muaOscFreqs(1,:),muaAreaPowers(1,:),'b*-',muaOscFreqs(2,:),muaAreaPowers(2,:),'ro--');axis tight; 
  xlabel('MUA freq [Hz]'); ylabel('Rhythm Power Area'); legend('E','I'); 
  if ~isempty(DenseFreq),vline(DenseFreq,'b'); end
  subplot(1,4,2); plot(val1,AvgSpikeRates(1,:),'b*-',val1,AvgSpikeRates(2,:),'ro--'); 
  xlabel(var1); ylabel('Avg Cell Spike Rate [Hz]'); legend('E','I'); axis tight;
  if ~isempty(DenseFreq),vline(StimShift,'b'); hline(RateShift,'b'); end
  subplot(1,4,4); plot(muaOscFreqs(1,:),AvgSpikeRates(1,:),'b*-',muaOscFreqs(2,:),AvgSpikeRates(2,:),'ro--',1:max(muaOscFreqs(:)),1:max(muaOscFreqs(:)),'k--');axis tight; 
  xlabel('MUA freq [Hz]'); ylabel('Avg Cell Spike Rate [Hz]'); legend('E','I','diag'); 
  if ~isempty(DenseFreq),vline(DenseFreq,'b'); hline(RateShift,'b'); end
  if nval1<10
    figure
    tmp=MUAs(:,1,:); ymin1=min(tmp(:)); ymax1=max(tmp(:));
    %tmp=MUApows(f<Fmax,1,:); ymin2=min(tmp(:)); ymax2=max(tmp(:));
    tmp=allrates(:,1,:); ymin3=min(tmp(:)); ymax3=max(tmp(:));
    nc=3;
    for stim=1:nval1
      subplot(nval1,nc,nc*(stim-1)+1);
      plot(t,MUAs(:,1,stim)); ylim([ymin1 ymax1]);
      xlabel('t'); ylabel('MUA');
      subplot(nval1,nc,nc*(stim-1)+2);
      plot(muaf,log10(MUApows(:,1,stim))); 
      xlim([0 Fmax]); xlabel('f'); ylabel('p'); ylim([-5 5]);
      vline(muaOscFreqs(1,stim),'k');
      subplot(nval1,nc,nc*(stim-1)+3);
      plot(tmins,allrates(:,1,stim)); ylim([ymin3 ymax3]);
      xlim([0 max(tmins)]); xlabel('t'); ylabel('r');
    end
  end
elseif nval1==1 && nval2>1
  x=val2;
  AvgSpikeRates=squeeze(AvgSpikeRates);
  muaOscFreqs=squeeze(muaOscFreqs);
  muaAreaPowers=squeeze(muaAreaPowers);
  rateOscFreqs=squeeze(rateOscFreqs);
  rateAreaPowers=squeeze(rateAreaPowers);
  MUAs=squeeze(MUAs);
  MUApows=squeeze(MUApows);
  allrates=squeeze(allrates);
  ratepows=squeeze(ratepows);
  NumSpikes=squeeze(NumSpikes);
  figure('position',[100 100 1650 500]);
  subplot(2,4,1); plot(x,muaOscFreqs(1,:),'b*-',x,muaOscFreqs(2,:),'ro-','linewidth',3); 
  xlabel(var2); ylabel('MUA freq [Hz]'); axis tight; 
  hold on; plot(min(x):max(x),min(x):max(x),'k--'); legend('E','I','diag');
  if InputType==1
    title(sprintf('EsomaType=%g,InputType=%g,tauI=%g,tauI=%g',EsomaType(1),InputType,tauI,tauI));
  else
    title(sprintf('EsomaType=%g,InputType=%g,tauI=%g,tauI=%g',EsomaType(1),InputType,tauI,tauI));
  end
  subplot(2,4,2); plot(x,AvgSpikeRates(1,:),'b*-',x,AvgSpikeRates(2,:),'ro-','linewidth',3);
  xlabel(var2); ylabel('Avg Spike Rate [Hz]'); axis tight;
  hold on; plot(min(x):max(x),min(x):max(x),'k--'); legend('E','I','diag');
  subplot(2,4,3); plot(muaOscFreqs(1,:),AvgSpikeRates(1,:),'b*-',muaOscFreqs(2,:),AvgSpikeRates(2,:),'ro-','linewidth',2);axis tight;
  xlabel('MUA freq [Hz]'); ylabel('Avg Spike Rate [Hz]'); 
  hold on; plot(1:max(muaOscFreqs(:)),1:max(muaOscFreqs(:)),'k--'); legend('E','I','diag');
  subplot(2,4,4); plot(muaOscFreqs(1,:),muaAreaPowers(1,:),'b*-',muaOscFreqs(2,:),muaAreaPowers(2,:),'ro-');axis tight; 
  xlabel('MUA freq [Hz]'); ylabel('MUA Power Area'); legend('E','I');
  subplot(2,4,5); plot(x,NumSpikes(1,:),'b*-',x,NumSpikes(2,:),'ro-');
  xlabel(var2); ylabel('avg #spikes per cell'); legend('E','I');
  subplot(2,4,6); plot(x,rateOscFreqs(1,:),'b*-',x,rateOscFreqs(2,:),'ro-'); 
  xlabel(var2); ylabel('rate spectral freq [Hz]'); axis tight;
  hold on; plot(min(x):max(x),min(x):max(x),'k--'); legend('E','I','diag');
  subplot(2,4,7); plot(rateOscFreqs(1,:),rateAreaPowers(1,:),'b*-',rateOscFreqs(2,:),rateAreaPowers(2,:),'ro-');axis tight; 
  xlabel('rate spectral freq [Hz]'); ylabel('rate power Area'); legend('E','I');
  
  if nval2<10
    figure
    tmp=MUAs(:,1,:); ymin1=min(tmp(:)); ymax1=max(tmp(:));
    tmp=allrates(:,1,:); ymin3=min(tmp(:)); ymax3=max(tmp(:));
    nc=1;
    for freq=1:nval2
      subplot(nval2,nc,nc*(freq-1)+1);
      plot(t,MUAs(:,1,freq)); ylim([ymin1 ymax1]);
      xlabel('t'); ylabel('MUA'); title(sprintf('f_0=%gHz | stim=%g',val2(freq),val1));
      if nc==1, continue; end
      subplot(nval2,nc,nc*(freq-1)+2);
      plot(muaf,log10(MUApows(:,1,freq))); 
      xlim([0 Fmax]); xlabel('f'); ylabel('p'); ylim([-5 5]);
      vline(muaOscFreqs(1,freq),'k');
      if nc==2, continue; end
      subplot(nval2,nc,nc*(freq-1)+3);
      plot(tmins,allrates(:,1,freq)); ylim([ymin3 ymax3]);
      xlim([0 max(tmins)]); xlabel('t'); ylabel('r');
      if nc==3, continue; end
      subplot(nval2,nc,nc*(freq-1)+4);
      plot(ratef,log10(ratepows(:,1,freq))); 
      xlim([0 min(max(ratef),Fmax)]); xlabel('f'); ylabel('p'); ylim([-5 5]);
      vline(rateOscFreqs(1,freq),'k');
    end
  end
elseif nval1>1 && nval2>1
  % resonance and frequency transformation
  figure('position',[100 10 900 950]);
  nc=3; linewidth=2;
  minasr=min(AvgSpikeRates(:)); maxasr=max(AvgSpikeRates(:));
  minmof=min(muaOscFreqs(:)); maxmof=max(muaOscFreqs(:));
  for stim=1:nval1
    x=val2;
    mof=squeeze(muaOscFreqs(:,stim,:));
    asr=squeeze(AvgSpikeRates(:,stim,:));    
    subplot(nval1,nc,nc*(stim-1)+1); plot(x,mof(1,:),'b*-',x,mof(2,:),'ro-','linewidth',linewidth);
    xlabel(var2); if stim==1, ylabel('MUA freq [Hz]'); end
    axis tight;
    hold on; plot(min(x):max(x),min(x):max(x),'k--');
    title(sprintf('%s=%g',var1,val1(stim)));
    subplot(nval1,nc,nc*(stim-1)+2); plot(x,asr(1,:),'b*-',x,asr(2,:),'ro-','linewidth',linewidth);
    xlabel(var2); axis tight;
    hold on; plot(min(x):max(x),min(x):max(x),'k--');
    if stim==1, ylabel('Avg Spike Rate'); end
    subplot(nval1,nc,nc*(stim-1)+3); plot(mof(1,:),asr(1,:),'b*-',mof(2,:),asr(2,:),'ro-','linewidth',linewidth);
    xlabel('MUA freq [Hz]'); if stim==1, ylabel('Avg Spike Rate [Hz]'); end
    axis tight; 
    hold on; plot(1:max(mof(:)),1:max(mof(:)),'k--');
    if stim==1
      title(sprintf('EsomaType=%g,InputType=%g,tauI=%g,tauI=%g',EsomaType(1),InputType,tauI,tauI));
      legend('E','I','diag');
    else
      title(sprintf('%s=%g',var1,val1(stim)));    
    end
  end
  
  if 1
    coh11 = squeeze(spkcoh(1,1,:,:)); % E1-E1, E2-E2, E1-E2
    coh22 = squeeze(spkcoh(2,2,:,:));
    coh12 = squeeze(spkcoh(1,2,:,:));
    rE1 = squeeze(AvgSpikeRates(1,:,:));
    rE2 = squeeze(AvgSpikeRates(2,:,:));
    rI = squeeze(AvgSpikeRates(3,:,:));
    fE1 = squeeze(muaOscFreqs(1,:,:));
    fE2 = squeeze(muaOscFreqs(2,:,:));
    fI = squeeze(muaOscFreqs(3,:,:));
    pE1 = squeeze(muaAreaPowers(1,:,:));
    pE2 = squeeze(muaAreaPowers(2,:,:));
    pI = squeeze(muaAreaPowers(3,:,:));
    
    figure('position',[250 130 1530 840]);
    subplot(2,3,1) % [Cxy]
    imagesc(coh11); xlabel(var2); ylabel(var1); title('coh(input)');
    set(gca,'xtick',1:nval2,'xticklabel',val2)
    set(gca,'ytick',1:nval1,'yticklabel',val1)
    subplot(2,3,2) % [rE]
    imagesc(rE1); xlabel(var2); ylabel(var1); title('rE(input)');
    set(gca,'xtick',1:nval2,'xticklabel',val2)
    set(gca,'ytick',1:nval1,'yticklabel',val1)
    subplot(2,3,3) % [rout]
    imagesc(rE2); xlabel(var2); ylabel(var1); title('rE(output)');
    set(gca,'xtick',1:nval2,'xticklabel',val2)
    set(gca,'ytick',1:nval1,'yticklabel',val1)
    subplot(2,3,4) % (Cxy,rout)
    plot(coh11,rE2,'.-'); xlabel('coh(input)'); ylabel('rE(output)');
    legend(cellfun(@num2str,num2cell(val2),'uni',0),'Location','NorthWest')
    subplot(2,3,5) % (rE,rout)
    plot(rE1,rE2,'.-'); xlabel('rE(input)'); ylabel('rE(output)');
    subplot(2,3,6) % (rE,Cxy)
    plot(rE1,coh11,'.-'); xlabel('rE(input)'); ylabel('coh(input)');
    set(gcf,'name',sprintf('gEI=%g, gIE=%g',gEI,gIE))
    
    % plot fE1, pE1
    
  end
  
%   figure('position',[100 10 900 950]);
%   nc=3; linewidth=2;
%   minasr=min(AvgSpikeRates(:)); maxasr=max(AvgSpikeRates(:));
%   minmof=min(muaOscFreqs(:)); maxmof=max(muaOscFreqs(:));
%   for stim=1:nval1
%     x=val2;
%     mof=squeeze(muaOscFreqs(:,stim,:));
%     asr=squeeze(AvgSpikeRates(:,stim,:));
%     subplot(nval1,nc,nc*(stim-1)+1); plot(x,mof(1,:)./x,'b*-',x,mof(2,:)./x,'ro-','linewidth',linewidth);
%     xlabel(var2); if stim==1, ylabel('MUA freq / f0 [Hz]'); end
%     axis tight; ylim([0 1.25]); hline(.5,'k');
%     title(sprintf('%s=%g',var1,val1(stim)));
%     subplot(nval1,nc,nc*(stim-1)+2); plot(x,asr(1,:)./x,'b*-',x,asr(2,:)./x,'ro-','linewidth',linewidth);
%     xlabel(var2); axis tight; ylim([0 1.25]); hline(.5,'k');
%     if stim==1, ylabel('Avg Spike Rate / f0'); end
%     subplot(nval1,nc,nc*(stim-1)+3); plot(mof(1,:),asr(1,:),'b*-',mof(2,:),asr(2,:),'ro-','linewidth',linewidth);
%     xlabel('MUA freq [Hz]'); if stim==1, ylabel('Avg Spike Rate [Hz]'); end
%     axis tight; 
%     hold on; plot(1:max(mof(:)),1:max(mof(:)),'k--');
%     if stim==1
%       title(sprintf('EsomaType=%g,InputType=%g (f0=%g),tauI=%g',EsomaType(1),InputType,f0,tauI));
%       legend('E','I','diag');
%     else
%       title(sprintf('%s=%g',var1,val1(stim)));    
%     end
%   end  

end

if 0
  %dat=squeeze(data(1).epochs.data(4,:,:));
  dat=squeeze(data(1).epochs.data(3,:,:));
  figure; plot(data(1).epochs.time,dat(:,1));
  ylim([-70 -60]);
end

if 0
  % post-hoc plots
  tauI=val1;
  f0=val2; npks=2;
  pks3=nan(nval1,npks);
  frq3=nan(nval1,npks);
  asr3=nan(nval1,npks);
  for stim=1:nval1
    mof=squeeze(muaOscFreqs(1,stim,:));
    asr=squeeze(AvgSpikeRates(1,stim,:));
    ht=prctile(asr(:),75);
    [pks,locs]=findpeaks(asr,'MinPeakHeight',ht);
    if numel(locs)>npks, locs=locs(1:npks); end
    pks3(stim,1:numel(locs))=f0(locs);
    frq3(stim,1:numel(locs))=mof(locs);
    asr3(stim,1:numel(locs))=asr(locs);
  end
  figure; 
  subplot(3,1,1); plot(tauI,pks3,'*-'); xlabel(var1); ylabel('input f0 [Hz]'); title('resonance frequencies'); axis tight
  subplot(3,1,2); plot(tauI,asr3,'*-'); xlabel(var1); ylabel('avg spike rate [Hz]'); axis tight  
  subplot(3,1,3); plot(tauI,frq3,'*-'); xlabel(var1); ylabel('MUAfreq [Hz]'); axis tight
  legend('1st peak','2nd peak')

  % similarly: sPIN_PFCcell-ci196_tauI15_sin6-40Hz_Estim2-7.mat
  % loop Estim, plot: (Estim,freqs)|tauI=15

  % repeat for sPIN_PFCcell-ci196_tauI10_sin2-60Hz_Estim2-7.mat
  % loop Estim, plot: (Estim,freqs)|tauI=10

  % also plot bandwidth? max/min freqs?
end

if 0
  figure;
  subplot(1,3,1); plot(val1,muaOscFreqs(1,:),'b*-'); xlabel(var1); ylabel('MUA freq [Hz]');
  axis tight; 
  if InputType==1
    title(sprintf('EsomaType=%g,InputType=%g,tauI=%g',EsomaType(1),InputType,tauI));
  else
    title(sprintf('EsomaType=%g,InputType=%g,tauI=%g',EsomaType(1),InputType,tauI));
  end
  subplot(1,3,2); plot(muaOscFreqs(1,:),muaAreaPowers(1,:),'b*-');axis tight; 
  xlabel('MUA freq [Hz]'); ylabel('Rhythm Power Area');
  subplot(1,3,3); plot(val1,AvgSpikeRates(1,:),'b*-'); 
  xlabel(var1); ylabel('Avg Cell Spike Rate [Hz]'); axis tight;
end


