% Purpose: explore electrophysiological properties of various cell models
cd /projectnb/crc-nak/sherfey/projects/ACC_rhythm_processing/model_files

%% SIMULATOR CONTROLS (what to do)
tspan=[0 1000];  % [ms], time limits on numerical integration
dt=.01;         % [ms], time step for numerical integration
solver='rk2';   % {options; 'rk4', 'rk2', 'rk1'}, solver to use
compile_flag=0; % whether to compile the simulation
verbose_flag=1; % whether to display log information
downsample_factor=10;
solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag,'downsample_factor',downsample_factor};
vary=[]; modifications=[];
% -------------------------------------------------------------------------

% Sodium:     (iNaF,iNaP)
% Potassium:  (iKDR,iKS,iM,iKCa,iAHP)
% Calcium:    (iHVA,iCaH,iCaT,CaBuffer)
% Nonspecific: ih,ileak
Cm=1.2; % Durstewitz 2002

taudCa=500; taurCa=80; CaRest = 50/1000; 	% umol/l, resting calcium concentration

% explore interactions among {CaBuffer,iHVA,iKCa}
base='dV/dt=(@current+Iinj*(t>100))/Cm; monitor functions;Cm=1.2; Iinj=10; V(0)=-65;';
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iHVA,iKCa}';
mods={'','gleak',.15;'','Eleak',-75;'','Iinj',20;'','gh',0;...
      '','gHVA',2.5;'','CaRest',50/1000;'','taurCa',[500];'','taudCa',80;...
      '','gKCa',.5;'','akca_scale',[1];'','bkca_scale',[1]};
vary={'','gHVA',2.5;'','gKCa',[.1 .5 .75];'','gNaF',[30 50 70]};
data=SimulateModel(ApplyModifications([base mechanisms],mods),'vary',vary,solver_options{:});
PlotData(data,'variable',{'V','cai','IKCa'});
PlotData(data,'variable','V');
PlotData(data(1),'variable',{'V','IHVA','ICaT','cai','IKCa','IAHP'});

% add iM
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iHVA,iKCa,iM}';
mods={'','gleak',.15;'','Eleak',-75;'','Iinj',20;'','gh',0;'','gM',1;...
      '','gHVA',2.5;'','CaRest',50/1000;'','taurCa',[500];'','taudCa',80;...
      '','gKCa',.5;'','akca_scale',[1];'','bkca_scale',[1]};
vary={'','gM',[0];'','gHVA',[1 2 3];'','gKCa',[0 .25 .5];'','taurCa',[500];'','taudCa',80};
data=SimulateModel(ApplyModifications([base mechanisms],mods),'vary',vary,solver_options{:});
PlotData(data(end-1),'variable',{'V','INaP','IHVA','ICaT','cai','IKCa','IAHP','IM'});
PlotData(data,'variable',{'V','cai','IKCa','IM'});

% explore interactions among {CaBuffer,iHVA,iAHP}
base='dV/dt=(@current+Iinj*(t>100))/Cm; monitor functions;Cm=1.2; Iinj=10; V(0)=-65;';
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iHVA,iAHP}';
mods={'','gleak',.15;'','Eleak',-75;'','Iinj',20;'','gh',0;...
      '','gHVA',2.5;'','CaRest',50/1000;'','taurCa',[500];'','taudCa',80;...
      '','gAHP',.5;'','aAHP_scale',1;'','gKDR',4};
vary={'','gHVA',[.25 .5 1];'','gAHP',[.1 .2 .3];'','gKDR',4.5;'','Iinj',25};
data=SimulateModel(ApplyModifications([base mechanisms],mods),'vary',vary,solver_options{:});
PlotData(data,'variable','V');
PlotData(data,'variable',{'cai','IAHP','V'});
PlotData(data(1),'variable',{'V','IHVA','ICaT','cai','IKCa','IAHP'});
PlotData(data(end),'variable',{'pop1_iAHP_f','pop1_iAHP_minf','pop1_iAHP_mtau','IAHP','pop1_iAHP_m'});

% explore interactions among {CaBuffer,iCaH,iAHP}


amps=-50:5:400; I=1:10:61;
amps=-50:50:400; I=1:10;
model=ApplyModifications([base mechanisms],mods);
data=ProbeCellProperties(model,'amplitudes',amps,'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(I))
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

base='dV/dt=(@current+Iinj*(t>100))/Cm; Cm=1.2; Iinj=10; V(0)=-65;';

% (iCaT,iKCa)
eqns=[base '{iNaF,iKDR,ileak,CaBuffer,iCaT,iKCa}'];
vary={'','gCaT',[0 1 2 5]; '','gKCa',[0 1 2 5]; '','Iinj',[10 20]};
PlotData(SimulateModel(eqns,'vary',vary,solver_options{:}));
% (iCaT,iAHP)
eqns=[base '{iNaF,iKDR,ileak,CaBuffer,iCaT,iAHP}'];
vary={'','gCaT',[0 1 2 5]; '','gAHP',[0 .01 .02 .05]; '','Iinj',[10 20]};
PlotData(SimulateModel(eqns,'vary',vary,solver_options{:}));
% (iHVA,iKCa)
eqns=[base '{iNaF,iKDR,ileak,CaBuffer,iHVA,iKCa}'];
vary={'','gHVA',[0 1 2 5]; '','gKCa',[0 1 2 5]; '','Iinj',[10 20]};
PlotData(SimulateModel(eqns,'vary',vary,solver_options{:}));
% (iHVA,iAHP)
eqns=[base '{iNaF,iKDR,ileak,CaBuffer,iHVA,iAHP}'];
vary={'','gHVA',[0 1 2 5]; '','gAHP',[0 .01 .02 .05]; '','Iinj',[10 20]};
PlotData(SimulateModel(eqns,'vary',vary,solver_options{:}));

% (iCaT,iKCa,ih)
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
% qualitative examination of a model set
eqns=['dV/dt=(@current+Iinj*(t>100))/Cm; Cm=1.2; Iinj=10; V(0)=-65;' mechanisms];
vary={'','gleak',.05; '','gh',[0 1 10 20 30 50]; '','gCaT',2; '','gKCa',1; '','Iinj',[5 7.5 10]};
PlotData(SimulateModel(eqns,'vary',vary,solver_options{:}));
% quantitative characterization of a particular model
eqns=['dV/dt=@current/Cm; Cm=1.2; V(0)=-65;' mechanisms];
mods={'','gleak',.05; '','gh',25; '','gCaT',2; '','gKCa',1};
model=ApplyModifications(eqns,mods);
data=ProbeCellProperties(model,'amplitudes',-50:5:500,'compile_flag',1,'verbose_flag',1,'solver','rk2');
stats=CalcCellProperties(data,'plot_flag',1); 
stats.pop1


% -----------------------------
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,ileak}; Cm=1.2; V(0)=-65';
data=SimulateModel(eqns,'vary',{'','Iinj',0:5:30})
PlotData(data)

eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,iNaP,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,iKS,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,iM,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,ih,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,iCaH,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,iCaT,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,CaBuffer,iHVA,ileak};gHVA=5; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,CaBuffer,iHVA,iKCa,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,CaBuffer,iHVA,iAHP,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,CaBuffer,iHVA,iAHP,iKCa,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iKDR,CaBuffer,iHVA,iCaT,iAHP,iKCa,ileak}; Cm=1.2; V(0)=-65';
eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iNaP,iKDR,iKS,iM,CaBuffer,iHVA,iCaH,iCaT,iAHP,iKCa,ileak,ih}; Cm=1.2; V(0)=-65';
data=SimulateModel(eqns,'vary',{'','Iinj',0:5:30})
PlotData(data)

eqns='dV/dt=(@current+Iinj*(t>50))/Cm; {iNaF,iNaP,iKDR,iKS,iM,ih,CaBuffer,iHVA,iCaH,iCaT,iAHP,iKCa,ileak}; Cm=1.2; V(0)=-65';
tic
data=SimulateModel(eqns,'vary',{'','Iinj',-10:10:80;'','Eleak',[-70 -60 -50]},'tspan',[0 1000])
toc
PlotData(data)

eqns='dV/dt=@current/Cm; {iNaF,iNaP,iKDR,iKS,iM,ih,CaBuffer,iHVA,iCaH,iCaT,iAHP,iKCa,ileak}; Cm=1.2; V(0)=-65';
data=ProbeCellProperties(eqns,'amplitudes',-50:50:400,'compile_flag',1,'verbose_flag',1);
stats=CalcCellProperties(data); 
stats.pop1

data=SimulateModel('dV/dt=@current-.4*(V+70)+I; {RSiNaF,RSiKDR}; V(0)=-65;','vary',{'','I',0:.2:1});
PlotData(data)
data=SimulateModel('dV/dt=@current-.4*(V+80)+I; {RSiNaF,RSiKDR,iCaT}; V(0)=-85;','vary',{'','I',0:5});
PlotData(data)
data=SimulateModel('dV/dt=@current-.4*(V+80)+I; {RSiNaF,RSiKDR,iCaT,iKCa,CaBuffer}; V(0)=-85;','vary',{'','I',[4 6];'','gKCa',0:5:20},'tspan',[0 500]);
PlotData(data)

% long AHP
data=SimulateModel('dV/dt=@current-.4*(V+70)+I; {RSiNaF,RSiKDR,iHVA,iAHP,CaBuffer}; V(0)=-65;','vary',{'','I',1;'','gHVA',[0 .025 .05 .1];'','gAHP',[0 10 20]},'tspan',[0 1000]);
data=SimulateModel('dV/dt=@current-.4*(V+70)+I; {RSiNaF,RSiKDR,iCaH,iKCa,CaBuffer}; V(0)=-65;','vary',{'','I',[1 2 10];'','gKCa',1},'tspan',[0 500]);
data=SimulateModel('dV/dt=@current-.4*(V+70)+I; {RSiNaF,RSiKDR,iCaT,iKCa,CaBuffer}; V(0)=-65;','vary',{'','I',0;'','gKCa',0:5:20},'tspan',[0 500]);
PlotData(data)

data=SimulateModel('dV/dt=@current-.4*(V+70)+I; {RSiNaF,RSiKDR,iCaH,iKCa,CaBuffer}; V(0)=-65;','vary',{'','I',[1 2 10];'','gKCa',1},'tspan',[0 500]);


mechanisms='{iNaF,iNaP,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
eqns=['dV/dt=(@current+Iinj)/Cm; Cm=1.2; Iinj=5; V(0)=-65;' mechanisms];
mods={'','gleak',.15;'','Eleak',-75;'','noise',5; 
  '','gh',.2; '','gCaT',1.25; '','gKCa',.5;'','gKDR',4.25;'','gNaF',55;'','gNaP',1};
model=ApplyModifications(eqns,mods);
data=SimulateModel(eqns,'vary',{'','Iinj',[20];'','gNaP',[0 1 2 5 10 20 50 100]})
PlotData(data)
data=ProbeCellProperties(model,'amplitudes',-50:5:400,'compile_flag',1,'verbose_flag',1,'solver','rk2');
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1
PlotData(data(50:60))

% -------------------------------------------------------------------------

% IPs clustered:
% {'AHP(time2trough)','RMP','Ih','median(ISIs)','ISI2/ISI3'}.

% Classes found:
% Class 1 has long AHP, slow adaptation, and dominates superficial layers.
% Class 2 has short AHP, rapid adaptation, and dominates deep layers.
% Class 3 has longest AHP, largest Ih sag, 
% Class 4 is rare and evenly distributed, with fastest/no AHP, elevated RMP
%         and abrupt stopping (few spikes, with ADP)
% Class 5 may have bigger Ih and slightly faster adaptation than Class 1, 
%         but otherwise similar IPs; is more concentrated in middle layers.

% Potential relation to Fiona's manual Groups:
% Group 1 <-> Class 3
% Group 2 <-> Class 1 (and maybe Class 5 in middle layers)
% Group 3 <-> Class 2
% Group 4 <-> Class 4

% --------------------------------------
% MODEL CELLS:
%   Base on classes 1, 2, and 4
%   Using IPs:
%     from clustering: {'AHP(time2trough)','RMP','Ih','median(ISIs)','ISI2/ISI3'}
%     significant diffs (p<.05): 
%       (1,2): {'ISI2/ISI3','median(ISIs)','AHP(time2trough)'}
%       (1,4): {'RMP','Ih','median(ISIs)','AHP(time2trough)'}
%       (2,4): {'RMP','Ih'}
% MODEL LAYERS:
%   Superficial: Class 1 dominates, with a little class 4
%   Deep layer: Classes 1, 2, and 4 contribute equally
% --------------------------------------

%{
Intrinsic Properties by Cell Class (mean/std,N)
                AHP(time2trough)	SpikeAmp        SpikeWidth        ThreshRate      RMP             1/ISI1          median(ISIs)		Ih                Ih-decay        ISI2/ISI3		
Class1 (n=21): 	(31  /12 ,21)			(66  /12 ,21)		(1.5 /0.28,21)		(1.8 /1.3,21)		(-75 /10 ,21)		(32  /13 ,21)		(11  /3.8,19)		(0.88/0.66,14)		(85  /63 ,14)   (0.9 /0.21,21)	
Class2 (n=13): 	(21  /12 ,13)			(68  /12 ,13)		(1.5 /0.32,13)		(1.7 /2  ,11)		(-75 /6.4,13)		(27  /15 ,13)		(22  /8.8,12)		(0.95/0.44,12)		(110 /67 ,12)		(0.56/0.13, 9)	
Class3 (n=7) : 	(51  /20 , 7)			(66  /8.4, 7)		(1.4 /0.17, 7)		(3.6 /1.5, 7)		(-71 /5.9, 7)		(37  /24 , 7)		(12  /7.3, 5)		(3.4 /1.1 , 7)		(130 /40 , 7)		(0.92/0.2 , 6)	
Class4 (n=6) : 	(9.8 /7.6, 3)			(59  /4.4, 3)		(1.5 /0.42, 3)		(2.6 /3.3, 3)		(-64 /3.3, 6)		(26  /17 , 6)		(31  /9.1, 6)		(2.3 /0.59, 6)		(92  /53 , 6)		(NaN /NaN , 0)	
Class5 (n=14): 	(39  /20 ,14)			(69  /4.4,14)		(1.4 /0.24,14)		(1.4 /1  ,13)		(NaN /NaN, 0)		(33  /29 ,14)		(14  /10 ,14)		(1.3 /0.95, 8)		(110 /48 , 8)		(0.73/0.18,11)	

Significant differences between cell classes:
#12-(Class1,Class2):
          ISI2/ISI3: p=2.4e-05
       median(ISIs): p=0.0019
   AHP(time2trough): p=0.034
#14-(Class1,Class4):
                RMP: p=0.00027
                 Ih: p=0.00079
       median(ISIs): p=0.003
   AHP(time2trough): p=0.018
#24-(Class2,Class4):
                RMP: p=9.5e-05
                 Ih: p=0.0013
   
%}



%% ACC L2/3/5 PY (rat ACd: LeBeau lab)
% Target properties (n=20, LeBeau lab):
% Homogeneous (std/mean < .4):
%  *IP5? AHP(5-20ms) [?]: 9.25e4 ? 2.79e4 (note: if not IP5, probably similar)
%  *IP6: SpikeWidth [ms]: 1.62 ? .279
%   IP8: RMP [mV]:      -74.0 ? 11.1
%   n/a: SpikeAmp [mV]:  65.8 ? 10.3
% Heterogeneous (std/mean > .4):
%   IP7: ThreshRate [Hz]: 2.41 ? 2.27 (L2/3a), 2.04 ? 1.18 (L3b/5)
%   IP9: StepRate [Hz]:  26.8 ? 12.0 (L2/3a), 29.3 ? 24.9 (L3b/5)
%   n/a: Ih:        [mV]: 1.79 ? 1.18
% Note: * means significantly different across layers (ttest2, p<.05).

% For IPs across all cells in all layers:
% cumulative explained variance: [23% 43% 63% 82% 88% 94% 100%]
% intrinsic property number:     [IP8 IP5 IP6 IP7 IP9 IP3 IP10]
% percent explained variance:    [23% 20% 20% 19% 6%  6%  6%  ]
% Cluster analysis identified 5 classes based on (IP8,IP5,IP6,IP7) (across all layers).
% IP5 and IP6 differed between layers; IP8 and IP9 were highly heterogeneous.
%{
  Class 1 cells (33% of recorded cells) had relatively narrow
  spikes (IP6) and higher spike rates in response to a depolarizing step (IP9). Class 2 
  (33%) was similar to Class 1 except for having a more depolarized resting potential 
  (IP8) and slightly shorter time to maximum deflection of the AHP (IP5). Class 3 (15%) 
  had a much longer-lasting AHP (IP5). Classes 4 (10%) and 5 (9%) both had high 
  resting potentials (IP8) and were differentiated by spike width (IP6) and threshold 
  firing rate (IP7).
  Analysis of the overlap between the different clustering approaches
  (manual versus iterative k-means) showed some congruence between neurons 
  included in Group 2 and Class 1 (38%) but the overlap was less for other pairs.

Hypothesis:
Class 1 = ACC PY L5b/L6 because (high StepRate, narrow spikes, low RMP)
Class 2 = ACC PY L5b/L6 because (similar to 1 except elevated RMP)
Class 3 = ACC PY L2/3   because (much longer-lasting AHP, low RMP)
Class 4 = ACC PY L6     because (narrow spikes, fast rate; elevated RMP, short AHP)
Class 5 = ACC PY L5     because (broad spikes, slow rate; elevated RMP, short AHP)

TODO: compare layers for cells in each class identified by cluster
analysis. calculate mean/std IPs for each class; use these values to
establish my cell model constraints with layer-specificity.
NEXT: CALCULATE MEAN/STD IPS FOR EACH CLASS

Intrinsic Properties by Cell Class (mean,std,N)
          AHP(time2trough)	SpikeAmp        SpikeWidth      ThreshRate        RMP             1/ISI1          median(ISIs)		Ih                Ih-decay          ISI2/ISI3		
Class1: 	(28  /8.9,20)			(72  /10 ,20)		(1.4 /0.2 ,20)	(1.3 /0.79,18)		(-81 /5.3,14)		(36  /25 ,20)		(18  /11 ,18)		(0.68/0.28,12)		(96  /71 ,12)     (1.6 /3.7  ,18)	
Class2: 	(23  /11 ,19)			(64  /7.7,19)		(1.4 /0.18,19)	(1.6 /0.82,19)		(-69 /4.8,19)		(28  /17 ,21)		(16  /8.9,21)		(1.5 /1.1 ,19)    (101 /54 ,19)     (0.83/0.26 ,18)	
Class3: 	(62  /11 , 9)			(69  /6.9, 9)		(1.5 /0.24, 9)	(1.8 /1.1 , 9)    (-80 /5.3, 4)		(26  /17 , 9)		(8.8 /2.9, 8)		(1.6 /0.81, 7)		(126 /53 , 7)     (0.8 /0.085, 7)	
Class4: 	(34  /23 , 6)			(59  /9.9, 6)		(1.5 /0.3 , 6)	(5.9 /0.98, 6)		(-70 /9.8, 6)		(40  /17 , 6)		(18  /10 , 4)		(4.1 /0.63, 4)		(87  /36 , 4)     (0.96/0.22 , 4)	
Class5: 	(24  /17 , 4)			(59  /6.3, 4)		(2.1 /0.19, 4)	(0.84/0.98, 3)		(-64 /6.2, 4)		(21  /9.3, 5)		(26  /8.7, 5)		(1.6 /0.44, 5)		(113 /65 , 5)     (0.77/0    , 1)	

          AHP(0-5ms)              AHP(5-20ms)             AHP(20-200ms)         AHP(total)              AHP(trough)		
Class1: 	(2.2e+04/5.2e+03,20)		(9.9e+04/2.2e+04,20)		(9.6e+05/3e+05,20)		(1.1e+06/3.3e+05,20)		(8.4 /2.3,20)	
Class2: 	(2.7e+04/5.9e+03,19)		(1.1e+05/2.1e+04,19)		(9.4e+05/2.1e+05,19)  (1.1e+06/2.3e+05,19)		(8.8 /2.4,19)	
Class3: 	(1.5e+04/5.7e+03, 9)		(7.3e+04/2.1e+04, 9)		(1.1e+06/1.4e+05, 9)	(1.2e+06/1.6e+05, 9)		(8.9 /1.7, 9)	
Class4: 	(1.7e+04/9.6e+03, 6)		(6.7e+04/2.9e+04, 6)		(7e+05/3.1e+05, 6)		(7.9e+05/3.5e+05, 6)		(4.9 /2.7, 6)	
Class5: 	(2.3e+04/4.7e+03, 4)		(9.4e+04/3.2e+04, 4)		(9.3e+05/4e+05, 4)		(1e+06/4.4e+05, 4)    	(7.4 /2.2, 4)	

%}

% Only network beta rhythms (15-30Hz) were present in L2/3/5 (in vitro) with kainate.

% Cell models:
% ...



%{
	- ion channel functions: (need reference for each mechanism model used)
		(CaT+AHP): spike-frequency adaptation and AHP
		(KA): spike-frequency adaptation; may affect presence of rhythm (Tallie statement)
		(M): spike-frequency adaptation; may control separation of beta/gamma (Fiona result)
		(h): controls gamma frequency and beta/gamma separation (Fiona result)
		(NaP): keeps the gamma going (Tallie result)
		(h+NaP): STOs (Horacio: Remme et al)
		(CaT): STOs in superior olive (Horacio: Manor et al)
	- note this in comments: (https://docs.google.com/spreadsheets/d/1emhIL3LNchq5Pz7GAd64u9tqsf9ysiZa1rt3nTENsCw/edit#gid=0)
		Eleak: IP5/8 up, IP7 down
		gcat: IP5/6/8/9 down, IP7 up; 
		gh: IP7/9 up, IP5/6/8 mixed effects; 
		gnaf: IP5/6 up, IP7/8/9 mixed effects; 
		where IP5=AHPtau [ms], IP6=spikewidth [ms], IP7=threshrate [Hz], IP8=RMP [mV], IP9=steprate [Hz]
	* idea: try to create cells with/without adaptation, AHP, and STOs, and with superficial vs deep IPs by varying (gCaT,gAHP,gNaP) in (Na,K,CaT,AHP,NaP,h)
		- note: choose mechanism models that can be used in the PFC model; record references in driver script
	- determine how many cells and network realizations are needed to get good stats reflecting the heterogeneity (note this in comments)
	* also: per Fiona experimental results, want gamma to involve and depend on cells with h- and NaP-currents, and to be present when there are INfast but not INslow without INfast
		- use these facts to tune baseline conductance levels for at least h and NaP
	- include h-current b/c there is a clear sag observed in Type A cells which have IPSPs that are rhythmic with the gamma rhythm; and NaP b/c necessary for rhythm persistence
%}

% -------------------------------------------------------------------------
%% ACC L6 PY (rat ACd: LeBeau lab)
% Target properties (n=20, LeBeau lab):
% Homogeneous (std/mean < .4):
%  *IP5? AHP(5-20ms) [?]: 1.04e5 ? 1.77e4 (smaller than L2/3/5)
%  *IP6: SpikeWidth [ms]: 1.38 ? .222     (smaller than L2/3/5)
%   IP8: RMP [mV]:      -74.4 ? 8.32
%   n/a: SpikeAmp [mV]:  64.7 ? 9.38
% Heterogeneous (std/mean > .4):
%   IP7: ThreshRate [Hz]: 1.52 ? 1.44
%   IP9: StepRate [Hz]:  30.1 ? 14.6      (maybe faster than L2/3/5)
%   n/a: Ih:        [mV]: 1.39 ? 1        (maybe smaller than L2/3/5)
% Note: * means significantly different across layers (ttest2, p<.05).

% Network beta (15-30Hz) and gamma (30-45Hz) rhythms were present in L5/6 (in vitro) with kainate.
% Speculation: L6 PY IPs may support faster activity more than the PY IPs in
% more superficial layers. Also, greater numbers of PV+ FS interneurons with
% faster inhibition (smaller time constants) could support the additional
% gamma seen in deep layers.

% Cell models:
% ...

% -------------------------------------------------------------------------
%% PFC L2/3 PY (monkey DLPFC)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...
% -------------------------------------------------------------------------
%% PFC L2/3 PY (rat PrL)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...
% -------------------------------------------------------------------------
%% PFC L5/6 PY (monkey DLPFC)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...
% -------------------------------------------------------------------------
%% PFC L5/6 PY (rat PrL)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...
% -------------------------------------------------------------------------
%% PFC/ACC FS IN (PV+)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...
% -------------------------------------------------------------------------
%% PFC/ACC LTS IN (CB+)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...
% -------------------------------------------------------------------------
%% PFC/ACC RSNP IN (CB+)
% Target properties:
% RMP: ... [Ref#]
% ...

% Cell models:
% ...

% -------------------------------------------------------------------------
% Prepare for Network Simulations
L='NaKMCa';%'PFC_L23';
nE=1;

% Tonic injected current parameters
Iapp=0;   % uA/cm2, amplitude of injected current
% Poisson noise parameters
baseline_rate=.1; % kHz, baseline poisson rate
gNOISE=.03;        % mS/cm2, noise scaling factor 
kick=1;            % mS/cm2, generic conductance kick per spike

input_parameters={'baseline_rate',baseline_rate,'kick',kick,'gNOISE',gNOISE,'Iapp',Iapp};

% Homogeneous biophysical parameters
Cm=1;         % uF/cm2, membrane capacitance
v_IC=-65;     % mV, voltage initial condition
v_IC_noise=1; % mV, scale of normally distributed IC noise
              % IC: V(0)=v_IC+v_IC_noise*randn
biophysical_parameters={'Cm',Cm,'v_IC',v_IC,'v_IC_noise',v_IC_noise};

% (gNa,gK,gNaP,gleak)  (note: Rleak=1/gleak)
% (ENa,EK,ECa,Eh)

biophysical_parameters=cat(2,{'gM',0,'gleak',.04,'Eleak',-65},biophysical_parameters);
% vary={L,'gM',9;L,'gCa',10:2:30;L,'gCan',0;L,'Iapp',10; L,'baseline_rate',0};
% vary={L,'gM',[5 7.5 10];L,'gCa',0;L,'Eleak',-80:5:-50;L,'Iapp',10; L,'baseline_rate',0};
vary={L,'gM',[0 10];L,'gCa',0;L,'Eleak',-80:5:-50;L,'Iapp',0; L,'baseline_rate',1;L,'gNOISE',[.01 .02]};

% Specify cell model
s=[];
s.pops=ACC_Ecell_specification(L);
s.pops.name=L;
s.pops.size=nE;
s.pops.parameters={input_parameters{:},biophysical_parameters{:}};

% Simulate
data=SimulateModel(s,'vary',vary,'modifications',modifications,solver_options{:});
PlotData(data,'plot_type','waveform')

%%

PlotData(data,'plot_type','power','xlim',[0 100])
PlotData(data,'plot_type','rastergram')


