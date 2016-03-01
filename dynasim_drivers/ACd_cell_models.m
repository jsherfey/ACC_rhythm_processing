% IPs clustered:
% {'AHP(time2trough)','RMP','Ih','median(ISIs)','ISI2/ISI3'}.
% See ~/projects/ACC_rhythm_processing/ACd_cell_classes.m for complete clustering results
% See ~/projects/ACC_rhythm_processing/ACC_clustering_all_IPs.m for cluster analysis
% See ~/projects/ACC_rhythm_processing/ACC_classes_by_layer.m for class comparison

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
%   Base on classes 1, 2, and 4.
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

Significant differences between cell classes 1, 2, and 4: (ttest2,'vartype','unequal')
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

rootoutdir='/projectnb/crc-nak/sherfey/projects/ACC_rhythm_processing/dynasim_studies';
model_path='/projectnb/crc-nak/sherfey/projects/ACC_rhythm_processing/model_files';
addpath(model_path);

%% Test cell models before fitting models to particular classes
cd(model_path);

% -------------------------------------------------------------------------
% Step 1: qualitative search (look for approx regimes)
% Method: use local SimulateModel with 'vary' and PlotData
% Goal: record mechanisms + parameter sets that loosely resemble desired classes 
% Record observations in 'notes' file in 'dynasim_studies/cells/test_sweeps'
% -------------------------------------------------------------------------

% Qualitative descriptions:
% Class 1 has long AHP, slow adaptation, and dominates superficial layers.
% Class 2 has short AHP, rapid adaptation, and dominates deep layers.
% Class 3 has longest AHP, largest Ih sag, 
% Class 4 is rare and evenly distributed, with fastest/no AHP, elevated RMP
%         and abrupt stopping (few spikes, with ADP)
% Class 5 may have bigger Ih and slightly faster adaptation than Class 1, 
%         but otherwise similar IPs; is more concentrated in middle layers.

tspan=[0 300];  % [ms], time limits on numerical integration
dt=.01;         % [ms], time step for numerical integration
solver='rk2';   % {options; 'rk4', 'rk2', 'rk1'}, solver to use
compile_flag=1; % whether to compile the simulation
verbose_flag=1; % whether to display log information
downsample_factor=10;
solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag,'downsample_factor',downsample_factor};
save_data_flag=1; save_plots_flag=1; plot_options={'visible','off'};

% Available mechanisms:
% Sodium:     (iNaF,iNaP, RSiNaF)
% Potassium:  (iKDR,iKS,iM,iKCa,iAHP, RSiKDR)
% Calcium:    (iHVA,iCaH,iCaT,CaBuffer)
% Nonspecific: (ih,ileak)

% long AHP
mechanisms='{RSiNaF,RSiKDR,iHVA,iAHP,CaBuffer,ileak}';
vary={'gHVA',[0 .025 .05 .1];'gAHP',[0 10 20]};
mechanisms='{RSiNaF,RSiKDR,iCaH,iKCa,CaBuffer,ileak}';
vary={'gKCa',1;'Iinj',[1 2 10]};
mechanisms='{RSiNaF,RSiKDR,iCaT,iKCa,CaBuffer,ileak}';
vary={'gKCa',0:5:20;'Iinj',[0]};
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
mods={'gCaT',2;'gKCa',.2;'gleak',.4;'Eleak',-75;'noise',0};
vary={'gh',[0 1 10 20 30 50];'Iinj',[5 7.5 10]};

% Class 2?:
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
mods={'gCaT',.1;'gKCa',2;'gleak',.4;'Eleak',-75;'noise',0};
vary={'gh',[1];'Iinj',[20];'akca_scale',[10];'bkca_scale',.1};
model=ApplyModifications(eqns,mods);
data=SimulateModel(model,'vary',vary,solver_options{:});
PlotData(data,'ylim',[-100 50]);


mechanisms='{iNaF,iKDR,ileak,iM}'; mods={'gNaF',50;'gKDR',4;'gleak',.4;'Eleak',-75;'Cm',1.2;'Iinj',35}; vary={'gM',0:.1:.5};
mechanisms='{iNaF,iKDR,ileak,iCaT}'; mods={'gNaF',50;'gKDR',4;'gleak',.4;'Eleak',-75;'Cm',1.2;'Iinj',35}; vary={'gCaT',0:.5:3};
mechanisms='{iNaF,iKDR,ileak,iHVA,CaBuffer}'; mods={'gNaF',50;'gKDR',4;'gleak',.4;'Eleak',-75;'Cm',1.2}; vary={'iHVA',0:.5:3;'Iinj',[20 30 40]};
mechanisms='{iNaF,iKDR,ileak,iCaH}'; mods={'gNaF',50;'gKDR',4;'gleak',.4;'Eleak',-75;'Cm',1.2}; vary={'gCaH',0:.5:3;'Iinj',[20 30 40]};

mechanisms='{RSiNaF,RSiKDR,ileak,iM}'; mods={'gNaF',200;'gKDR',20;'gleak',.4;'Eleak',-75;'Cm',1.2}; vary={'Iinj',-10:5:30;'gM',1};
mechanisms='{RSiNaF,RSiKDR,ileak,iCaT}'; mods={'gNaF',200;'gKDR',20;'gleak',.4;'Eleak',-75;'Cm',1.2}; vary={'Iinj',-10:10:30;'gCaT',[0:.1:.5]};
%mechanisms='{RSiNaF,RSiKDR,ileak,iHVA,CaBuffer}'; mods={'gNaF',200;'gKDR',20;'gleak',.4;'Eleak',-75;'Cm',1.2}; vary={'iHVA',0:.5:3;'Iinj',[20 30 40]};

% generic execution
eqns=['dV/dt=(@current+Iinj*(t>100)+noise*rand(1,Npop))/Cm; V(0)=-75; Cm=1.2; Iinj=10; noise=0; monitor functions;' mechanisms];
model=ApplyModifications(eqns,mods);
prefix=['cell_' get_strID(mechanisms,'mech')];
study_dir=fullfile(rootoutdir,'cells','test_sweeps',get_strID(mechanisms,'mech'),[get_strID(mods,'mods') '__' get_strID(vary,'vary')]);
if save_plots_flag==1, plot_functions=@PlotData; else plot_functions=[]; end
data=SimulateModel(model,'vary',vary,solver_options{:},'plot_functions',plot_functions,'plot_options',plot_options,'study_dir',study_dir,'prefix',prefix,'save_results_flag',save_plots_flag,'save_data_flag',save_data_flag);
PlotData(data,'ylim',[-100 50]); file=fullfile(study_dir,'waveforms.jpg'); set(gcf,'PaperPositionMode','auto'); print(gcf,file,'-djpeg');

if 0 % quantitatively assess a particular model
  index=8; m=data(index).model; amps=-50:25:550;
  d=ProbeCellProperties(ApplyModifications(m,{'Iinj',0}),'amplitudes',amps,solver_options{:});%'compile_flag',1,'verbose_flag',1,'solver','rk2');
  PlotData(d(unique(round(linspace(1,length(d),10)))))
  stats=CalcCellProperties(d,'plot_flag',0);
  stats.pop1
  file=fullfile(study_dir,sprintf('sim%g_ProbeCell_AHPtime%g_RMP%g_Ih%g_ISImedian%g_AR%g_APamp%g_APdur%g_ISI%g_FRmin%g.jpg',index,stats.pop1.AHP_time2trough,stats.pop1.RMP,stats.pop1.Ih_abssag,stats.pop1.ISI_median,stats.pop1.AR23,stats.pop1.AP_amp,stats.pop1.AP_dur,stats.pop1.ISI1,stats.pop1.FR_min));
  print(gcf,file,'-djpeg');
end  
% -----

tspan=[0 300];  % [ms], time limits on numerical integration
dt=.01;         % [ms], time step for numerical integration
solver='rk2';   % {options; 'rk4', 'rk2', 'rk1'}, solver to use
compile_flag=1; % whether to compile the simulation
verbose_flag=1; % whether to display log information
downsample_factor=10;
solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag,'downsample_factor',downsample_factor};
save_data_flag=1; save_plots_flag=0; plot_options={'visible','off'};

mods1={'gNaF',50;'gKDR',4;'gleak',.4;'Eleak',-75;'Cm',1.2};
mods2={'gNaF',200;'gKDR',20;'gleak',.4;'Eleak',-75;'Cm',1.2};
gvals=[0 .05 .1 .2 .5 1 2 5 10 20 50]; gvals2=[0 .05 .1 .5 1 5 10];
Ivals=[-5 0 5 10 30 50]; RSIvals=[-2 0 2 3 4]; I0=35; RSI0=3;
common_dir=fullfile(rootoutdir,'cells','test_sweeps');
sets={
  '{iNaF,iKDR,ileak,iKS}', mods1, {'gKS',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,ih}', mods1, {'gh',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,iM}', mods1, {'gM',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,iNaP}', mods1, {'gNaP',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,CaBuffer,iHVA}', mods1, {'gHVA',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,CaBuffer,iCaH}', mods1, {'gCaH',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,CaBuffer,iCaT}', mods1, {'gCaT',gvals;'Iinj',Ivals};
  '{iNaF,iKDR,ileak,CaBuffer,iHVA,iAHP}', cat(1,mods1,{'Iinj',I0}), {'gHVA',gvals;'gAHP',gvals2};
  '{iNaF,iKDR,ileak,CaBuffer,iCaH,iKCa}', cat(1,mods1,{'Iinj',I0}), {'gCaH',gvals;'gKCa',gvals2};
  '{iNaF,iKDR,ileak,CaBuffer,iCaH,iAHP}', cat(1,mods1,{'Iinj',I0}), {'gCaH',gvals;'gAHP',gvals2};
  '{iNaF,iKDR,ileak,CaBuffer,iCaT,iKCa}', cat(1,mods1,{'Iinj',I0}), {'gCaT',gvals;'gKCa',gvals2};
  '{iNaF,iKDR,ileak,CaBuffer,iCaT,iAHP}', cat(1,mods1,{'Iinj',I0}), {'gCaT',gvals;'gAHP',gvals2};
  '{RSiNaF,RSiKDR,ileak,iKS}', mods2, {'gKS',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,ih}', mods2, {'gh',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,iM}', mods2, {'gM',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,iNaP}', mods2, {'gNaP',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iHVA}', mods2, {'gHVA',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iCaH}', mods2, {'gCaH',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iCaT}', mods2, {'gCaT',gvals;'Iinj',RSIvals};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iHVA,iAHP}', cat(1,mods2,{'Iinj',RSI0}), {'gHVA',gvals;'gAHP',gvals2};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iCaH,iKCa}', cat(1,mods2,{'Iinj',RSI0}), {'gCaH',gvals;'gKCa',gvals2};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iCaH,iAHP}', cat(1,mods2,{'Iinj',RSI0}), {'gCaH',gvals;'gAHP',gvals2};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iCaT,iKCa}', cat(1,mods2,{'Iinj',RSI0}), {'gCaT',gvals;'gKCa',gvals2};
  '{RSiNaF,RSiKDR,ileak,CaBuffer,iCaT,iAHP}', cat(1,mods2,{'Iinj',RSI0}), {'gCaT',gvals;'gAHP',gvals2};
};

for s=1:size(sets,1)
  fprintf('SET %g of %g\n',s,size(sets,1));
  mechanisms=sets{s,1};
  mods=sets{s,2};
  vary=sets{s,3};
  eqns=sprintf('dV/dt=(@current+Iinj*(t>100&t<%g)+noise*rand(1,Npop))/Cm; V(0)=-75; Cm=1.2; Iinj=10; noise=0; monitor functions;%s',tspan(2)-50,mechanisms);
  model=ApplyModifications(eqns,mods);
  prefix=['cell_' get_strID(mechanisms,'mech')];
  study_dir=fullfile(rootoutdir,'cells','test_sweeps',get_strID(mechanisms,'mech'),[get_strID(mods,'mods') '__' get_strID(vary,'vary')]);
  if save_plots_flag==1, plot_functions=@PlotData; else plot_functions=[]; end
  data=SimulateModel(model,'vary',vary,solver_options{:},'plot_functions',plot_functions,'plot_options',plot_options,'study_dir',study_dir,'prefix',prefix,'save_results_flag',save_plots_flag,'save_data_flag',save_data_flag);
  PlotData(data,'ylim',[-100 50]); set(gcf,'PaperPositionMode','auto');
  file=fullfile(study_dir,'waveforms.jpg'); print(gcf,file,'-djpeg');
  file=fullfile(common_dir,sprintf('%s__%s__%s__waveforms.jpg',get_strID(mechanisms,'mech'),get_strID(mods,'mods'),get_strID(vary,'vary'))); print(gcf,file,'-djpeg');
end

% -------------------------------------------------------------------------
% Step 2: quantitative assessment of single IP = f(param)
% Method: use experiments and/or local/cluster SimulateModel with experiment 
%         and analysis, varying one parameter at a time; plot (param,IPs)
% -------------------------------------------------------------------------
if 0
  
  % Quantitative properties:
  % TARGETS: [AHP_time2trough]   [RMP]    [Ih_abssag] [ISI_median]    [AR23]
  % Class 1      19-43ms       -85 to -65  .2-1.5mV     7-15ms        .7-1.1
  % Class 2      9-33ms        -81 to -69  .5-1.4mV     13-31ms       .4-.7
  % Class 3      31-71ms       -77 to -65  2.3-4.5mV    5-19ms        .7-1.1
  % Class 4      2-17ms        -67 to -61  1.7-2.9mV    22-40ms        n/a
  % Class 5      19-59ms           ??      .3-2.3mV     4-24ms        .5-.9
  % Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
  %          59-68mV   1.49-1.53ms   31-39ms     1-4Hz

  % quantitative characterization of a particular model
  eqns=['dV/dt=(@current+noise*randn)/Cm; Cm=1.2; noise=5; V(0)=-65;' mechanisms];
  %eqns=['dV/dt=@current/Cm; Cm=1.2; V(0)=-65;' mechanisms];
  mods={'','gleak',.15;'','Eleak',-75;'','noise',5; 
    '','gh',.2; '','gCaT',1.25; '','gKCa',.5;'','gKDR',4.25;'','gNaF',55};

  mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iNaP,iKCa}';
  eqns=['dV/dt=(@current+noise*randn)/Cm; Cm=1.2; noise=5; V(0)=-65;' mechanisms];
  mods={'','gleak',.15;'','Eleak',-75;'','noise',5; 
    '','gh',.15; '','gCaT',1.5; '','gKCa',.5;'','gKDR',4.25;'','gNaF',55;'','gNaP',1};

  mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iNaP,iKCa,iHVA}';
  mods={'','gleak',.15;'','Eleak',-75;'','noise',0;'','gHVA',5;'','vtauHVA',500;
  '','gh',.15; '','gCaT',0; '','gKCa',1;'','gKDR',4.25;'','gNaF',55;'','gNaP',0};

  mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iAHP}';
  mods={'','taurCa',80/7;'','gleak',.15; '','gh',.2; '','gCaT',1.5; '','gAHP',.1;'','Eleak',-75;'','noise',5};

  mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iHVA,iKCa}';
  mods={'','gleak',.15; '','gh',.2; '','gHVA',2.5; '','gKCa',.5;'','Eleak',-75;'','noise',0};

  eqns=['dV/dt=(@current+noise*randn)/Cm; Cm=1.2; noise=5; V(0)=-65;' mechanisms];
  amps=-50:5:400; I=1:10:61;
  amps=-50:50:400; I=1:10;
  model=ApplyModifications(eqns,mods);
  data=ProbeCellProperties(model,'amplitudes',amps,'compile_flag',1,'verbose_flag',1,'solver','rk2');
  PlotData(data(I))
  stats=CalcCellProperties(data,'plot_flag',0); 
  stats.pop1

  % Class 1 model (without match to common targets):
  mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
  eqns=['dV/dt=(@current+noise*randn)/Cm; Cm=1.2; noise=5; V(0)=-65;' mechanisms];
  mods={'','gleak',.15; '','gh',.2; '','gCaT',1.5; '','gKCa',.2;'','Eleak',-75;'','noise',5};
  model=ApplyModifications(eqns,mods);
  data=ProbeCellProperties(model,'amplitudes',-50:10:400,'compile_flag',0,'verbose_flag',1,'solver','rk2');
  PlotData(data(unique(round(linspace(1,length(data),10)))))
  stats=CalcCellProperties(data,'plot_flag',0); 
  stats.pop1
end



%% Class 1: ACC L2/3/5 PY (rat ACd: LeBeau lab) (~ Group 2)
%{
  Distribution across layers (note: similar to Group 2):
   superficial  middle    deep
      43%        33%      24%
%}
%{
Intrinsic Properties (mean +/- std) sorted by Heterogeneity (=std/mean)
Class1 (n=21):
HOMOGENEOUS PROPERTIES (std/mu<.4):
(0.13):              RMP: -74.9 +/- 9.96            (n=21)
(0.18):         SpikeAmp: 65.6 +/- 11.7            (n=21)
(0.19):       SpikeWidth: 1.49 +/- 0.276           (n=21)
(0.23):        ISI2/ISI3: 0.899 +/- 0.209           (n=21)
(0.33):     median(ISIs): 11.5 +/- 3.78            (n=19)
(0.39): AHP(time2trough): 30.7 +/- 12              (n=21)
HETEROGENEOUS PROPERTIES:
(0.41):           1/ISI1: 31.9 +/- 12.9            (n=21)
(0.72):       ThreshRate: 1.82 +/- 1.31            (n=21)
(0.74):         Ih-decay:   85 +/- 62.9            (n=14)
(0.75):               Ih: 0.885 +/- 0.66            (n=14)
%}
%{
Class 1 (wrt Class 2): has 
  slower adaptation (.9>.5)  <- difference consistent w/ sup vs deep PrL (but not IL) 
                                (http://www.neuroelectro.org/article/40619/)
  shorter ISIs (11<22) 
  longer AHPs (31>21)

Class 1 (wrt Class 4): has 
  lower RMP (-75<-64)        <- difference consistent w/ sup vs deep Ctx PY
                                (see L2/3: http://www.neuroelectro.org/neuron/110 and L5/6: http://www.neuroelectro.org/neuron/111)
  less Ih sag (.9<2.3)
  shorter ISIs (11<31)
  longer AHPs (31>10)
%}

%            AHP duration              hyperpol-sag   <ISI>      adaptation
% TARGETS: [AHP_time2trough]   [RMP]   [Ih_abssag] [ISI_median]    [AR23]
%          'AHP(time2trough)'  'RMP'      'Ih'    'median(ISIs)' 'ISI2/ISI3'
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16          .7-1.1
% Class 2      11-32ms       -81 to -69  .5-1.4mV    13-30         .5-.7
% Class 4      2-17ms        -67 to -60  1.7-2.9mV   22-39        (stops)
% Class 3      30-70ms       -77 to -65  2.4-4.4mV   15-19         .7-1.1
% Class 5      18-58ms           ??      .3-2mV      5-25          .5-.9
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%         'SpikeAmp' 'SpikeWidth'  'ISI1'   'ThreshRate'
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz
   
% Most important IPs to examine:
%   - AHP_time2trough     'AHP(time2trough)'   
%   - ISI_median          'median(ISIs)'       
%   - AR23                'ISI2/ISI3'          
%   - RMP                 'RMP'
%   - Ih_abssag           'Ih'
% Others
%   - AP_amp              'SpikeAmp'          (59-68mV)
%   - AP_dur              'SpikeWidth'        (1.49-1.53ms)
%   - ISI1                'ISI1' (1/'1/ISI1') (31-39ms)
%   - FR_min              'ThreshRate'        (1-4Hz)

% -------------------------------------------------------------------------
% Step 1: qualitative search (look for approx regimes)
% Method: use local SimulateModel with 'vary' and PlotData
% -------------------------------------------------------------------------
cd /projectnb/crc-nak/sherfey/projects/ACC_rhythm_processing/model_files
tspan=[0 500];  % [ms], time limits on numerical integration
dt=.01;         % [ms], time step for numerical integration
solver='rk2';   % {options; 'rk4', 'rk2', 'rk1'}, solver to use
compile_flag=1; % whether to compile the simulation
verbose_flag=1; % whether to display log information
downsample_factor=10;
solver_options={'tspan',tspan,'dt',dt,'solver',solver,'compile_flag',compile_flag,'verbose_flag',verbose_flag,'downsample_factor',downsample_factor};
vary=[]; modifications=[];

% Sodium:     (iNaF,iNaP)
% Potassium:  (iKDR,iKS,iM,iKCa,iAHP)
% Calcium:    (iHVA,iCaH,iCaT,CaBuffer)
% Nonspecific: ih,ileak
Cm=1.2; % Durstewitz 2003

%            AHP duration              hyperpol-sag   <ISI>      adaptation
% TARGETS: [AHP_time2trough]   [RMP]   [Ih_abssag] [ISI_median]    [AR23]
%          'AHP(time2trough)'  'RMP'      'Ih'    'median(ISIs)' 'ISI2/ISI3'
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16          .7-1.1
% Class 2      11-32ms       -81 to -69  .5-1.4mV    13-30         .5-.7
% Class 4      2-17ms        -67 to -60  1.7-2.9mV   22-39        (stops)
% Class 3      30-70ms       -77 to -65  2.4-4.4mV   15-19         .7-1.1
% Class 5      18-58ms           ??      .3-2mV      5-25          .5-.9
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%         'SpikeAmp' 'SpikeWidth'  'ISI1'   'ThreshRate'
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz

% TARGETS: [AHP_time2trough]   [RMP]  [Ih_abssag]   [ISI_median]    [AR23]
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16ms        .7-1.1
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz
   
% Class 1 model (Group 2) version 1:
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
%eqns=['dV/dt=(@current+noise*randn)/Cm; Cm=1.2; noise=5; V(0)=-65;monitor functions;' mechanisms];
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1.2; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'gleak',.15;'gh',.025;'gCaT',1.5;'gKCa',.2;'Eleak',-75;'noise',5};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'gM',[.1];'gh',[.006 .007 .008 .009];'gAHP',[.4];'Iinj',[-2]});
% PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',(-10:10:400),'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1
% MODEL:    [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
% Class 1      3ms (refine)     -74       .6mV          8ms        .94 (ARif=.88)
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%          81mV       1.31ms      f(Inp)     1-2Hz

% Class 1 model (Group 2) version 2 (better):
mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa,iM}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1.2; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'akdr_scale',.75;'amnaf_scale',2;'gKDR',10;'gNaF',35;'gM',3;'gleak',.15;'gh',.025;'gCaT',0;'gKCa',.2;'Eleak',-75;'Iinj',15};
model=ApplyModifications(eqns,mods);
data=SimulateModel(model,solver_options{:},'vary',{'Iinj',[0 5 10 15 20 25]});
PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',[-10 -7.5 -5 0:10:400],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1
% MODEL:    [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
% Class 1      11.6ms (refine)     -75       .6mV         f(Inp)    .91 (ARif=.5)
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%          94mV       1.7ms      f(Inp)       8Hz

% Observation:
% (akdr_scale < amnaf_scale) -> gradual rebound after repolarization phase (Class 1)
% (akdr_scale > amnaf_scale) -> sharp rebound after repolarization phase (Class 2)

% Heterogeneity:
% gh: 0-.1?  (Ih: 0.885 +/- 0.66)
% Eleak: -85 to -65

% -------------------------------------------------------------------------
% Step 2: quantitative assessment of single IP = f(param); plot (param,IPs)
% -------------------------------------------------------------------------

paramdist=@(mu,sigma,num_vals)normrnd(mu,sigma,[1 num_vals]);

% HETEROGENEOUS gh:

num_cells=10;
sag_mu=0.885; sag_sd=0.66; level=2*(sag_sd/sag_mu);
gh_mu=.025; gh_sigma=level*gh_mu;
gh=max(0,paramdist(gh_mu,gh_sigma,num_cells));

mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1.2; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'size',num_cells;'gh',gh;'gCaT',1.5;'gKCa',.2;'gleak',.15;'Eleak',-75;'noise',5};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'gM',[.1];'gh',[.006 .007 .008 .009];'gAHP',[.4];'Iinj',[-2]});
% PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',(-10:10:400),'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

sim_sag_mu=mean(stats.pop1.Ih_abssag);
sim_sag_sd=std(stats.pop1.Ih_abssag);
fprintf('Exp Ih mu: %g, sd: %g\n',sag_mu,sag_sd);
fprintf('Sim Ih mu: %g, sd: %g\n',sim_sag_mu,sim_sag_sd);

figure; 
subplot(1,3,1); hist(gh)
subplot(1,3,2); hist(stats.pop1.Ih_abssag);
subplot(1,3,3); plot(gh,stats.pop1.Ih_abssag,'o-');

% plot all IPs, one per row
params=gh; param_label='gh';
IP_labels={'AHP_time2trough','RMP','Ih_abssag','ISI_median','AR23'};
IP_lims={[0 40],[-90 -60],[0 2],[0 20],[.5 1.5]};
figure; cnt=0; num_IPs=length(IP_labels);
for i=1:num_IPs
  subplot(num_IPs,2,cnt+1); 
  hist(stats.pop1.(IP_labels{i})); xlim(IP_lims{i});
  xlabel(IP_labels{i});
  subplot(num_IPs,2,cnt+2); 
  plot(params,stats.pop1.(IP_labels{i}),'o-'); 
  xlabel(param_label); ylabel(IP_labels{i}); ylim(IP_lims{i});
  cnt=cnt+2;
end

% HETEROGENEOUS Eleak:
num_cells=10;
RMP_mu=-74.9; RMP_sd=9.96; level=.75*(RMP_sd/RMP_mu);
Eleak_mu=-75; Eleak_sigma=level*Eleak_mu;
Eleak=paramdist(Eleak_mu,Eleak_sigma,num_cells);

mechanisms='{iNaF,iKDR,ileak,ih,CaBuffer,iCaT,iKCa}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1.2; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'size',num_cells;'gh',.025;'gCaT',1.5;'gKCa',.2;'gleak',.15;'Eleak',Eleak;'noise',5};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'gM',[.1];'gh',[.006 .007 .008 .009];'gAHP',[.4];'Iinj',[-2]});
% PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',(-10:10:400),'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

sim_RMP_mu=mean(stats.pop1.RMP);
sim_RMP_sd=std(stats.pop1.RMP);
fprintf('Exp RMP mu: %g, sd: %g\n',RMP_mu,RMP_sd);
fprintf('Sim RMP mu: %g, sd: %g\n',sim_RMP_mu,sim_RMP_sd);

figure; 
subplot(1,3,1); hist(Eleak)
subplot(1,3,2); hist(stats.pop1.RMP);
subplot(1,3,3); plot(Eleak,stats.pop1.RMP,'o-');

% plot all IPs, one per row
params=Eleak; param_label='Eleak';
IP_labels={'AHP_time2trough','RMP','Ih_abssag','ISI_median','AR23'};
IP_lims={[0 40],[-90 -60],[0 2],[0 20],[.5 1.5]};
figure; cnt=0; num_IPs=length(IP_labels);
for i=1:num_IPs
  subplot(num_IPs,2,cnt+1); 
  hist(stats.pop1.(IP_labels{i})); xlim(IP_lims{i});
  xlabel(IP_labels{i});
  subplot(num_IPs,2,cnt+2); 
  plot(params,stats.pop1.(IP_labels{i}),'o-'); 
  xlabel(param_label); ylabel(IP_labels{i}); ylim(IP_lims{i});
  cnt=cnt+2;
end

%% Class 2: ACC L5/6 PY (rat ACd: LeBeau lab) (~ Group 3)
%{
  Distribution across layers:
   superficial  middle    deep
      23%        23%      46%
  Note: Fiona Group 3 is not present in superficial layers...
%}
%{
Intrinsic Properties (mean +/- std) sorted by Heterogeneity (=std/mean)
Class2 (n=13):
HOMOGENEOUS PROPERTIES (std/mu<.4):
(0.085):              RMP: -75.4 +/- 6.44            (n=13)
(0.17):         SpikeAmp: 67.7 +/- 11.8            (n=13)
(0.22):       SpikeWidth: 1.49 +/- 0.321           (n=13)
(0.24):        ISI2/ISI3: 0.565 +/- 0.134           (n=9)
HETEROGENEOUS PROPERTIES:
(0.4):     median(ISIs): 21.7 +/- 8.76            (n=12)
(0.46):               Ih: 0.95 +/- 0.437           (n=12)
(0.54): AHP(time2trough): 21.4 +/- 11.7            (n=13)
(0.58):           1/ISI1: 26.7 +/- 15.4            (n=13)
(0.59):         Ih-decay:  113 +/- 66.5            (n=12)
(1.2):       ThreshRate: 1.69 +/- 1.99            (n=11)
%}

%            AHP duration              hyperpol-sag   <ISI>      adaptation
% TARGETS: [AHP_time2trough]   [RMP]   [Ih_abssag] [ISI_median]    [AR23]
%          'AHP(time2trough)'  'RMP'      'Ih'    'median(ISIs)' 'ISI2/ISI3'
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16          .7-1.1
% Class 2      11-32ms       -81 to -69  .5-1.4mV    13-30         .5-.7
% Class 4      2-17ms        -67 to -60  1.7-2.9mV   22-39        (stops)
% Class 3      30-70ms       -77 to -65  2.4-4.4mV   15-19         .7-1.1
% Class 5      18-58ms           ??      .3-2mV      5-25          .5-.9
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%         'SpikeAmp' 'SpikeWidth'  'ISI1'   'ThreshRate'
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz

mechanisms='{iNaF,iKDR,ileak,cadyn,can,iAHP,M,ih}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'hV1NaF',23.1;'hV2NaF',25.1;'akdr_scale',1;'gAHP',.054;'gleak',.017;'gM',.084;'gcan',.0056;'Eleak',-66;'taurCa',1.93;'taudCa',200/7;'gKDR',6;'gNaF',50.25;'CaRest',5e-5;'aAHP_scale',1e6};

% Class 2 model (Group 3)
mechanisms='{iNaF,iKDR,ileak,cadyn,can,iAHP,M,ih}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'hV1NaF',23.1;'hV2NaF',25.1;'gh',.001;'akdr_scale',1;'gAHP',.4;'gleak',.017;'gM',.1;'gcan',.0056;'Eleak',-70;'taurCa',1.93;'taudCa',200/7;'gKDR',6;'gNaF',50.25;'CaRest',5e-5;'aAHP_scale',1e6};
model=ApplyModifications(eqns,mods);
%data=SimulateModel(model,solver_options{:},'vary',{'gM',[.1];'gh',[.001];'gAHP',[.4];'Iinj',[15]});
%PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',(-10:10:400),'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1
% MODEL:    [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
% Class 2      26ms              -73mV     .97mV    14ms        .97 (ARif=.89)
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%          84mV       1.45ms       f(Inp)     1-2Hz

% Experiment:
% RMP: -75.4 +/- 6.44            (n=13)
% Ih: 0.95 +/- 0.437           (n=12)

% Model heterogeneity:
% gh: 0-.002
% Eleak: -80 to -70

paramdist=@(mu,sigma,num_vals)normrnd(mu,sigma,[1 num_vals]);

% HETEROGENEOUS gh:

num_cells=10;
sag_mu=0.95; sag_sd=0.437; level=2*(sag_sd/sag_mu);
gh_mu=.001; gh_sigma=level*gh_mu;
gh=max(0,paramdist(gh_mu,gh_sigma,num_cells));

mechanisms='{iNaF,iKDR,ileak,cadyn,can,iAHP,M,ih}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'size',num_cells;'hV1NaF',23.1;'hV2NaF',25.1;'gh',gh;'akdr_scale',1;'gAHP',.4;'gleak',.017;'gM',.1;'gcan',.0056;'Eleak',-70;'taurCa',1.93;'taudCa',200/7;'gKDR',6;'gNaF',50.25;'CaRest',5e-5;'aAHP_scale',1e6};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'gM',[.1];'gh',[.006 .007 .008 .009];'gAHP',[.4];'Iinj',[-2]});
% PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',(-10:10:400),'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

sim_sag_mu=mean(stats.pop1.Ih_abssag);
sim_sag_sd=std(stats.pop1.Ih_abssag);
fprintf('Exp Ih mu: %g, sd: %g\n',sag_mu,sag_sd);
fprintf('Sim Ih mu: %g, sd: %g\n',sim_sag_mu,sim_sag_sd);

figure; 
subplot(1,3,1); hist(gh)
subplot(1,3,2); hist(stats.pop1.Ih_abssag);
subplot(1,3,3); plot(gh,stats.pop1.Ih_abssag,'o-');

% plot all IPs, one per row
params=gh; param_label='gh';
IP_labels={'AHP_time2trough','RMP','Ih_abssag','ISI_median','AR23'};
IP_lims={[0 40],[-90 -60],[0 2],[0 20],[.5 1.5]};
figure; cnt=0; num_IPs=length(IP_labels);
for i=1:num_IPs
  subplot(num_IPs,2,cnt+1); 
  hist(stats.pop1.(IP_labels{i})); xlim(IP_lims{i});
  xlabel(IP_labels{i});
  subplot(num_IPs,2,cnt+2); 
  plot(params,stats.pop1.(IP_labels{i}),'o-'); 
  xlabel(param_label); ylabel(IP_labels{i}); ylim(IP_lims{i});
  cnt=cnt+2;
end

% HETEROGENEOUS Eleak:
num_cells=10;
RMP_mu=-75.4; RMP_sd=6.44; level=.75*(RMP_sd/RMP_mu);
Eleak_mu=-70; Eleak_sigma=level*Eleak_mu;
Eleak=paramdist(Eleak_mu,Eleak_sigma,num_cells);

mechanisms='{iNaF,iKDR,ileak,cadyn,can,iAHP,M,ih}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'size',num_cells;'hV1NaF',23.1;'hV2NaF',25.1;'gh',.001;'akdr_scale',1;'gAHP',.4;'Eleak',Eleak;'gleak',.017;'gM',.1;'gcan',.0056;'taurCa',1.93;'taudCa',200/7;'gKDR',6;'gNaF',50.25;'CaRest',5e-5;'aAHP_scale',1e6};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'gM',[.1];'gh',[.006 .007 .008 .009];'gAHP',[.4];'Iinj',[-2]});
% PlotData(data)
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',(-10:10:400),'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

sim_RMP_mu=mean(stats.pop1.RMP);
sim_RMP_sd=std(stats.pop1.RMP);
fprintf('Exp RMP mu: %g, sd: %g\n',RMP_mu,RMP_sd);
fprintf('Sim RMP mu: %g, sd: %g\n',sim_RMP_mu,sim_RMP_sd);

figure; 
subplot(1,3,1); hist(Eleak)
subplot(1,3,2); hist(stats.pop1.RMP);
subplot(1,3,3); plot(Eleak,stats.pop1.RMP,'o-');

% plot all IPs, one per row
params=Eleak; param_label='Eleak';
IP_labels={'AHP_time2trough','RMP','Ih_abssag','ISI_median','AR23'};
IP_lims={[0 40],[-90 -60],[0 2],[0 20],[.5 1.5]};
figure; cnt=0; num_IPs=length(IP_labels);
for i=1:num_IPs
  subplot(num_IPs,2,cnt+1); 
  hist(stats.pop1.(IP_labels{i})); xlim(IP_lims{i});
  xlabel(IP_labels{i});
  subplot(num_IPs,2,cnt+2); 
  plot(params,stats.pop1.(IP_labels{i}),'o-'); 
  xlabel(param_label); ylabel(IP_labels{i}); ylim(IP_lims{i});
  cnt=cnt+2;
end

data(26).pop1_TONIC % Iinj=16
PlotData(data(26))


%% Class 3:

data=SimulateModel('dV/dt=@current-.4*(V+80)+I; {RSiNaF,RSiKDR,iCaT,iKCa,CaBuffer}; V(0)=-85;','vary',{'','I',[4 6];'','gKCa',0:5:20},'tspan',[0 500]);
data=SimulateModel('dV/dt=@current-.4*(V+80)+I; {RSiNaF,RSiKDR,iCaT,iKCa,CaBuffer}; V(0)=-85;','vary',{'','I',[4 6];'','gKCa',0:10:50},'tspan',[0 500]);
PlotData(data)

%            AHP duration              hyperpol-sag   <ISI>      adaptation
% TARGETS: [AHP_time2trough]   [RMP]   [Ih_abssag] [ISI_median]    [AR23]
%          'AHP(time2trough)'  'RMP'      'Ih'    'median(ISIs)' 'ISI2/ISI3'
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16          .7-1.1
% Class 2      11-32ms       -81 to -69  .5-1.4mV    13-30         .5-.7
% Class 3      30-70ms       -77 to -65  2.4-4.4mV   15-19         .7-1.1
% Class 4      2-17ms        -67 to -60  1.7-2.9mV   22-39        (stops)
% Class 5      18-58ms           ??      .3-2mV      5-25          .5-.9
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%         'SpikeAmp' 'SpikeWidth'  'ISI1'   'ThreshRate'
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz

% Class 3 model (Group 1):
% (RSiKDR,RSiNaF,iM,ileak) w/ (gM=10-50,Iinj=5-10) | Cm=1.2, gleak=.4, Eleak=-75, gKDR=20, gNaF=200
mechanisms='{RSiKDR,RSiNaF,iM,ileak,ih}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1; Iinj=2; V(0)=-65;monitor functions;' mechanisms];
mods={'gNaF',200;'gKDR',20;'gleak',.4;'Eleak',-75;'Cm',1.2;'gM',20;'gh',1.5};
% mods={'gNaF',125;'gKDR',5;'gleak',.4;'Eleak',-75;'Cm',1.2;'gM',20;'gh',1.5};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'Iinj',5:10});
% PlotData(data,'ylim',[-100 50])
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',[-10 -7.5 -5 0:10:400],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

% MODEL:    [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
% Class 3      9ms [x]          -70mV     2.5mV       37ms [x]     .89 (ARif=.25)
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%          103mV       .42ms [x]    f(Inp)     1-5Hz

%% Class 4: ACC L5/6 "reset" (rat ACd: LeBeau lab) (~ Group 4)
%{
  Distribution across layers:
   superficial  middle    deep
      33%        33%      33%
  Note: twice as many Group 4 cells in deep layer
%}
%{
Intrinsic Properties (mean +/- std) sorted by Heterogeneity (=std/mean)
Class4 (n=6):
HOMOGENEOUS PROPERTIES (std/mu<.4):
(0.052):              RMP: -63.9 +/- 3.34            (n=6)
(0.074):         SpikeAmp: 58.8 +/- 4.37            (n=3)
(0.26):               Ih: 2.27 +/- 0.59            (n=6)
(0.28):       SpikeWidth: 1.53 +/- 0.421           (n=3)
(0.3):     median(ISIs): 30.5 +/- 9.08            (n=6)
HETEROGENEOUS PROPERTIES:
(0.58):         Ih-decay: 92.2 +/- 53.5            (n=6)
(0.67):           1/ISI1: 25.7 +/- 17.3            (n=6)
(0.78): AHP(time2trough):  9.8 +/- 7.64            (n=3)
(1.3):       ThreshRate: 2.57 +/- 3.28            (n=3)
(NaN):        ISI2/ISI3:  NaN +/- NaN             (n=0)
%}
%{
Class 4 has higher RMP (-64>-75) and greater Ih (2.3>.9) than classes 1 and 2.
%}

%            AHP duration              hyperpol-sag   <ISI>      adaptation
% TARGETS: [AHP_time2trough]   [RMP]   [Ih_abssag] [ISI_median]    [AR23]
%          'AHP(time2trough)'  'RMP'      'Ih'    'median(ISIs)' 'ISI2/ISI3'
% Class 1      20-40ms       -85 to -65  .3-1.5mV    8-16          .7-1.1
% Class 2      11-32ms       -81 to -69  .5-1.4mV    13-30         .5-.7
% Class 3      30-70ms       -77 to -65  2.4-4.4mV   15-19         .7-1.1
% Class 4      2-17ms        -67 to -60  1.7-2.9mV   22-39        (stops)
% Class 5      18-58ms           ??      .3-2mV      5-25          .5-.9
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%         'SpikeAmp' 'SpikeWidth'  'ISI1'   'ThreshRate'
%          59-68mV   1.49-1.53ms   31-39ms     1-4Hz

% Class 4 model (Group 4):
% (iNaF,iKDR,iHVA,iAHP,ileak,CaBuffer) w/ (gHVA~2, gAHP~10) 	| Cm=1.2, gleak=.4, Eleak=-75, gKDR=4, gNaF=50, Iinj=35

mechanisms='{iNaF,iKDR,ileak,CaBuffer,iHVA,iAHP,ih}';
eqns=['dV/dt=(@current+Iinj*(t>50))/Cm; Cm=1.2; Iinj=2; V(0)=-75;monitor functions;' mechanisms];
mods={'gNaF',75;'gKDR',4;'gleak',.4;'Eleak',-65;'Cm',1.2;'gHVA',2;'gAHP',10;'gh',10};
model=ApplyModifications(eqns,mods);
% data=SimulateModel(model,solver_options{:},'vary',{'Iinj',-10:10:40});
% PlotData(data,'ylim',[-100 50])
modl=ApplyModifications(model,{'Iinj',0});
data=ProbeCellProperties(modl,'amplitudes',[-10 -7.5 -5 0:10:400],'tspan',[0 750],'compile_flag',1,'verbose_flag',1,'solver','rk2');
PlotData(data(unique(round(linspace(1,length(data),10)))),'ylim',[-100 50])
stats=CalcCellProperties(data,'plot_flag',0); 
stats.pop1

% MODEL:    [AHP_time2trough]   [RMP]   [Ih_abssag]  [ISI_median]    [AR23]
% Class 4      4.5ms          -62mV     1.7mV       11ms [x]     (stops)
% Common:  [AP_amp]    [AP_dur]    [ISI1]     [FR_min]
%          82mV       1.5ms       f(Inp)     1-2Hz

% only classes 1 and 2 have heterogeneity, and only in gh and Eleak.

% Group 5 model:
% (iNaF,iKDR,iHVA,iAHP,ileak,CaBuffer) w/ (gHVA~2, gAHP~10) 	| Cm=1.2, gleak=.4, Eleak=-75, gKDR=4, gNaF=50, Iinj=35
