

model='dv/dt=(@current-gleak*(v+70))/Cm; gleak=.04; Cm=2; {iNa,iK}';
Rleak=1/.04; Cm=2; taum_model=Rleak*Cm;
data=ProbeCellProperties(model,'num_repetitions',1);

model='dv/dt=(@current-gleak*(v+70))/Cm; gleak=.1; Cm=1; {iNa,iK}';
Rleak=1/.1; Cm=1; taum_model=Rleak*Cm;
data=ProbeCellProperties(model,'num_repetitions',2);
stats=CalcCellProperties(data); stats.pop1

model='dv/dt=(@current-gleak*(v+70)+5*randn)/Cm; gleak=.1; Cm=1; {iNa,iK}';
data=ProbeCellProperties(model,'num_repetitions',10);
stats=CalcCellProperties(data); stats.pop1

model='dv/dt=(@current-gleak*(v+70)+5*randn)/Cm; gleak=.1; Cm=1; {iNa,iK}';
data=ProbeCellProperties(model,'num_repetitions',2);
stats=CalcCellProperties(data); stats.pop1

Rleak=1/.1; Cm=1; taum_model=Rleak*Cm;
R_in=stats.pop1.R_in; tau_m=stats.pop1.tau_m;
[Rleak,R_in],[taum_model,tau_m]

model='dv[2]/dt=(@current-gleak*(v+70)+5*randn(1,N_pop))/Cm; gleak=.1; Cm=1; {iNa,iK}';
data=ProbeCellProperties(model,'num_repetitions',2);
stats=CalcCellProperties(data); stats.pop1

% relationship to my names for Tallie's IPs:
% FR_step = steprate
% FR_min = threshrate
% Ih_sag = Ih
% AP_dur = spikewidth
% AP_amp = SpikeAmp

%% TEST (@ProbeCellProperties,@CalcCellProperties) w/ vary on cluster & local
model='dv/dt=(@current-gleak*(v+70)+5*randn)/Cm; gleak=.1; Cm=1; {iNa,iK}';
% Cluster test
study_dir='test_IPs_cluster';
pop='pop1'; Cm_vals=[1 2 3];
vary={pop,'Cm',Cm_vals};
experiment=@ProbeCellProperties; analysis=@CalcCellProperties;
[~,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir,...
  'experiment',experiment,'analysis_functions',analysis,...
  'cluster_flag',1,'save_data_flag',1,'tspan',[0 1500],'verbose_flag',1);
!qstat -u sherfey
!cat ~/batchdirs/test_IPs/pbsout/sim_job1.out

% Local test
study_dir='test_IPs_local';
vary={pop,'Cm',[1 2 3]};
experiment=@ProbeCellProperties; analysis=@CalcCellProperties;
[~,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir,...
  'experiment',experiment,'analysis_functions',analysis,...
  'cluster_flag',0,'save_data_flag',1,'tspan',[0 1500],'verbose_flag',1);

% Check results
stats=ImportResults(study_dir,analysis);
tau_m =arrayfun(@(x)x.(pop).tau_m,stats);
RMP   =arrayfun(@(x)x.(pop).RMP,stats);
%Cm_vals=[stats.([pop '_Cm'])];
figure; 
subplot(2,1,1); plot(Cm_vals,tau_m,'o-'); xlabel('Cm'); ylabel('tau_m [ms]');
subplot(2,1,2); plot(Cm_vals,RMP,'o-'); xlabel('Cm'); ylabel('RMP [mV]');

data=ImportData(study_dir);
dat=SelectData(data,'varied',{'pop1','Cm',1});
dat=SelectData(dat,'varied',{'pop1','TONIC',[1 4]});
PlotData(dat);
PlotFR(data)

% Cluster test saving results but not simulated data
study_dir='test_IPs_cluster_2';
pop='pop1'; Cm_vals=[1 2 3]; vary={pop,'Cm',Cm_vals};
experiment=@ProbeCellProperties; analysis=@CalcCellProperties;
[~,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir,...
  'experiment',experiment,'analysis_functions',analysis,...
  'cluster_flag',1,'save_data_flag',0,'tspan',[0 1500],'verbose_flag',1);
!qstat -u sherfey
!cat ~/batchdirs/test_IPs_cluster_2/pbsout/sim_job1.out
stats=ImportResults(study_dir,analysis);
tau_m =arrayfun(@(x)x.(pop).tau_m,stats);
figure; plot(Cm_vals,tau_m,'o-'); xlabel('Cm'); ylabel('tau_m [ms]');

% Cluster test with population (size>1) and num_repetitions>1
model='dv[2]/dt=(@current-gleak*(v+70)+5*randn(1,N_pop))/Cm; gleak=.1; Cm=1; {iNa,iK}';
study_dir='test_IPs_cluster_3';
pop='pop1'; Cm_vals=[1 2 3]; vary={pop,'Cm',Cm_vals};
experiment=@ProbeCellProperties; analysis=@CalcCellProperties;
experiment_options={'num_repetitions',2};
[~,studyinfo]=SimulateModel(model,'vary',vary,'study_dir',study_dir,...
  'experiment',experiment,'experiment_options',experiment_options,'analysis_functions',analysis,...
  'cluster_flag',1,'save_data_flag',0,'tspan',[0 1500],'verbose_flag',1,'solver','rk2');
!qstat -u sherfey
!cat ~/batchdirs/test_IPs_cluster_3/pbsout/sim_job1.out
stats=ImportResults(study_dir,analysis);
stats(1).experiment_options
stats(1).pop1


%% Complete example using the cluster:
% turn off tonic input and poisson noise
spec_=ApplyModifications(spec,{'(E,If,Is)','Iapp',0;'(E,If,Is)','gNOISE',0});
% define biophysical parameter ranges to explore
vary={'E','Eleak',[-80 -70 -60]}; % 'E','gCaT',[0 1];'E','gAHP',[0 1]
% specify experiment and post-experiment analysis to probe cell properties
experiment=@ProbeCellProperties;
analysis=@CalcCellProperties;
analysis_options={'option1',option1,'option2',option2};
% specify plots to generate after each experiment
plot_functions={@PlotData,@PlotData,@PlotData};
plot_options={{'plot_type','waveform'},{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
% run all simulations on cluster
[~,studyinfo]=SimulateModel(spec_,'vary',vary,'experiment',experiment,...
  'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
  'plot_functions',plot_functions,'plot_options',plot_options,...
  solver_options{:},cluster_options{:},output_options{:});
% import and plot all data across experiments
data=ImportData(studyinfo); % [1 x num_sims]
PlotData(data,'variable','E_v'); % plot voltage response
PlotData(data,'variable','E_INPUT'); % plot inputs
% import cell properties across experiments
stats=ImportResults(studyinfo,analysis); % [num_experiments x num_calls]
  % i.e., (ephys properties for all cells) for each cell model
% select and plot data for one experiment (corresponding to one element of stats)
dat=SelectData(data,'varied',{'E','Eleak',-80});
dat=SelectData(data,'varied',stats(1).modifications(1,:));
PlotData(dat,'variable','E_v'); stats(1)
% parameters varied across experiments (corresponding to stats)
[P1,L1,V1]=CollectVariedParams(stats);
  % P = all_values [num_sims x num_varied]
  % L = param names {1 x num_varied}
  % V = unique values {1 x num_varied}
% parameters varied within and across experiments (corresponding to data)
[P2,L2,V2]=CollectVariedParams(data);

% to support this workflow, need to:
% - [DONE] modify SimulateModel: debug saving data from experiments; save all data
%   returned by experiment to a single file, with one entry in
%   studyinfo.simulations giving modifications based on 'vary' option and not
%   those imposed by the experiment. for each simulation iteration, update
%   tmpdata(:).simulator_options.modifications using the modifications for
%   the given iteration.
% - [DONE] modify ImportData: adjust data indexing to support concatenating data
%   files containing multiple data sets.
% - [DONE] modify AnalyzeData to allow numel(data)>1, pass all data to the
%   analysis function, and add 'varied' to the results structure based on
%   data(1).simulator_options.modifications (instead of data.varied). also
%   store stats(i).modifications=data(1).simulator_options.modifications.

% note: data(i).varied contains info on all parameters varied within and
% across experiments. data(i).simulator_options.modifications and
% studyinfo.simulations(j).modifications contains only the parameters
% varied across experiments.

% when to use experiments: (copy to tutorial)
% if you have an analysis that needs to be applied to the data from a set
% of simulations and you want to examine how the analysis results change as
% some component of the model is varied. that is, for every version of the
% model, an experiment produces multiple data sets, the analysis collapses
% them into a single result, and one result is saved for each version of
% the model. example) an experiment that varies the amplitude of an input
% across a set of simulations, followed by an analysis that computes
% system properties based on responses across the drives (eg, slope of 
% f/I curve), and then you want to examine how system properties change
% with model components (eg, how membrane capacitance affects the slope of
% the f/I curve).

% test experiment wrapped by varied components
[data,studyinfo]=SimulateModel('dv/dt=@current+Iext; {iNa,iK}; Iext=0',...
  'experiment',@ProbeFI,'vary',{'','Iext',[0 10]},'verbose_flag',1);
PlotData(data);

% test experiment wrapped by varied components and data saved to disk
[data,studyinfo]=SimulateModel('dv/dt=@current+Iext; {iNa,iK}; Iext=0',...
  'experiment',@ProbeFI,'vary',{'','Iext',[0 10]},'study_dir','test_exp','save_data_flag',1,'verbose_flag',1);
d=ImportData(studyinfo);
PlotData(d)

% cluster test experiment wrapped by varied components and data saved to disk
[data,studyinfo]=SimulateModel('dv/dt=@current+Iext; {iNa,iK}; Iext=0','cluster_flag',1,...
  'experiment',@ProbeFI,'vary',{'','Iext',[0 10]},'study_dir','test_exp','save_data_flag',1,'verbose_flag',1);
d=ImportData(studyinfo);
PlotData(d)

% 1. [DONE] Create minimal ProbeCellProperties and CalcCellProperties for workflow dev
% 2. [DONE] Edit and test serial AnalyzeData manually:
%    [DONE] stats=AnalyzeData(ProbeCellProperties(m),@CalcCellProperties)
%    [DONE] stats=AnalyzeData(SimulateModel(m,'experiment',@ProbeCellProperties))
% 3. [DONE] Test cluster AnalyzeData through SimulateModel
%    [DONE] SimulateModel(m,'experiment',@ProbeCellProperties,'analysis_functions',@CalcCellProperties,'cluster_flag',1);
% 4. Create full ProbeCellProperties and CalcCellProperties
% 5. Use to fit cell models to IPs from Tallie's data and literature

% test serial AnalyzeData manually:
data=ProbeCellProperties(s,'verbose_flag',1,'solver','rk2');
PlotData(data,'variable','E_v','ylim',[-150 50]);
stats = CalcCellProperties(data)
stats = AnalyzeData(data,@CalcCellProperties,'plot_flag',0)

data=SimulateModel(s,'experiment',@ProbeCellProperties,'verbose_flag',1,'solver','rk2')
stats=AnalyzeData(data,@CalcCellProperties,'plot_flag',0)

% test cluster AnalyzeData through SimulateModel saving data and results
[~,studyinfo]=SimulateModel(s,'experiment',@ProbeCellProperties,'verbose_flag',1,'solver','rk2',...
  'analysis_functions',@CalcCellProperties,'cluster_flag',1,...
  'study_dir','test_exp','save_data_flag',1,'tspan',[0 250]);
data=ImportData(studyinfo)
stats=ImportResults(studyinfo,@CalcCellProperties)

% test cluster AnalyzeData through SimulateModel with 'vary' and saving data and results
vary={'E','Cm',[.5 2]};
[~,studyinfo]=SimulateModel(s,'experiment',@ProbeCellProperties,'vary',vary,'verbose_flag',1,'solver','rk2',...
  'analysis_functions',@CalcCellProperties,'cluster_flag',1,...
  'study_dir','test_exp','save_data_flag',1,'tspan',[0 250]);
data=ImportData(studyinfo)
PlotData(data,'variable','E_v','ylim',[-150 50]);
stats=ImportResults(studyinfo,@CalcCellProperties)
figure
subplot(1,2,1); plot(stats(1).amplitudes,stats(1).E_v_mean_per_amp); ylim([-100 -50])
subplot(1,2,2); plot(stats(2).amplitudes,stats(2).E_v_mean_per_amp); ylim([-100 -50])

% test cluster AnalyzeData through SimulateModel saving results only
vary={'E','Cm',[.5 2]};
[~,studyinfo]=SimulateModel(s,'experiment',@ProbeCellProperties,'vary',vary,'verbose_flag',1,'solver','rk2',...
  'analysis_functions',@CalcCellProperties,'cluster_flag',1,...
  'study_dir','test_exp','save_data_flag',0,'tspan',[0 250]);
stats=ImportResults(studyinfo,@CalcCellProperties)
figure
subplot(1,2,1); plot(stats(1).amplitudes,stats(1).E_v_mean_per_amp); ylim([-100 -50])
subplot(1,2,2); plot(stats(2).amplitudes,stats(2).E_v_mean_per_amp); ylim([-100 -50])

% full example (saving results and data)
vary={'E','Cm',[.5 2]};
experiment=@ProbeCellProperties;
analysis_functions=@CalcCellProperties;
analysis_options={'plot_flag',0};
plot_functions={@PlotData,@PlotData,@PlotData};
plot_options={{'plot_type','waveform','ylim',[-100 0]},{'plot_type','rastergram'},{'plot_type','power','xlim',[0 100]}};
[~,studyinfo]=SimulateModel(s,'experiment',experiment,'vary',vary,'verbose_flag',1,'solver','rk2',...
  'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
  'plot_functions',plot_functions,'plot_options',plot_options,...
  'cluster_flag',1,'study_dir','test_exp','save_data_flag',1,'tspan',[0 250]);
data=ImportData(studyinfo);
PlotData(data,'variable','E_pulse')               % plot inputs
PlotData(data,'variable','E_v','ylim',[-150 50]); % plot responses
stats=ImportResults(studyinfo,analysis_functions);
figure
subplot(1,2,1); plot(stats(1).amplitudes,stats(1).E_v_mean_per_amp); ylim([-100 -50])
subplot(1,2,2); plot(stats(2).amplitudes,stats(2).E_v_mean_per_amp); ylim([-100 -50])
stats(1)
dat=SelectData(data,'varied',stats(1).modifications(1,:)); 
%dat=SelectData(data,'varied',{'E','Cm',2});
%dat=SelectData(data,'varied',{'E_Cm',2});
PlotData(dat,'variable','E_v','ylim',[-150 50]); 
[P1,L1,V1]=CollectVariedParams(stats);
[P2,L2,V2]=CollectVariedParams(data);

% full example (saving results only)
vary={'E','Cm',[.5 2]};
experiment=@ProbeCellProperties;
experiment_options={'amplitudes',0:2:18};
analysis_functions=@CalcCellProperties;
analysis_options={'plot_flag',0};
plot_functions={@PlotData,@PlotData,@PlotData,@PlotData};
plot_options={
  {'plot_type','waveform','variable','E_pulse'},
  {'plot_type','waveform','variable','v','ylim',[-100 0]},
  {'plot_type','rastergram'},
  {'plot_type','power','xlim',[0 100]}};
[~,studyinfo]=SimulateModel(s,'vary',vary,'cluster_flag',1,'tspan',[0 250],'solver','rk2',...
  'experiment',experiment,'experiment_options',experiment_options,...
  'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
  'plot_functions',plot_functions,'plot_options',plot_options,...
  'study_dir','test_exp','save_data_flag',0,'verbose_flag',1);
stats=ImportResults(studyinfo,analysis_functions);
figure
subplot(1,2,1); plot(stats(1).amplitudes,stats(1).E_v_mean_per_amp); ylim([-100 -50])
subplot(1,2,2); plot(stats(2).amplitudes,stats(2).E_v_mean_per_amp); ylim([-100 -50])
stats(1)
[P1,L1,V1]=CollectVariedParams(stats);







