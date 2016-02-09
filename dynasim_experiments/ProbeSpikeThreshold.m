function [threshold_amps,threshold_V,threshold_rate]=ProbeSpikeThreshold(model,varargin)
%% [threshold_amps,threshold_V]=ProbeSpikeThreshold(model,varargin)
% purpose: run experiment delivering varying levels of tonic input to a
% single cell model. determine minimum drive to produce spike.
% inputs:
%   model - DynaSim model structure
%   options:
%     'amplitudes' - tonic input amplitudes to test
%     'detection_threshold' - threshold for detecting spikes
%     (any options for SimulateModel)
%     (any options for CalcFR)
% outputs:
%   threshold_amps: array of minimum tonic drive strengths for spiking
% 
% Example:
% % define single cell model
% eqns='dv/dt=@current+input; input=0; {iNa,iK}';
% % perform threshold detection experiment
% [thresh_amp,thresh_v] = ProbeSpikeThreshold(eqns)
% % ---------------
% % confirm result:
% data=SimulateModel(s,'vary',{'','input',0:thresh_amp+2});
% % plot waveforms
% PlotData(data)
% % add line at voltage threshold
% line(xlim,[thresh_v thresh_v])
% % plot spike rates vs input amplitude
% PlotFR(data) % threshold amp: 5
% % add line at amplitude threshold
% line([thresh_amp thresh_amp],ylim);
% 
% see also: SimulateModel, ProbeRhythmThreshold

options=CheckOptions(varargin,{...
  'target','ODE1',[],...
  'amplitudes',0:10:100,[],...
  'detection_threshold',5,[],...
  'skip_time',10,[],... % time [ms] to exclude from detection
  'tspan',[0 200],[],...
  },false);

model=CheckModel(model);
orig_model=model;

npops=length(model.specification.populations);
pop_sizes=[model.specification.populations.size];
if npops>1 || any(pop_sizes)>1
  error('experiment only works on a single cell model.');
end

% add input to model
modifications=cell(npops,3);
for i=1:npops
  name=model.specification.populations(i).name;
  % prepare list of modifications to add tonic drive to all populations in model
  modifications(i,:)={name,'equations',['cat(' options.target ',+TONIC)']};
end
% apply modifications
model=ApplyModifications(model,modifications);

threshold_amps=nan;
threshold_V=nan;
threshold_rate=nan;

% vary amplitude
fprintf('testing amplitudes: [%s]\n',num2str(options.amplitudes));
for i=1:length(options.amplitudes)
  amp=options.amplitudes(i);
  data=SimulateModel(model,'modifications',{name,'TONIC',amp},varargin{:});
  data=SelectData(data,'time_limits',[options.skip_time options.tspan(2)]);
  % check for spiking
  dat=data.(data.labels{1});
  if any(dat(:)>options.detection_threshold)
    fprintf('found spike at amplitude=%g\n',amp);
    % check if this is the minimum positive drive
    if amp==min(options.amplitudes(options.amplitudes>0))
      % test smaller inputs
      fprintf('testing smaller inputs...\n');
      options.amplitudes=options.amplitudes/10;
      keyvals=Options2Keyval(options);
      [threshold_amps,threshold_V,threshold_rate]=ProbeSpikeThreshold(orig_model,keyvals{:});
      if ~isnan(threshold_amps)
        % use the returned threshold_amps
        return;
      end
    end
    % use this threshold_amps
    fprintf('returning threshold amplitude = %g\n',amp);
    threshold_amps=amp;
    data=CalcFR(data,'bin_size',30,'bin_shift',10);
    threshold_rate=nanmean(data.(data.results{1}),1);
    return;
  end
  threshold_V=max(dat(:));  
end

