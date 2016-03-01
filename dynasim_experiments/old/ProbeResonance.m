function [data,stats]=ProbeResonance(model,varargin)
% Experiment: data=ProbeResonance(model,'option',value,...)
% Deliver one rhythmic input, varying the frequency of the inputs over
% range in order to determine what frequency generates the maximal
% response. The input may be a direct current injection (sinusoid, noisy 
% sinusoid, rectified sinusoid, or gaussian smoothed pulses) or 
% nonhomogeneous Poisson process filtered through double-exponential AMPA synapses.
% see also: ProbeTwoRhythms

% Check options
options=CheckOptions(varargin,{...
  'input_type','poisson',{'sin','noisy_sin','rectified_sin','poisson','noisy_rectified_sin'},...
  'target',[],[],... % target populations; default: first population only
  'gINPUT',1,[],... % max amplitude of applied current
  'f',5:5:60,[],... % Hz, modulation frequency (0: homogeneous poisson process)
  'num_repetitions',1,[],... % number of realizations of the input (set to >1 for noisy inputs if stats are desired)
  'dc',0,[],... % kHz
  'ac',1,[],... % kHz
  'baseline',.1,[],... % kHz
  'tau',2,[],... % ms
  'bin_size',30,[],... % ms
  'bin_shift',10,[],... % ms
  'plot_flag',1,[],...
  },false);

options.realization=1:options.num_repetitions;
% Check model
model=CheckModel(model);
% Check targets
pop_names={model.specification.populations.name};
if isempty(options.target)
  options.target=pop_names{1};
elseif ischar(options.target)
  options.target={options.target};
end
npops=length(pop_names);

% Add input to model and set 'vary' (f,realization; other options)
modifications={};
vary={};
for i=1:npops
  name=model.specification.populations(i).name;
  N=model.specification.populations(i).size;
  if ~ismember(name,options.target)
    % add nothing
    continue;
  end
  % prepare list of modifications to add tonic drive to all populations in model
  switch options.input_type
    case 'poisson'
      % hack: get ID to add to avoid conflict from changing number varied
      param_names={'gINPUT','dc','ac','tau','baseline','f','realization'};
      ID=[];
      for p=1:length(param_names)
        if numel(options.(param_names{p}))>1
          ID=[ID p];
        end
      end
      % add AMPA synapse for exponentially filtered modulated poisson process
      s=sprintf('s=get_input(''%s'',Npop,T,f,dc,ac,tau,xc,baseline,phase); xc=.5; phase=0',options.input_type);
      eqn_mods=sprintf('cat(ODE1,-gINPUT*s(k,:).*(X-0)+0*%g; %s)',sum(ID),s);
      modifications(end+1,:)={name,'equations',eqn_mods};
      for p=1:length(param_names)
        param=param_names{p};
        if numel(options.(param))==1
          modifications(end+1,:)={name,param,options.(param)};
        else
          vary(end+1,:)={name,param,options.(param)};
        end
      end
    case {'sin','noisy_sin','rectified_sin','noisy_rectified_sin'}
      % LEGACY:
      if i==1, mod_cnt=1; vary_cnt=1; end
      % add +sin() to equations and f to parameters
      I=sprintf('I=get_input(''%s'',Npop,T,f)',options.input_type);
      eqn_mods=sprintf('cat(ODE1,+gINPUT*I(k,:); %s)',I);
      modifications(mod_cnt+0,:)={name,'equations',eqn_mods};
      modifications(mod_cnt+1,:)={name,'gINPUT',options.gINPUT};
      % prepare specification to vary tonic drive across simulations
      vary(vary_cnt+0,:)={name,'f',options.f};
      vary(vary_cnt+1,:)={name,'realization',1:options.num_repetitions};
      mod_cnt=mod_cnt+2;
      vary_cnt=vary_cnt+2;
    otherwise
      error('input_type not recognized.');
  end
end

fprintf('Running experiment: %s\n',mfilename);

% apply modifications to effectively add experimental apparatus to model
model=ApplyModifications(model,modifications);

% execute experimental protocol by varying parameters across simulations
data=SimulateModel(model,'vary',vary,varargin{:});

%% post-processing
% Calculate statistics across frequencies and realizations
if nargout>1 || options.plot_flag==1
  fprintf('Postprocessing: calculating resonance and spectral properties...\n');
  
  realization_param=[options.target{1} '_realization'];
  freq_param=[options.target{1} '_f'];

  [stats,data]=CalcPopulationStats(data,...
    'sweep_parameter',freq_param,...       % max over f
    'repetition_parameter',realization_param,... % mean/std over realizations
    varargin{:});
end

if options.plot_flag && options.num_repetitions>1
  % ---------------------------------------------------
  % Resonance Plots
  % ---------------------------------------------------
  if isfield(stats,'FR_variables') && numel(stats.(stats.FR_variables{1}).repetition_sets.varied)==1
    var=stats.FR_variables{1};
    if isfield(stats.(var),'sweep_sets')
      N=stats.num_repetitions;
      sweep_parameter=stats.(var).sweep_sets.sweep_parameter;
      f=stats.(var).sweep_sets.(sweep_parameter);
      p=stats.(var).sweep_sets.repetition_sets.parameters;
      FRmu=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_mu;
      FRsd=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_sd;
      fcmu=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_mu;
      fcsd=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_sd;

      % plot Firing Rate Tuning Curves
      colors='bgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmky';
      figure('position',[200 340 1500 550])
      subplot(1,2,1);
      for i=1:length(p)
        color=colors(i);
        x=f;
        y=FRmu(:,i);
        e=FRsd(:,i)/sqrt(N);
        plot_CI(x,y,e,color);
        hold on
      end
      xlabel('Input Frequency')
      ylabel([strrep(var,'_','\_') ' (population mean) [Hz]']);  
      % add legend
      pname=strrep(stats.(var).repetition_sets.varied{1},'_','\_');
      legend(cellfun(@(x)[pname '=' num2str(x)],num2cell(p),'uni',0))
      % add num_repetitions
      xmin=min(xlim); xmax=max(xlim);
      ymin=min(ylim); ymax=max(ylim);
      text_xpos=xmin+.05*(xmax-xmin);
      text_ypos=ymin+.9*(ymax-ymin);
      text(text_xpos,text_ypos,['num\_realizations=' num2str(stats.num_repetitions)]);

      % plot Resonance Frequency with confidence intervals
      subplot(1,2,2);
      color='b';
      x=p;
      y=fcmu;
      e=fcsd/sqrt(N);
      plot_CI(x',y',e',color);
      xlabel(strrep(stats.(var).repetition_sets.varied{1},'_','\_'))
      ylabel('Resonance Frequency [Hz]');
      % add num_repetitions
      xmin=min(xlim); xmax=max(xlim);
      ymin=min(ylim); ymax=max(ylim);
      text_xpos=xmin+.05*(xmax-xmin);
      text_ypos=ymin+.9*(ymax-ymin);
      text(text_xpos,text_ypos,['num\_realizations=' num2str(stats.num_repetitions)]);
    end
  end
  % ---------------------------------------------------
  % Spectral Plots
  % ---------------------------------------------------
  colors='bgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmky';
  var=stats.MUA_variables{1};
  pow_freqs=stats.(var).frequency;
  Pxx_mu=stats.(var).repetition_sets.Pxx_mu;
  Pxx_sd=stats.(var).repetition_sets.Pxx_sd;
  fMUA_mu=stats.(var).repetition_sets.PeakFreq_mu;
  fMUA_sd=stats.(var).repetition_sets.PeakFreq_sd;
  p_values=stats.(var).repetition_sets.parameters;
  p_names=stats.(var).repetition_sets.varied;
  N=stats.num_repetitions;
  
  figure('position',[200 340 1500 550])
  subplot(1,2,1);
  for i=1:size(Pxx_mu,2)
    if length(p_names)==1
      color=colors(i);
    else
      try
        color=colors(p_values(i,1)==p);
      catch
        color=colors(i);
      end
    end
    x=pow_freqs;
    y=Pxx_mu(:,i);
    e=Pxx_sd(:,i)/sqrt(N);
    plot_CI(x,y,e,color);
    hold on
  end
  xlim([0 100]);
  xlabel('Spectral Frequency [Hz]')
  ylabel(strrep(var,'_','\_'));    
  % generate legend
  legend_str={};
  for i=1:size(p_values,2)
    str='';
    for j=1:size(p_values,1)
      str=[str sprintf('%s=%g,',p_names{j},p_values(j,i))];
    end
    legend_str{end+1}=strrep(str(1:end-1),'_','\_');
  end
  legend(legend_str);
  
  % plot Spectral Peak Frequency with confidence intervals
  subplot(1,2,2);
  color='b';
  if length(p_names)==1
    x=p_values;
  else
    x=(1:length(fMUA_mu))';
  end
  y=fMUA_mu;
  e=fMUA_sd/sqrt(N);
  plot_CI(x,y,e,color);
  xlabel(strrep([stats.(var).repetition_sets.varied{:}],'_','\_'))
  ylabel('Peak Spectral Frequency [Hz]');
  % add num_repetitions
  xmin=min(xlim); xmax=max(xlim);
  ymin=min(ylim); ymax=max(ylim);
  text_xpos=xmin+.05*(xmax-xmin);
  text_ypos=ymin+.9*(ymax-ymin);
  text(text_xpos,text_ypos,['num\_realizations=' num2str(stats.num_repetitions)]);

end


% %% Post-processing
% stats=[];
% % Calculate statistics across frequencies and realizations
% if nargout>1 || options.plot_flag==1
%   fprintf('Postprocessing: calculating resonance and spectral properties...\n');
%   % First collect information on what was varied
%   varied=data(1).varied;
%   num_varied=length(varied); % number of model components varied across simulations
%   num_sims=length(data); % number of data sets (one per simulation)
%   params_all=zeros(num_sims,num_varied); % values for each simulation
%   % loop over varied components and collect values
%   for i=1:num_varied
%     params_all(:,i)=[data.(varied{i})]; % values for each simulation
%       % each row of params_all contains the varied values for a single simulation.
%       % [num_sims x num_varied]
%   end
%   num_freqs=length(options.f);
%   % realization parameter
%   realization_param=[options.target{1} '_realization'];
%   if num_freqs>1
%     % firing rate analysis
%     data=CalcFR(data,'bin_size',options.bin_size,'bin_shift',options.bin_shift);
%     % frequency parameter
%     freq_param=[options.target{1} '_f'];
%     freq_idx=(ismember(varied,freq_param)); % index of freq param
%     other_idx=(~ismember(varied,freq_param)); % indices of all other varied params
%     % get varied params without freq parameter
%     params_sweeps=unique(params_all(:,other_idx),'rows','stable');
%       % each row of params_sweeps contains the varied values besides input
%       % frequency (each row has num_freqs corresponding sims).
%       % [params_sweeps] = (num_sims/num_freqs) x (num_varied-1)
%     num_sweeps=size(params_sweeps,1); % = (num_sims/num_freqs)
%     % get indices to sims for each sweep
%     [~,LOCB]=ismember(params_all(:,other_idx),params_sweeps,'rows');
%       % LOCB is an array indicating which sets of params varied match the
%       % subsets in params_sweeps
%     sweep_indices=nan(num_sweeps,num_freqs);
%     for i=1:num_sweeps
%       sweep_indices(i,:)=find(LOCB==i);
%         % each row of sweep_indices contains indices into data for one
%         % sweep across input frequencies.
%     end
%     if options.num_repetitions>1
%       % collect params for sweep sets (multiple realizations of a sweep)
%       varied_=varied(other_idx);
%       other_idx_sweeps=(~ismember(varied_,realization_param));
%       varied_sweeps=varied_(other_idx_sweeps);
%       params_sweeps_sets=unique(params_sweeps(:,other_idx_sweeps),'rows','stable');
%       num_sweeps_sets=size(params_sweeps_sets,1);
%       [~,LOCB]=ismember(params_sweeps(:,other_idx_sweeps),params_sweeps_sets,'rows');
%       sweep_set_indices=nan(num_sweeps_sets,options.num_repetitions);
%       for i=1:num_sweeps_sets
%         sweep_set_indices(i,:)=find(LOCB==i);
%       end
%     end
%     % get list of variables whose activity should be used to measure resonance
%     FR_fields=data(1).results(~cellfun(@isempty,regexp(data(1).results,'.*_FR$')));
%     FR_fields=setdiff(FR_fields,'time_FR');
%     stats.FR_variables=FR_fields;
%     % Calculate resonance frequency over each frequency sweep
%     for v=1:length(FR_fields)
%       var=FR_fields{v};
%       FR=zeros(num_sweeps,num_freqs);
%       fc=zeros(num_sweeps,1);
%       for i=1:num_sweeps
%         dat=data(sweep_indices(i,:));
%         FR(i,:)=cellfun(@(x)mean(x(:)),{dat.(var)});
%         fc(i,1)=options.f(FR(i,:)==max(FR(i,:)));
%       end
%       stats.(var).pop_FR=FR;                % [num_sweeps x num_freqs] where num_sweeps=(num_sims/num_freqs)
%       stats.(var).ResonanceFreq=fc;         % [num_sweeps x 1]
%       stats.(var).parameters=params_sweeps; % [num_sweeps x (num_varied-1)]
%       stats.(var).varied=varied(other_idx); % [1 x (num_varied-1)]
%       if options.num_repetitions>1        
%         % calculate mean/std resonance frequency
%         FRmu=zeros(num_sweeps_sets,num_freqs);
%         FRsd=zeros(num_sweeps_sets,num_freqs);
%         fcmu=zeros(num_sweeps_sets,1);
%         fcsd=zeros(num_sweeps_sets,1);
%         for i=1:num_sweeps_sets
%           inds=sweep_set_indices(i,:);
%           FRmu(i,:)=mean(FR(inds,:),1);
%           FRsd(i,:)=std(FR(inds,:),[],1);
%           fcmu(i,:)=mean(fc(inds,1),1);
%           fcsd(i,:)=std(fc(inds,1),[],1);          
%         end
%         stats.(var).repetition_sets.pop_FR_mu=FRmu;                % [num_sweep_sets x num_freqs], where num_sweep_sets=(num_sweeps/num_repetitions)
%         stats.(var).repetition_sets.pop_FR_sd=FRsd;                % [num_sweep_sets x num_freqs]
%         stats.(var).repetition_sets.ResonanceFreq_mu=fcmu;         % [num_sweep_sets x 1]
%         stats.(var).repetition_sets.ResonanceFreq_sd=fcsd;         % [num_sweep_sets x 1]
%         stats.(var).repetition_sets.parameters=params_sweeps_sets; % [num_sweep_sets x (num_varied-2)]
%         stats.(var).repetition_sets.varied=varied_sweeps;    % [1 x (num_varied-2)]
%       else
%         stats.(var).repetition_sets.pop_FR_mu=FR;               
%         stats.(var).repetition_sets.pop_FR_sd=zeros(size(FR));          
%         stats.(var).repetition_sets.ResonanceFreq_mu=fc;        
%         stats.(var).repetition_sets.ResonanceFreq_sd=zeros(size(fc));        
%         stats.(var).repetition_sets.parameters=params_sweeps;
%         stats.(var).repetition_sets.varied=varied(other_idx);            
%       end
%     end
%   end
%   stats.parameters=params_all;            % [num_sims x num_varied]
%   stats.varied=varied;                    % [1 x num_varied]
%   stats.frequencies=options.f;            % [1 x num_freqs]
%   stats.num_repetitions=options.num_repetitions;
%   stats.options=options;
%   
%   % spectral analysis (fMUA, fSUA; area power)
%   data=CalcPower(data);    
%   if options.num_repetitions>1 % Calculate mean/std over realizations
%     % Collect information on what was varied besides 'realization'
%     % Collect params for realization sets (multiple realizations of a point in search space)
%     other_idx_realizations=(~ismember(varied,realization_param));
%     varied_realizations=varied(other_idx_realizations);
%     params_realizations=unique(params_all(:,other_idx_realizations),'rows','stable');
%     num_repetition_sets=size(params_realizations,1);
%     [~,LOCB]=ismember(params_all(:,other_idx_realizations),params_realizations,'rows');
%     realization_set_indices=nan(num_repetition_sets,options.num_repetitions);
%     for i=1:num_repetition_sets
%       realization_set_indices(i,:)=find(LOCB==i);
%     end
%     % get list of variables whose power should be averaged over realizations
%     MUA_fields=data(1).results(~cellfun(@isempty,regexp(data(1).results,'.*_MUA$')));
%     SUA_fields=data(1).results(~cellfun(@isempty,regexp(data(1).results,'.*_SUA$')));
%     stats.MUA_variables=MUA_fields;
%     stats.SUA_variables=SUA_fields;
%     % calculate mean/std fMUA and fSUA
%     for v=1:length(MUA_fields)
%       var=MUA_fields{v};
%       pow_freqs=data(1).(var).frequency;
%       num_pow_freqs=length(pow_freqs);
%       PxxMUAmu=zeros(num_repetition_sets,num_pow_freqs);
%       PxxMUAsd=zeros(num_repetition_sets,num_pow_freqs);
%       fMUAmu=zeros(num_repetition_sets,1);
%       fMUAsd=zeros(num_repetition_sets,1);
%       paMUAmu=zeros(num_repetition_sets,1);
%       paMUAsd=zeros(num_repetition_sets,1);
%       for i=1:num_repetition_sets
%         dat=data(realization_set_indices(i,:));
%         fMUA=arrayfun(@(x)x.E_v_Power_MUA.PeakFreq,dat); % [1 x num_repetitions]
%         fMUAmu(i)=mean(fMUA,2);
%         fMUAsd(i)=std(fMUA,[],2);
%         paMUA=arrayfun(@(x)x.E_v_Power_MUA.PeakArea,dat);
%         paMUAmu(i)=mean(paMUA,2);
%         paMUAsd(i)=std(paMUA,[],2);
%         PxxMUA=arrayfun(@(x)x.E_v_Power_MUA.Pxx,dat,'uni',0);
%         PxxMUA=cat(2,PxxMUA{:}); % [num_pow_freqs x num_repetitions]
%         PxxMUAmu(i,:)=mean(PxxMUA,2);
%         PxxMUAsd(i,:)=std(PxxMUA,[],2);
%       end
%       PxxMUAmu(isnan(PxxMUAmu))=0; PxxMUAsd(isnan(PxxMUAsd))=0;
%       fMUAmu(isnan(fMUAmu))=0;     fMUAsd(isnan(fMUAsd))=0;
%       paMUAmu(isnan(paMUAmu))=0;   paMUAsd(isnan(paMUAsd))=0;
%       stats.(var).repetition_sets.Pxx_mu=PxxMUAmu;     % num_repetition_sets x num_pow_freqs, where num_repetition_sets=(num_sims/num_repetitions)
%       stats.(var).repetition_sets.Pxx_sd=PxxMUAsd;     % num_repetition_sets x num_pow_freqs
%       stats.(var).repetition_sets.PeakFreq_mu=fMUAmu;  % num_repetition_sets x 1
%       stats.(var).repetition_sets.PeakFreq_sd=fMUAsd;  % num_repetition_sets x 1
%       stats.(var).repetition_sets.PeakArea_mu=paMUAmu;  % num_repetition_sets x 1
%       stats.(var).repetition_sets.PeakArea_sd=paMUAsd;  % num_repetition_sets x 1
%       stats.(var).frequency=pow_freqs';                 % 1 x num_pow_freqs
%       stats.(var).parameters=params_realizations;       % num_repetition_sets x (num_varied-1), excluding 'realization' from varied
%       stats.(var).varied=varied_realizations;           % 1 x (num_varied-1)
%     end
%   end
% end
