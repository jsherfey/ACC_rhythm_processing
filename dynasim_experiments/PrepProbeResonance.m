function [model,vary]=PrepProbeResonance(model,varargin)
% Experiment: data=ProbeResonance(model,'option',value,...)
% Deliver one rhythmic input, varying the frequency of the inputs over
% range in order to determine what frequency generates the maximal
% response. The input may be a direct current injection (sinusoid, noisy 
% sinusoid, rectified sinusoid, or gaussian smoothed pulses) or 
% nonhomogeneous Poisson process filtered through double-exponential AMPA synapses.
% see also: ProbeTwoRhythms, CalcResonanceStats

% Check options
options=CheckOptions(varargin,{...
  'input_type','poisson',{'sin','noisy_sin','rectified_sin','poisson','noisy_rectified_sin'},...
  'target',[],[],... % target populations; default: first population only
  'gINPUT',1,[],... % max amplitude of applied current
  'f',5:5:60,[],... % Hz, modulation frequency (0: homogeneous poisson process)
  'num_repetitions',1,[],... % number of realizations of the input (set to >1 for noisy inputs if stats are desired)
  'DC',0,[],... % kHz
  'AC',1,[],... % kHz
  'baseline',.1,[],... % kHz
  'tau',2,[],... % ms
  },false);

options.repetition=1:options.num_repetitions;
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

% Add input to model and set 'vary' (f,repetition; other options)
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
      param_names={'gINPUT','DC','AC','tau','baseline','f','repetition'};
      ID=[];
      for p=1:length(param_names)
        if numel(options.(param_names{p}))>1
          ID=[ID p];
        end
      end
      % add AMPA synapse for exponentially filtered modulated poisson process
      s=sprintf('s=get_input(''%s'',Npop,T,f,DC,AC,tau,ones(1,Npop),baseline,phase); phase=0',options.input_type);
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
      vary(vary_cnt+1,:)={name,'repetition',1:options.num_repetitions};
      mod_cnt=mod_cnt+2;
      vary_cnt=vary_cnt+2;
    otherwise
      error('input_type not recognized.');
  end
end

% apply modifications to effectively add experimental apparatus to model
model=ApplyModifications(model,modifications);

% fprintf('Running experiment: %s\n',mfilename);
% 
% % execute experimental protocol by varying parameters across simulations
% data=SimulateModel(model,'vary',vary,varargin{:});
% 
% %% post-processing
% % Calculate statistics across frequencies and realizations
% if nargout>1 || options.plot_flag==1
%   fprintf('Postprocessing: calculating resonance and spectral properties...\n');
%   
%   realization_param=[options.target{1} '_realization'];
%   freq_param=[options.target{1} '_f'];
% 
%   [stats,data]=CalcResonanceStats(data,...
%     'sweep_parameter',freq_param,...       % max over f
%     'repetition_parameter',realization_param,... % mean/std over realizations
%     varargin{:});
% end
% 
% if options.plot_flag && options.num_repetitions>1
%   % ---------------------------------------------------
%   % Resonance Plots
%   % ---------------------------------------------------
%   if isfield(stats,'FR_variables') && numel(stats.(stats.FR_variables{1}).repetition_sets.varied)==1
%     var=stats.FR_variables{1};
%     if isfield(stats.(var),'sweep_sets')
%       N=stats.num_repetitions;
%       sweep_parameter=stats.(var).sweep_sets.sweep_parameter;
%       f=stats.(var).sweep_sets.(sweep_parameter);
%       p=stats.(var).sweep_sets.repetition_sets.parameters;
%       FRmu=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_mu;
%       FRsd=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_sd;
%       fcmu=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_mu;
%       fcsd=stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_sd;
% 
%       % plot Firing Rate Tuning Curves
%       colors='bgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmky';
%       figure('position',[200 340 1500 550])
%       subplot(1,2,1);
%       for i=1:length(p)
%         color=colors(i);
%         x=f;
%         y=FRmu(:,i);
%         e=FRsd(:,i)/sqrt(N);
%         plot_CI(x,y,e,color);
%         hold on
%       end
%       xlabel('Input Frequency')
%       ylabel([strrep(var,'_','\_') ' (population mean) [Hz]']);  
%       % add legend
%       pname=strrep(stats.(var).repetition_sets.varied{1},'_','\_');
%       legend(cellfun(@(x)[pname '=' num2str(x)],num2cell(p),'uni',0))
%       % add num_repetitions
%       xmin=min(xlim); xmax=max(xlim);
%       ymin=min(ylim); ymax=max(ylim);
%       text_xpos=xmin+.05*(xmax-xmin);
%       text_ypos=ymin+.9*(ymax-ymin);
%       text(text_xpos,text_ypos,['num\_realizations=' num2str(stats.num_repetitions)]);
% 
%       % plot Resonance Frequency with confidence intervals
%       subplot(1,2,2);
%       color='b';
%       x=p;
%       y=fcmu;
%       e=fcsd/sqrt(N);
%       plot_CI(x',y',e',color);
%       xlabel(strrep(stats.(var).repetition_sets.varied{1},'_','\_'))
%       ylabel('Resonance Frequency [Hz]');
%       % add num_repetitions
%       xmin=min(xlim); xmax=max(xlim);
%       ymin=min(ylim); ymax=max(ylim);
%       text_xpos=xmin+.05*(xmax-xmin);
%       text_ypos=ymin+.9*(ymax-ymin);
%       text(text_xpos,text_ypos,['num\_realizations=' num2str(stats.num_repetitions)]);
%     end
%   end
%   % ---------------------------------------------------
%   % Spectral Plots
%   % ---------------------------------------------------
%   colors='bgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmkybgrcmky';
%   var=stats.MUA_variables{1};
%   pow_freqs=stats.(var).frequency;
%   Pxx_mu=stats.(var).repetition_sets.Pxx_mu;
%   Pxx_sd=stats.(var).repetition_sets.Pxx_sd;
%   fMUA_mu=stats.(var).repetition_sets.PeakFreq_mu;
%   fMUA_sd=stats.(var).repetition_sets.PeakFreq_sd;
%   p_values=stats.(var).repetition_sets.parameters;
%   p_names=stats.(var).repetition_sets.varied;
%   N=stats.num_repetitions;
%   
%   figure('position',[200 340 1500 550])
%   subplot(1,2,1);
%   for i=1:size(Pxx_mu,2)
%     if length(p_names)==1
%       color=colors(i);
%     else
%       try
%         color=colors(p_values(i,1)==p);
%       catch
%         color=colors(i);
%       end
%     end
%     x=pow_freqs;
%     y=Pxx_mu(:,i);
%     e=Pxx_sd(:,i)/sqrt(N);
%     plot_CI(x,y,e,color);
%     hold on
%   end
%   xlim([0 100]);
%   xlabel('Spectral Frequency [Hz]')
%   ylabel(strrep(var,'_','\_'));    
%   % generate legend
%   legend_str={};
%   for i=1:size(p_values,2)
%     str='';
%     for j=1:size(p_values,1)
%       str=[str sprintf('%s=%g,',p_names{j},p_values(j,i))];
%     end
%     legend_str{end+1}=strrep(str(1:end-1),'_','\_');
%   end
%   legend(legend_str);
%   
%   % plot Spectral Peak Frequency with confidence intervals
%   subplot(1,2,2);
%   color='b';
%   if length(p_names)==1
%     x=p_values;
%   else
%     x=(1:length(fMUA_mu))';
%   end
%   y=fMUA_mu;
%   e=fMUA_sd/sqrt(N);
%   plot_CI(x,y,e,color);
%   xlabel(strrep([stats.(var).repetition_sets.varied{:}],'_','\_'))
%   ylabel('Peak Spectral Frequency [Hz]');
%   % add num_repetitions
%   xmin=min(xlim); xmax=max(xlim);
%   ymin=min(ylim); ymax=max(ylim);
%   text_xpos=xmin+.05*(xmax-xmin);
%   text_ypos=ymin+.9*(ymax-ymin);
%   text(text_xpos,text_ypos,['num\_realizations=' num2str(stats.num_repetitions)]);
% 
% end
% 

