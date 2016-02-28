function [model,vary]=PrepProbeTwoRhythms(model,varargin)
% Experiment: data=ProbeTwoRhythms(model,'option',value,...)
% Deliver two rhythmic inputs, each to >=50% of the target population(s),
% with optional overlap, varying the frequency of the inputs over a grid
% (f1,f2). The number of spikes and total effective strength of the inputs
% are balanced across all cells. The input may be a direct current
% injection (sinusoid, noisy sinusoid, rectified sinusoid, or gaussian
% smoothed pulses) or nonhomogeneous Poisson process filtered through
% double-exponential AMPA synapses. The purpose is to assess how the
% spectral components of the inputs are reflected in the target network.
% Additionally, the onset phase of the two inputs can be varied, mimicking
% experimental stimulus onset asynchrony or upstream differences in 
% neural activity transmission delays.
% see also: ProbeResonance

% Check options
options=CheckOptions(varargin,{...
  'input_type','poisson',{'sin','noisy_sin','rectified_sin','poisson','noisy_rectified_sin'},...
  'target',[],[],... % target populations; default: first population only
  'SOA',0,[],... % stimulus onset asynchrony
  'gINPUT',1,[],... % max amplitude of applied current
  'f1',35,[],... % Hz, modulation frequency of first stimulus
  'f2',45,[],... % Hz, modulation frequency of second stimulus
  'num_repetitions',1,[],... % number of repetitions of the input (set to >1 for noisy inputs if stats are desired)
  'dc',0,[],... % kHz
  'ac',1,[],... % kHz
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

% convert SOA to phase1 and phase2
% ...

% Add two inputs to model (each to 50% of target populations) and set 'vary' (f1,f2,repetition)
modifications={};
vary={};
for i=1:npops
  name=model.specification.populations(i).name;
  N=model.specification.populations(i).size;
  if ~ismember(name,options.target)
    % add nothing
    continue;
  end
  % define kernels
  % input to first 50% of cells
  K1=zeros(1,N);
  K1(1:ceil(N/2))=1;
  % input to second 50% of cells
  K2=zeros(1,N);
  K2(ceil(N/2)+1:N)=1;
  % prepare list of modifications to add tonic drive to all populations in model
  switch options.input_type
    case 'poisson'
      % add AMPA synapse for exponentially filtered modulated poisson process
      s1=sprintf('s1=get_input(''%s'',Npop,T,f1,dc,ac,tau,K1,baseline,phase1); phase1=0',options.input_type);
      s2=sprintf('s2=get_input(''%s'',Npop,T,f2,dc,ac,tau,K2,baseline,phase2); phase2=0',options.input_type);
      eqn_mods=sprintf('cat(ODE1,-gINPUT*(s1(k,:)+s2(k,:)).*(X-0); %s; %s)',s1,s2);
      modifications(end+1,:)={name,'equations',eqn_mods};
      modifications(end+1,:)={name,'K1',K1};
      modifications(end+1,:)={name,'K2',K2};
      param_names={'gINPUT','dc','ac','tau','baseline','f1','f2','repetition'};
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
      % add +sin() to equations and (f1,f2) to parameters
      I1=sprintf('I1=get_input(''%s'',Npop,T,f1)',options.input_type);
      I2=sprintf('I2=get_input(''%s'',Npop,T,f2)',options.input_type);
      eqn_mods=sprintf('cat(ODE1,+gINPUT*(K1.*I1(k,:)+K2.*I2(k,:)); %s; %s)',I1,I2);
      modifications(mod_cnt+0,:)={name,'equations',eqn_mods};
      modifications(mod_cnt+1,:)={name,'gINPUT',options.gINPUT};
      modifications(mod_cnt+2,:)={name,'K1',K1};
      modifications(mod_cnt+3,:)={name,'K2',K2};
      % prepare specification to vary tonic drive across simulations
      vary(vary_cnt+0,:)={name,'f1',options.f1};
      vary(vary_cnt+1,:)={name,'f2',options.f2};
      vary(vary_cnt+2,:)={name,'repetition',1:options.num_repetitions};
      mod_cnt=mod_cnt+4;
      vary_cnt=vary_cnt+3;
    otherwise
      error('input_type not recognized.');
  end
end

% apply modifications to effectively add experimental apparatus to model
model=ApplyModifications(model,modifications);

% fprintf('Running experiment: %s\n',mfilename);
% 
% % execute experimental protocol by varying parameters (eg, f1,f2) across simulations
% data=SimulateModel(model,'vary',vary,varargin{:});
% 
% %% Post-processing: calculate fMUA and <FR> over repetitions
% if nargout>1
%   fprintf('Post-processing: analyzing response properties...\n');
%   N=model.specification.populations(1).size;
%   variable=data(1).labels{1};
%   dat1=SelectData(data,'roi',{variable,[1 ceil(N/2)]});
%   dat2=SelectData(data,'roi',{variable,[ceil(N/2)+1 N]});
%   repetition_parameter=[options.target{1} '_repetition'];
%   stats(1)=CalcResonanceStats(dat1,'variable',variable,'repetition_parameter',repetition_parameter,varargin{:});
%   stats(2)=CalcResonanceStats(dat2,'variable',variable,'repetition_parameter',repetition_parameter,varargin{:});  
%   clear dat1 dat2
% end
% 
% % post-process downsampling before returning data
% dsfact=options.post_downsample_factor;
% if dsfact>1
%   for i=1:length(data)
%     for j=1:length(data(i).labels)
%       data(i).(data(i).labels{j})=data(i).(data(i).labels{j})(1:dsfact:end,:);
%       data(i).(data(i).labels{j})=single(data(i).(data(i).labels{j}));
%     end
%   end
% end
