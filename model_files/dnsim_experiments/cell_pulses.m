function [allX,stats] = cell_pulses(spec,varargin)
% purpose: characterize subthreshold response of single cells

parms = mmil_args2parms( varargin, ...
       {  ...
          'minamp',-7,[],...
          'maxamp',7,[],...
          'damp',1,[],...
          'pulsetime',200,[],...
          'baseline',100,[],...
          'verbose',0,[],...
          'coder',1,[],...
          'dt',.01,[],...
          'solver','euler',[],...
          'dsfact',1,[],...
          'var','V',[],...
          'placeholder','current',[],...
          'plot_flag',0,[],...
          'override',[],[],...
       }, false);
     
% parms.minamp=-7; 
% parms.maxamp=7; 
% parms.damp=1;
% parms.pulsetime=200;
% parms.baseline=100;
% parms.coder=1;
% parms.dt=.01;
% parms.solver='euler';
% parms.dsfact=1;
% parms.verbose=0;
% parms.var='V';
% parms.placeholder='current';
% parms.plot_flag=1;
% parms.override=[];

ncells=length(spec.nodes);

% stimuli
amps=parms.minamp:parms.damp:parms.maxamp;
amps=unique([-amps amps]);
namps=length(amps);

% time
t_on=parms.baseline;
t_off=parms.baseline+parms.pulsetime;
t_end=t_off+parms.pulsetime;
t=0:parms.dt:t_end;
t=t(1:parms.dsfact:end);
ntimes=length(t);

% preallocation
allX=zeros(ntimes,namps,ncells);

%% step inputs
tstart=tic;
for cellind=1:ncells % this cell type (node)
  fprintf('processing cell type %g of %g\n',cellind,ncells);

  var=[spec.nodes(cellind).label '_' parms.var]; % state variable to process

  for stepnum=1:namps % this pulse

    % prepare input and add to dynamics
    amp=amps(stepnum); % this pulse amplitude
    inpstring=sprintf('%s+(tmp_stim*(t>tmp_onset & t<tmp_offset))',parms.placeholder);
    dynamics=strrep(spec.nodes(cellind).dynamics,parms.placeholder,inpstring);

    % prepare model for this run
    model=[];
    model.nodes.label=spec.nodes(cellind).label;
    model.nodes.multiplicity=1;
    model.nodes.dynamics=dynamics;
    model.nodes.mechanisms=spec.nodes(cellind).mechanisms;
    model.nodes.parameters=spec.nodes(cellind).parameters;
    model.nodes.connections=[];

    override={model.nodes.label,'parameters','tmp_onset',t_on;
              model.nodes.label,'parameters','tmp_offset',t_off;
              model.nodes.label,'parameters','tmp_stim',amp};
    if ~isempty(parms.override)
      override=cat(1,override,parms.override);
    end
    
    [data,model] = runsim(model,'timelimits',[0 t_end],'override',override,'timesurfer_flag',0,'dt',parms.dt,'dsfact',parms.dsfact,'SOLVER',parms.solver,'coder',parms.coder,'verbose',parms.verbose);
    X=data.(var);
    
    allX(:,stepnum,cellind)=X;
    
  end % end pulse loop
end % end cell type loop
toc(tstart)

if parms.plot_flag
  figure('position',[140 60 1700 850])
  for cellind=1:ncells
    subplot(ncells,3,1+3*(cellind-1));
    plot(t,allX(:,amps<=0,cellind)); axis tight; xlabel('time (ms)'); ylabel('response amplitude');
    title(spec.nodes(cellind).label);
    subplot(ncells,3,2+3*(cellind-1));
    plot(t,allX(:,amps>=0,cellind)); axis tight; xlabel('time (ms)'); ylabel('response amplitude');
    title(spec.nodes(cellind).label);
    subplot(ncells,3,3+3*(cellind-1));
    imagesc(t,amps,allX(:,:,cellind)'); axis xy; colorbar; colormap(1-gray)
    xlabel('time (ms)'); ylabel('input amplitude');  title(spec.nodes(cellind).label);
  end
end

% process results
% fprintf('\nProcessing results...\n');
threshold=0;
tind=find(t>=t_on & t<t_off);
Veq=nan(namps,ncells);
Rss=nan(namps,ncells);

ISI1=nan(namps,ncells);
ISIss=nan(namps,ncells);
ISIadapt=nan(namps,ncells);
PowAreaSpiking=nan(namps,ncells);
PowAreaSubthresh=nan(namps,ncells);
PowFreqSpiking=nan(namps,ncells);
PowFreqSubthresh=nan(namps,ncells);

taum=nan(namps,ncells);
Reff=nan(1,ncells);
fIslope=nan(1,ncells);
for cellind=1:ncells 
  for stepnum=1:namps
    X=allX(:,stepnum,cellind);
    % check for spikes in [t_on .9(t_off-t_on)]
    if any(X(tind)>threshold)
      [~,ind] = findpeaks(double(X(tind)),'MinPeakHeight',threshold);
      if length(ind)>1
        % calc from mean ISIs
        Rss(stepnum,cellind)=1/(mean(diff(t(tind(ind))))/1000);
      else
%         Rss(stepnum,cellind)=length(ind)/(parms.pulsetime/1000);
      end  
      % if spike in 50-100%: calc adaptation ratio (rSS/r0)
      % ...
    else
      % if no spike: calc Veq (over 50-100%)
      Veq(stepnum,cellind)=mean(X(t>=(t_on+.5*(t_off-t_on)) & t<t_off));
      if amps(stepnum)~=0
        % calc effective membrane time constant (time to Vmax(1/e))
        bl=X(tind(1)-1); % baseline value is point immediately before stim onset
        V=X-bl; % shift baseline to V=0
        if amps(stepnum)<0
          V=abs(V); % invert values so that deflection is positive
        end
        Vmax=V(tind(end)); % max is final point during stim period before offset
        Vthresh=Vmax*(1/exp(1)); % voltage at 1/e
        Vdecay=V(tind(end):end); % voltage decay after stim offset
        tthresh=t(tind(end)+find(Vdecay<Vthresh,1,'first')-1); % time at which decay crosses 1/e
        taum(stepnum,cellind)=tthresh-t_off; % time to drop to 1/e after stim offset
        %figure; plot(t,V); hline(Vthresh,'k'); vline(tthresh,'k');
      end
    end
  end
  % calc Reff (from slope of amp/Veq)
  if any(~isnan(Veq(:,cellind)))
    sel=~isnan(Veq(:,cellind));
    P=polyfit(amps(sel),Veq(sel,cellind)',1);
    Reff(cellind)=P(1);
  end
  if any(~isnan(Rss(:,cellind)))
    sel=~isnan(Rss(:,cellind));
    P=polyfit(amps(sel),Rss(sel,cellind)',1);
    fIslope(cellind)=P(1);
  end
end

stats.amps=amps;
stats.Veq=Veq;
stats.Rss=Rss;
stats.Reff=Reff;
stats.fIslope=fIslope;
stats.Vthresh=nanmax(stats.Veq,[],1);
stats.taum=taum;
stats.taum_median=nanmedian(taum,1);
stats.model=model;
stats.time=t;


