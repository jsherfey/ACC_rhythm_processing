function [allX,stats] = net_tuning(spec,varargin)

parms = mmil_args2parms( varargin, ...
       {  ...
          'minamp',7,[],...
          'maxamp',7,[],...
          'damp',1,[],...
          'minfreq',2,[],...
          'maxfreq',60,[],...
          'dfreq',2,[],...
          'targets',[],[],...
          'pulsetime',1000,[],...
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
  
% parms.targets=[];
% parms.minamp=3; % 2
% parms.maxamp=3; % 2
% parms.damp=1;
% parms.minfreq=5; % 5
% parms.maxfreq=60;% 50
% parms.dfreq=5; % 5
% parms.pulsetime=1000;
% parms.baseline=100;
% parms.verbose=0;
% parms.coder=1;
% parms.dt=.01;
% parms.solver='euler';
% parms.dsfact=10;
% parms.var='V';
% parms.placeholder='current';
% parms.plot_flag=1;
% parms.override=[];

if isempty(parms.targets)
  parms.targets={spec.nodes(1).label};
end
if ~iscell(parms.targets)
  parms.targets={parms.targets};
end

ncelltypes=length(spec.nodes);
ntargets=length(parms.targets);

% stimulus parameters
amps=parms.minamp:parms.damp:parms.maxamp;
namps=length(amps);
freqs=parms.minfreq:parms.dfreq:parms.maxfreq;
nfreqs=length(freqs);

% expand to specify amp and freq for each desired cell type for all sims
allfreqs=repmat(freqs,[namps 1]);
allfreqs=allfreqs(:);
allamps=repmat(amps,[1 nfreqs])';
nvariations=length(allamps);
for i=2:ntargets
  allamps=repmat(allamps,[ones(1,i-1) nvariations]);
  allfreqs=repmat(allfreqs,[ones(1,i-1) nvariations]);
end
nsims=numel(allamps);

% prepare model specification (add input for sweeping drives)
model=spec;

inpstring=sprintf('%s+(tmp_stim*(t>tmp_onset & t<tmp_offset)*((1+sin(2*pi*tmp_freq*t/1000))/2))',parms.placeholder);
% inpstring=sprintf('%s+(tmp_stim*(t>tmp_onset & t<tmp_offset))',parms.placeholder);
for target=1:ntargets
  cellind=find(strcmp(parms.targets{target},{spec.nodes.label}));
  dynamics=strrep(spec.nodes(cellind).dynamics,parms.placeholder,inpstring);
  model.nodes(cellind).dynamics=dynamics;
end

% time
t_on=parms.baseline;
t_off=parms.baseline+parms.pulsetime;
t_end=t_off;%+parms.pulsetime;
t=0:parms.dt:t_end;
t=t(1:parms.dsfact:end);
ntimes=length(t);

tind=find(t>=t_on & t<t_off);
threshold=0;

% preallocation
allX=nan([ntimes,size(allamps),ncelltypes]);
Ravg=nan([size(allamps) ncelltypes]);

%% vary drives
tstart=tic;
for sim=1:nsims
  fprintf('processing sim %g of %g\n',sim,nsims);
  % prepare parameters for this simulation
  clear inds
  [inds{1:ntargets}]=ind2sub(size(allamps),sim);
  override={};
  for target=1:ntargets
%     amp=amps(inds{target});
%     freq=freqs(inds{target});
    amp=allamps(inds{target});
    freq=allfreqs(inds{target});
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_onset',t_on};
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_offset',t_off};
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_stim',amp};
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_freq',freq};
  end  
  if ~isempty(parms.override)
    override=cat(1,override,parms.override);
  end  
  data=runsim(model,'timelimits',[0 t_end],'override',override,'timesurfer_flag',0,'dt',parms.dt,'dsfact',parms.dsfact,'SOLVER',parms.solver,'coder',parms.coder,'verbose',parms.verbose);
  
  % process results
  for cellind=1:ncelltypes % cell types
    var=[model.nodes(cellind).label '_' parms.var]; % state variable to process
    allX(:,inds{:},cellind)=nanmean(data.(var),2);
    n=model.nodes(cellind).multiplicity;
    rsum=0;
    for thiscell=1:n
      X=data.(var)(:,thiscell);
      if any(X(tind)>threshold)
        [~,ind] = findpeaks(double(X(tind)),'MinPeakHeight',threshold);
        if length(ind)>1
          rsum=rsum+(1/(mean(diff(t(tind(ind))))/1000));
        end
      end
    end % cell loop
    Ravg(inds{:},cellind)=rsum/n;
  end % cell type loop
  
end % end loop over simulations
toc(tstart)

% Process results (compute filter center frequency fc and bandwidth BW)
center_freq=nan([namps,ncelltypes]);
BW=nan([namps,ncelltypes]);
Q=nan([namps,ncelltypes]);
slope=nan([namps,ncelltypes]);

% for sim=1:nsims
%   clear inds
%   [inds{1:ntargets}]=ind2sub(size(allamps),sim);
%   for cellind=1:ncelltypes
%     X=allX(:,inds{:},cellind);

for cellind=1:ncelltypes % cell types
  rtmp=reshape(Ravg(:,:,cellind),[namps nfreqs]);
  for ampind=1:namps
    x=freqs;
    y=rtmp(ampind,:);
    peakthresh=prctile(y,65);
    [ypeaks,ind]=findpeaks(y,'MinPeakHeight',peakthresh);
    if ~isempty(ind)
      ind=ind(ypeaks==max(ypeaks));
      ind=ind(1);
    else
      ind=find(y==max(y),1,'first');
    end
    fc=x(ind); a=y(ind);
    % find corner frequencies (at -3dB/50% level) and amplitudes
    cornerthresh=min(y)+.5*range(y);
    [ind,fcross,ycross]=crossing(y,x,cornerthresh);
    if length(fcross)>1
      iind1=(find(fcross<fc,1,'last'));  
      if isempty(iind1)
        % no crossings below center
        f1=x(1);
        b1=y(1);
      else
        ind1=ind(iind1);
        f1=fcross(iind1); 
        b1=y(ind1);
      end        
      iind2=(find(fcross>fc,1,'first')); 
      if isempty(iind2)
        % no crossings above center
        f2=x(end);
        b2=y(end);
      else
        ind2=ind(iind2);
        f2=fcross(iind2); 
        b2=y(ind2);      
      end
      % filter characteristics
      s1=(a-b1)/(fc-f1);
      s2=(a-b2)/(f2-fc);
      BW(ampind,cellind)=f2-f1;
      Q(ampind,cellind)=fc/(f2-f1);
      slope(ampind,cellind)=(s1+s2)/2;
      center_freq(ampind,cellind)=fc;
    end
  end
end

if parms.plot_flag
  figure('position',[140 60 1000 850]); nc=1;
  for cellind=1:ncelltypes % cell types
    subplot(ncelltypes,nc,1+nc*(cellind-1));
    rtmp=reshape(Ravg(:,:,cellind),[namps nfreqs]);
    plot(freqs,rtmp,'.-'); axis tight
    title(sprintf('<r%s>',spec.nodes(cellind).label));
    legend(cellfun(@(x)sprintf('%s.amp=%g',parms.targets{1},x),num2cell(amps),'uni',0),'location','eastoutside');
    for ampind=1:namps
      vline(center_freq(ampind,cellind),'k');
    end
  end
  drawnow;
end

stats.amps=amps;
stats.namps=namps;
stats.freqs=freqs;
stats.nfreqs=nfreqs;
stats.Ravg=Ravg;              % (amp*freq) x ..targets.. x celltype
stats.filter_fc=center_freq;  % amp x celltype
stats.filter_BW=BW;           % amp x celltype
stats.filter_Q=Q;             % amp x celltype
stats.filter_slope=slope;     % amp x celltype





