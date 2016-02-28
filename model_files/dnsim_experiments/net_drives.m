function [allX,stats] = net_drives(spec,varargin)

parms = mmil_args2parms( varargin, ...
       {  ...
          'minamp',-7,[],...
          'maxamp',7,[],...
          'damp',1,[],...
          'targets',[],[],...
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
  
% parms.minamp=-3; 
% parms.maxamp=3; 
% parms.damp=1;
% parms.targets=[];
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

if isempty(parms.targets)
  parms.targets={spec.nodes(1).label};
end
if ~iscell(parms.targets)
  parms.targets={parms.targets};
end

ncells=length(spec.nodes);
ntargets=length(parms.targets);

% stimuli
amps=parms.minamp:parms.damp:parms.maxamp;
% amps=unique([-amps amps]);
namps=length(amps);
% expand to specify amp for each desired cell type for all sims
allamps=[];
allamps(:,1)=amps;
for i=2:ntargets
  allamps=repmat(allamps,[ones(1,i-1) namps]);
end
nsims=numel(allamps);

% prepare model specification (add input for sweeping drives)
model=spec;
inpstring=sprintf('%s+(tmp_stim*(t>tmp_onset & t<tmp_offset))',parms.placeholder);
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
allX=nan([ntimes,size(allamps),ncells]);
Ravg=nan([size(allamps) ncells]);

%% vary drives
tstart=tic;
for sim=1:nsims
  fprintf('processing sim %g of %g\n',sim,nsims);
  % prepare parameters for this simulation
  clear inds
  [inds{1:ntargets}]=ind2sub(size(allamps),sim);
  override={};
  for target=1:ntargets
    amp=amps(inds{target});
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_onset',t_on};
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_offset',t_off};
    [override(end+1,1:4)]={parms.targets{target},'parameters','tmp_stim',amp};
  end  
  if ~isempty(parms.override)
    override=cat(1,override,parms.override);
  end  
  data=runsim(model,'timelimits',[0 t_end],'override',override,'timesurfer_flag',0,'dt',parms.dt,'dsfact',parms.dsfact,'SOLVER',parms.solver,'coder',parms.coder,'verbose',parms.verbose);
  
  % process results
  for cellind=1:ncells
    var=[model.nodes(cellind).label '_' parms.var]; % state variable to process
    allX(:,inds{:},cellind)=nanmean(data.(var),2);
    n=model.nodes(cellind).multiplicity;
    for thiscell=1:n
      X=data.(var)(:,thiscell);
      rsum=0;
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

% Process MUA results
fprintf('\nProcessing results...');
Fmin=7;
Fmax=500;
PxxSmooth=1;
powthreshprc=95;
Fwin=5;
OscFreq=nan([size(allamps),ncells]);
AreaPower=nan([size(allamps),ncells]);
for sim=1:nsims
  clear inds
  [inds{1:ntargets}]=ind2sub(size(allamps),sim);
  for cellind=1:ncells
    X=allX(:,inds{:},cellind);
    [P,f,Pstats]=js_pwelch(X(tind),t(tind)/1000,Fmin,Fmax,PxxSmooth,powthreshprc,Fwin);
    if sim==1 && cellind==1
      allP=nan([length(f),size(allamps),ncells]);
    end
    allP(:,inds{:},cellind)=P;
    OscFreq(inds{:},cellind)=Pstats.OscFreq;
    AreaPower(inds{:},cellind)=Pstats.AreaPower;
    if parms.plot_flag, figure; subplot(2,1,1); plot(t,X); axis tight; subplot(2,1,2); plot(f,P,'.-'); xlim([0 100]); vline(Pstats.OscFreq,'k'); title(sprintf('fMUA=%g, area=%g',Pstats.OscFreq,Pstats.AreaPower)); pause(.5); close; end
  end
end
fprintf('done.\n');

if parms.plot_flag
  figure('position',[140 60 1000 850]); nc=3;
  for cellind=1:ncells
    if ntargets==2
      subplot(ncells,nc,1+nc*(cellind-1));
      imagesc(amps,amps,Ravg(:,:,cellind)); axis tight; axis square; axis xy
      title(sprintf('<r%s>',spec.nodes(cellind).label));
      xlabel([parms.targets{2} ' input amplitude']); 
      ylabel([parms.targets{1} ' input amplitude']);
      colorbar; caxis([min(Ravg(:)) max(Ravg(:))]);
      colormap(1-gray); vline(0,'b'); hline(0,'b');
      subplot(ncells,nc,2+nc*(cellind-1));
      imagesc(amps,amps,OscFreq(:,:,cellind)); axis tight; axis square; axis xy
      title(sprintf('<fMUA\\_%s>',spec.nodes(cellind).label));
      xlabel([parms.targets{2} ' input amplitude']); 
      ylabel([parms.targets{1} ' input amplitude']);
      colorbar; caxis([min(OscFreq(:)) max(OscFreq(:))]);
      colormap(1-gray); vline(0,'b'); hline(0,'b');
      subplot(ncells,nc,3+nc*(cellind-1));
      imagesc(amps,amps,AreaPower(:,:,cellind)); axis tight; axis square; axis xy
      title(sprintf('<AreaPower\\_%s>',spec.nodes(cellind).label));
      xlabel([parms.targets{2} ' input amplitude']); 
      ylabel([parms.targets{1} ' input amplitude']);
      colorbar; caxis([min(AreaPower(:)) max(AreaPower(:))]);
      colormap(1-gray); vline(0,'b'); hline(0,'b');
    end
  end
end
% figure; plot(Ravg(:,:,1),OscFreq(:,:,1))

stats.amps=amps;
stats.Ravg=Ravg;
stats.time=t;
stats.freq=f;
stats.allPower=allP;
stats.OscFreq=OscFreq; % target drives x celltypes
stats.AreaPower=AreaPower;


