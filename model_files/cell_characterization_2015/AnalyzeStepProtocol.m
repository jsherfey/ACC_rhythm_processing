% Purpose: characterize simulated cell response to step inputs (depol, hyperol) and tonic thresh input
% inputs: data, spec

%{
section = step of a given size
step = pulse of a given section
%}
tic

plot_flag=1;
if ~exist('datatype','var')
  datatype='sim';
end

switch datatype
  case 'sim'
    model=buildmodel(spec,'verbose',0);
    protocol='iStepProtocol';
    varind=find(strcmp([spec.nodes(1).label '_v'],{data.sensor_info.label}));
    % extract data
    t=data.epochs.time;
    V=data.epochs.data(varind,:);
    i=strcmp(protocol,model.nodes(1).mechanisms);
    p=model.nodes(1).mechs(i).params;
    dt=t(2)-t(1);
    Fs=round(1/dt);
    % get stimulation protocol info
    prestep_time_offset=0;   % shift sections and artifacts so onset artifact can be detected by tallie's cell characterization scripts
    P=getCharParms(model,'sim','cellind',1,'mechanism',protocol,'prestep_time_offset',prestep_time_offset,'time',t,'cellfld','nodes');
    % add artifact to intracellular voltage
    if 0
      sz=200/1000;                % artifact size, [V]
      prestep_time_offset=-.01;   % shift sections and artifacts so onset artifact can be detected by tallie's cell characterization scripts
      P=getCharParms(model,'sim','cellind',1,'mechanism',protocol,'prestep_time_offset',prestep_time_offset,'time',t,'cellfld','nodes');
      funcs = fieldnames(P);
      field=mean(V,3)/1000; % pop average
      intra=V(:,:,1)/1000;  % 1 cell
      for f=1:2 % make hyperpol artifacts up=>down; depol down=>up
        tstart=P.(funcs{f}).sections_start_sec - prestep_time_offset;
        tinf=tstart+((0:p.nsteps*p.nsections-1)*p.isi)/1000;
        tind=round(tinf/dt);
        intra(tind)=sz*((-1)^(f-1)); % onset artifact
        tinf=tinf+p.steptime/1000;
        tind=round(tinf/dt);
        intra(tind)=-sz*((-1)^(f-1)); % offset artifact
        %figure('position',[25 420 1550 250]); plot(t,intra);
      end
    end
    Y=double(V);
    nsect=p.nsections;
    nstep=p.nsteps;
    stepL=p.steptime/1000;
    stepI=p.stepsize;
  case 'exp'
    V=e1*1000;
    dt=t(2)-t(1);
    Fs=round(1/dt);
    prestep_time_offset=0;   % shift sections and artifacts so onset artifact can be detected by tallie's cell characterization scripts
    P = getCharParms(depinfo,'exp');
    tmp = getCharParms(hypinfo,'exp');
    P.CharHyperpolStepTA = tmp.CharHyperpolStepTA;
    tmp = getCharParms(toninfo,'exp');
    P.CharDepolTonicSpikesTA = tmp.CharDepolTonicSpikesTA;
    t0=P.CharDepolStepTA.sections_start_sec;
    tL=P.CharDepolStepTA.sections_length_sec;
    if t0==P.CharDepolTonicSpikesTA.sections_start_sec
      P.CharDepolTonicSpikesTA.sections_start_sec = t0+tL*length(P.CharDepolStepTA.sections_label_num);
      P.CharDepolTonicSpikesTA.sections_length_sec = floor(t(end)-P.CharDepolTonicSpikesTA.sections_start_sec);
    end
    Y=double(V);
    stepL=steptime;  
    stepI=P.CharDepolStepTA.sections_label_num(1);
end

if plot_flag
  figure('position',[150 100 1550 850]); 
  subplot(4,1,1); plot(t,V); axis tight; title([datatype ' data']);
end

%% Hyperpolarizing steps

%{
OUTPUTS:
o1.Baseline_mV
o1.step_sections
o1.sections_label_num
o1.Resistance_Mohms

INPUTS:
% {t,intra}
%             offset_voltage: 0
%     tonic_injected_current: 0
%         sections_label_num: [1 2 3 4]
%         sections_start_sec: 0.4900
%        sections_length_sec: 1
%         baseline_start_sec: 0.2500
%        baseline_length_sec: 0.2500

(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec,plot_flag)

%}

switch datatype
  case 'sim'
    FinPad1=200;
    FinPad2=70;    
  case 'exp'
    FinPad1=500;
    FinPad2=250;    
end

if isfield(P,'CharHyperpolStepTA') && nsect>0
  t0=P.CharHyperpolStepTA.sections_start_sec;
  L =P.CharHyperpolStepTA.sections_length_sec;
  bt0=P.CharHyperpolStepTA.baseline_start_sec;
  btL=P.CharHyperpolStepTA.baseline_length_sec;
  baseline=max(1,floor(Fs*bt0))+(0:floor(Fs*btL)-1);
  sections_start_sec = t0+(0:nsect-1)*L;
  sections_start_ind = floor(Fs*sections_start_sec);
  %mloc_off=sections_start_ind+floor(Fs*stepL);
  npts = floor(Fs*L);
  if isnan(baseline), baseline=2; end

  Baseline_mV = mean(Y(baseline));
  sections_label_num = P.CharHyperpolStepTA.sections_label_num;

  sections = nan(npts,nsect);
  End_mV = nan(1,nsect);
  Resistance_Mohms = nan(1,nsect);
  for s=1:nsect
    sections(:,s)=Y(sections_start_ind(s)+(0:npts-1));
    if 1
      switch datatype
        case 'sim'
          movavgwinsize = 5;
          sel=1:ceil(1.1*Fs*stepL); % exclude post-step section interval
          dy = moving(diff(moving(sections(sel,s),movavgwinsize)),movavgwinsize);
          [~,mloc_off] = max(dy);
          end_sect = sections(max(1,mloc_off-(Fs/4)):mloc_off,s);
          [~,lend] = max(end_sect);
          lend = mloc_off-(Fs/4)+lend;
          StepFin1 = max(1,round(lend-FinPad1));
          StepFin2 = max(1,round(lend-FinPad2));
          End_mV(s) = mean(sections(max(1,StepFin1:StepFin2),s));
        case 'exp'
          tmp=-sections(:,s);
          [pks,locs] = findpeaks(tmp,'MinPeakHeight',prctile(tmp,99.9),'MinPeakDistance',round(Fs*stepL/2));
          lend = locs(end);
          StepFin1 = max(1,round(lend-FinPad1));
          StepFin2 = max(1,round(lend-FinPad2));
          End_mV(s) = mean(sections(max(1,StepFin1:StepFin2),s));   
          %figure; plot(tmp); hold on; plot(locs,pks,'r*'); vline(StepFin1,'k'); vline(StepFin2,'k');
      end    
      VDiff = End_mV(s)-Baseline_mV;
      I_ForDiff = sections_label_num(s);
      Resistance_Mohms(s) = VDiff/I_ForDiff;    
    end
  end

  o1.Baseline_mV = Baseline_mV;
  o1.sections_label_num = sections_label_num;
  o1.step_sections = sections;
  o1.Resistance_Mohms = Resistance_Mohms;
  o1.End_mV = End_mV;

  if plot_flag
    subplot(4,1,2); plot((0:size(sections,1)-1)*dt,sections); axis tight
    hline(o1.Baseline_mV,'r'); vline(StepFin1/Fs,'k'); vline(StepFin2/Fs,'k');
    for k=1:nsect, if ~isnan(End_mV(k)), hline(End_mV(k),'k'); end; end
  end
end

%% Depolarizing steps

%{
OUTPUTS:
o2.Baseline_mV
o2.Spikes_ISIs
o2.Spikes_InstFreq
o2.sections_label_num
o2.y_NoSpike_sect
o2.Spikes_ISI_median

INPUTS:
% {t,intra}
%             offset_voltage: 0
%     tonic_injected_current: 0
%         sections_label_num: [1 2 3 4]
%         sections_start_sec: 4.4900
%        sections_length_sec: 1
%         baseline_start_sec: 0.2500
%        baseline_length_sec: 0.2500
%                step_length: 0.4000

(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec,step_length,plot_flag)

%}

if isfield(P,'CharDepolStepTA') && nsect>0
  t0=P.CharDepolStepTA.sections_start_sec;
  L =P.CharDepolStepTA.sections_length_sec;
  bt0=P.CharDepolStepTA.baseline_start_sec;
  btL=P.CharDepolStepTA.baseline_length_sec;

  baseline=max(1,floor(Fs*bt0))+(0:floor(Fs*btL)-1);
  sections_start_sec = t0+(0:nsect-1)*L;
  sections_start_ind = floor(Fs*sections_start_sec);
  sections_stop_ind = sections_start_ind+floor(Fs*stepL);
  npts = floor(Fs*L);
  if isnan(baseline), baseline=2; end

  Baseline_mV        = mean(Y(baseline));
  sections_label_num = P.CharDepolStepTA.sections_label_num;
  StepBeg1=1000;
  StepBeg2=2000;
  clear Spikes_InstFreq Spikes_ISIs
  sections = nan(npts,nsect);
  nospikes = ones(1,nsect);
  [Spikes_ISIs{1:nsect}] = deal(nan);
  [Spikes_InstFreq{1:nsect}] = deal(nan);
  Spikes_ISI_median = nan(1,nsect);
  End_mV = nan(1,nsect);
  Beg_mV = nan(1,nsect);
  Resistance_Mohms = nan(1,nsect);
  for s=1:nsect
    sections(:,s)=Y(sections_start_ind(s)+(0:npts-1));
    [pks,locs] = findpeaks(sections(:,s),'MinPeakDistance',10,'MinPeakHeight',0);
    % only consider the first step
    sel=locs<(2*stepL*Fs);
    pks=pks(sel);
    locs=locs(sel);
    if ~isempty(pks)
      nospikes(s)=0;
      if numel(pks)>1
        Spikes_ISIs{s}=diff(locs)/Fs;
        Spikes_ISI_median(s)=median(Spikes_ISIs{s});
        Spikes_InstFreq{s}=1./Spikes_ISIs{s};
      end
    else
      if 1 % compute voltage response for I/V curve
      switch datatype
        case 'sim'
          movavgwinsize = 5;
          sel=1:ceil(1.1*Fs*stepL); % exclude post-step section interval
          dy = moving(diff(moving(sections(sel,s),movavgwinsize)),movavgwinsize);
          [~,mloc_off] = max(dy);
          end_sect = sections(max(1,mloc_off-(Fs/4)):mloc_off,s);
          [~,lend] = max(end_sect);
          lend = mloc_off-(Fs/4)+lend;
          StepFin1 = max(1,round(lend-FinPad1));
          StepFin2 = max(1,round(lend-FinPad2));
          End_mV(s) = mean(sections(max(1,StepFin1:StepFin2),s));
          Beg_mV(s) = mean(sections(min(size(sections,1),StepBeg1:StepBeg2),s));
        case 'exp'
          tmp=-sections(:,s);
          [pks,locs] = findpeaks(tmp,'MinPeakHeight',prctile(tmp,99.9),'MinPeakDistance',round(Fs*stepL/2));
          lbeg = locs(end);
          StepBeg1 = max(1,round(lbeg+FinPad2));
          StepBeg2 = max(1,round(lbeg+FinPad1));
          Beg_mV(s) = mean(sections(min(size(sections,1),StepBeg1:StepBeg2),s));
          tmp=sections(:,s);
          [pks,locs] = findpeaks(tmp,'MinPeakHeight',prctile(tmp,99.9),'MinPeakDistance',round(Fs*stepL/2));        
          lend = locs(end);
          StepFin1 = max(1,round(lend-FinPad1));
          StepFin2 = max(1,round(lend-FinPad2));
          End_mV(s) = mean(sections(max(1,StepFin1:StepFin2),s));   
          %figure; plot(tmp); hold on; plot(locs,pks,'r*'); vline(StepBeg1,'k'); vline(StepBeg2,'k');
      end
        VDiff = End_mV(s)-Baseline_mV;
        I_ForDiff = sections_label_num(s);
        Resistance_Mohms(s) = VDiff/I_ForDiff;    
      end    
    end
  end

  o2.Baseline_mV        = Baseline_mV;
  o2.sections_label_num = sections_label_num;
  o2.step_sections      = sections;
  o2.Spikes_ISIs        = Spikes_ISIs;
  o2.Spikes_InstFreq    = Spikes_InstFreq;
  o2.Spikes_ISI_median  = Spikes_ISI_median;
  o2.y_NoSpike_sect = sections(:,nospikes==1);
  o2.y_Spiking_sect = sections(:,nospikes==0);
  o2.Resistance_Mohms = Resistance_Mohms;
  o2.End_mV = End_mV;
  o2.Beg_mV = Beg_mV;

  if plot_flag
    cnt=0;
    tmpdat=sections;
    for i=1:nsect
      if any(sections(:,i)>0)
        tmpdat(:,i)=tmpdat(:,i)+50*cnt;
        cnt=cnt+1;
      end
    end
    subplot(4,1,3); plot((0:size(tmpdat,1)-1)*dt,tmpdat); axis tight
    hline(o2.Baseline_mV,'r'); vline(StepBeg1/Fs,'k'); vline(StepBeg2/Fs,'k');
    for k=1:nsect, if ~isnan(Beg_mV(k)), hline(Beg_mV(k),'k'); end; end
  end
end

%% Tonic depolarization

%{
OUTPUTS:
o3.Spike_Width.sect_1
o3.SpikeRate_persec.sect_1
o3.SpikeAHP_HalfWidth_ms.sect_1

INPUTS:
% {t,field,intra}
%                bpFiltParms: [5 80]
%                      Notch: []
%                    PlotsYN: 'n'
%             offset_voltage: 0
%     tonic_injected_current: 0
%         sections_label_num: 1
%         sections_start_sec: 8.5000
%        sections_length_sec: 91.5000

(x,yF,yIC,bpFiltParms,Notch,PlotsYN,offset_voltage,tonic_injected_current,varargin)

%}

if isfield(P,'CharDepolTonicSpikesTA')
  t0 = P.CharDepolTonicSpikesTA.sections_start_sec;
  L  = P.CharDepolTonicSpikesTA.sections_length_sec;

  prepad=200; postpad=599; maxspikes=20;
  sections_start_ind = max(1,floor(Fs*t0));
  sections_stop_ind = min(length(Y),sections_start_ind+floor(Fs*L));
  section = Y(sections_start_ind:sections_stop_ind);
  [pks,locs] = findpeaks(section,'MinPeakDistance',10,'MinPeakHeight',0);
  if length(locs)>maxspikes
    locs=locs(end-maxspikes+1:end);
  end
  if ~isempty(pks) 
    % spike rates
    if numel(pks)>1
      Spikes_ISI=diff(locs)/Fs;
      Spikes_InstFreq=1/nanmedian(Spikes_ISI);
      Spikes_ISI_median=median(Spikes_ISI);
    else
      Spikes_InstFreq=[nan];    
    end
    % spike morphology
    nspikes=min(maxspikes,length(pks));
    so=nan(prepad+postpad+1,nspikes);
    amp=nan(1,nspikes);
    Spike_amp_mV=nan(1,nspikes);
    Spike_width_ms=nan(1,nspikes);
    ahp_halfWidth=nan(1,nspikes);
    cnt=1;
    for k=1:nspikes
      try
        if (locs(k)-prepad)>0 && (locs(k)+postpad)<length(section)
          so(:,k) = section(locs(k)-prepad:locs(k)+postpad);
          amp(k) = max(so(:,k));
          [~,SpikeSt_i] = max(diff(diff(so(:,k))));
          SpikeSt_pk = so(SpikeSt_i,k);
          Spike_hh = ((amp(k)-SpikeSt_pk)/2)+SpikeSt_pk;
          search_sect = so(SpikeSt_i:end,k);
          cross_i = crossing(search_sect,1:length(search_sect),Spike_hh);
          Spike_width_ms(k) = (cross_i(2)-cross_i(1))/Fs;
          Spike_amp_mV(k)=amp(k)-SpikeSt_pk;

          [~,SpikeEn_i] = crossing(search_sect,1:length(search_sect),search_sect(1));
          search_sect2 = section(round(locs(k)-prepad+SpikeSt_i+SpikeEn_i(2)-2):end);
          if Fs/2<length(search_sect2)
            ahp_base = min(search_sect2(1:Fs/2));
          else
            ahp_base = min(search_sect2(1:end));
          end
          baseline = nanmean(section(locs(k)-prepad:end));
          ahp_hh = ((baseline-ahp_base)/2)+ahp_base;
          cross_i_ahp = crossing(downsample(search_sect2,10),1:length(downsample(search_sect2,10)),ahp_hh);
          if length(cross_i_ahp)>1
            ahp_halfWidth(k) = ((cross_i_ahp(2)-cross_i_ahp(1))*10)/Fs;
          end
        end
      catch
        break;
      end
    end
  else
    Spike_width_ms=nan;
    Spikes_InstFreq=nan;
    ahp_halfWidth=nan;
    nspikes=0;
  end
  o3.SpikeRate_persec.sect_1 = Spikes_InstFreq;
  o3.Spike_Width.sect_1 = Spike_width_ms;
  o3.SpikeAHP_HalfWidth_ms.sect_1 = ahp_halfWidth;
  o3.Num_Spikes = nspikes;

  if plot_flag
    subplot(4,1,4); plot((0:length(section)-1)*dt,section); axis tight; 
    ylim([min(ylim) 0]);%xlim([0 min(5,max(xlim))]);
  end
end

%% ephys metrics. see: getcharacteristics()

% o1,o2,o3

%{
OUTPUTS:
% - resting membrane potential
% - instant spike frequency on depol step
% - spike width at half height
% - spike rate at threshold
% - time to maximal AHP trough
% - maximal AHP trough size
% - difference between 3rd and 2rd spike ISI on depol step
%}
if nsect>0
  % resting membrane potential [mV]
  Vrest = (o1.Baseline_mV+o2.Baseline_mV)/2;

  % instant spike frequency on depol step [Hz]
  %spikefreq = nanmax(cellfun(@nanmax,o2.Spikes_InstFreq));
  %spikefreq = nanmax(cellfun(@mean,o2.Spikes_InstFreq));
  %spikefreq = nanmean(cellfun(@mean,o2.Spikes_InstFreq));
  spikefreq = nanmean(o2.Spikes_InstFreq{end});
else
  Vrest=nan; spikefreq=nan;
end

% spike width at half height [ms]
spikewidth = nanmean(o3.Spike_Width.sect_1)*1000;

% spike rate at threshold [Hz]
% threshfreq = nanmean(o3.SpikeRate_persec.sect_1);
threshfreq = nanmin(o3.SpikeRate_persec.sect_1);

% time to maximal AHP trough [ms]
AHPtime = nanmean(o3.SpikeAHP_HalfWidth_ms.sect_1)*1000;

AHPsize = nan;
isidiff23 = nan;

ephys_metrics = [AHPtime AHPsize spikewidth spikefreq isidiff23 threshfreq Vrest]';
ephys_labels = {'AHPtime','AHPsize','spikewidth','spikefreq','isidiff23','threshfreq','Vrest'};

if nsect>0

  %% I/V. see: extractfeatures()

  sel=~isnan(o1.End_mV);
  I = -stepI*o1.sections_label_num(sel);     % hyperpolarizing step inputs
  v = o1.End_mV(sel);                             % voltage step response
  %sel=~isnan(o2.End_mV);
  %v = [v o2.End_mV(sel)];                         % voltage step response
  sel=~isnan(o2.Beg_mV);
  v = [v o2.Beg_mV(sel)];                         % voltage step response
  I = [I stepI*o2.sections_label_num(sel)];  % depolarizing step inputs
  [I,ind]=sort(I);
  v=v(ind);
  % [I;v]

  % undo scaling since already in sections_label_num
  if strcmp(datatype,'exp'), I=I/stepI; end

  IVP = polyfit(I,v,1);
  IVm = IVP(1); % slope
  IVb = IVP(2); % intercept

  if plot_flag
    figure('position',[1240 20 450 940]); 
    %subplot(2,1,1); plot(I,v,'b-o',I,polyval(IVP,I),'r*--','markersize',10); xlabel('current'); ylabel('voltage'); legend('simdata','fit'); title('I/V curve');
    %subplot(3,1,1); plot(I,v,'b-o',I,polyval(IVP,I),'r*--','markersize',10); xlabel('current'); ylabel('voltage'); legend('simdata','fit'); title('I/V curve');
    subplot(2,1,1); plot(I,v,'b-o','markersize',10); xlabel('current'); ylabel('voltage'); title('I/V curve (final step)');
  end

  %% f/I. see: extractfeatures()
  % f/I from median rate
  Idep = stepI*o2.sections_label_num;
  fdep = 1./o2.Spikes_ISI_median;
  fdep(isnan(fdep))=0;
  % undo scaling since already in sections_label_num
  if strcmp(datatype,'exp'), Idep=Idep/stepI; end
  % fit f/I curve
  fIP = polyfit(Idep,fdep,1);
  fIm = fIP(1); % slope
  fIb = fIP(2); % intercept

  % f/I from initial ISI
  fdep0=zeros(1,length(Idep));
  sel=~cellfun(@(x)all(isnan(x)),o2.Spikes_InstFreq);
  fdep0(sel) = cellfun(@(x)x(1),o2.Spikes_InstFreq(sel));
  fIP0 = polyfit(Idep,fdep0,1);

  % steady-state f/I
  fdepSS=zeros(1,length(Idep));
  fdepSS(sel) = cellfun(@(x)x(end),o2.Spikes_InstFreq(sel));
  fIPSS = polyfit(Idep,fdepSS,1);


  if plot_flag
    %subplot(2,1,2); plot(Idep,fdep,'b-o',Idep,polyval(fIP,Idep),'r*--','markersize',10); xlabel('current'); ylabel('frequency [Hz]'); legend('simdata','fit'); title('f/I curve');
    %subplot(3,1,2); plot(Idep,fdep0,'b-o',Idep,polyval(fIP0,Idep),'r*--','markersize',10); xlabel('current'); ylabel('frequency [Hz]'); legend('simdata','fit'); title('f/I curve (initial)');
    %subplot(3,1,3); plot(Idep,fdepSS,'b-o',Idep,polyval(fIP0,Idep),'r*--','markersize',10); xlabel('current'); ylabel('frequency [Hz]'); legend('simdata','fit'); title('f/I curve (steady state)');
    subplot(2,1,2); plot(Idep,fdep0,'b-o',Idep,fdepSS,'r*--',Idep,fdep,'g-','markersize',10); xlabel('current'); ylabel('frequency [Hz]'); legend('1st ISI','SS','median'); title('f/I curves (first step)');
  end

  [I;v]
  [Idep;fdep]
  
end

% ephys_labels,ephys_metrics
% [I;v]
% [Idep;fdep]

toc

for i=1:length(ephys_labels), fprintf('%s: %g\n',ephys_labels{i},ephys_metrics(i)); end

%{
hyper=1; depol=2; tonic=3;
o1 = result{hyper}; 
o2 = result{depol};
o3 = result{tonic};

% getcharacteristics()
o1.Baseline_mV
o2.Baseline_mV
o2.Spikes_ISIs
o2.Spikes_InstFreq
o3.Spike_Width.sect_1
o3.SpikeRate_persec.sect_1
o3.SpikeAHP_HalfWidth_ms.sect_1

% I/V
o1.step_sections
o1.sections_label_num
o2.sections_label_num
o2.y_NoSpike_sect

% f/I
o2.Spikes_ISI_median
o2.Spikes_InstFreq

maybe:
o1.Resistance_Mohms

%}