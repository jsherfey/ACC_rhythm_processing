for i=1:numLayers
  if i==1
    InputMatrixE{1}=zeros(nE,nt);     % define stimulus (dim1:A/B, dim2:C/D)
    InputMatrixIs{1}=zeros(nIs,nt);   % define rule (r1,r2) - dendritic inhib
    InputMatrixIf{1}=zeros(nI,nt);     % define rule (r1,r2) - somatic inhib
  else
    InputMatrixE{i}=zeros(nE2,nt);
    InputMatrixIf{i}=zeros(nI2,nt);
    InputMatrixIs{i}=zeros(nIs2,nt);
  end
end
if include_thalamus_flag
  InputMatrixTC=zeros(nTC,nt);
  InputMatrixRE=zeros(nRE,nt);
end


bumpscaleE=1/10;
bumpscaleI=1/10;
bumpscaleIs=1/10;
% tip:
  % bumpscale=1: nearly tonic
  % bumpscale<<1: much noisier

if 0 % tonic inputs
  InputMatrixE{1}(SDIM1,nearest(t,S1time(1)):nearest(t,S1time(2)))=StimAmpE*1;
  InputMatrixE{1}(SDIM2,nearest(t,S2time(1)):nearest(t,S2time(2)))=StimAmpE;%.15(soma)%1(dend)
  InputMatrixIf{1}(NOTRULE,nearest(t,Context1Time(1)):nearest(t,Context1Time(2)))=RuleAmpI*1;
  InputMatrixIs{1}(Context1,nearest(t,Context1Time(1)):nearest(t,Context1Time(2)))=RuleAmpI;
  InputMatrixE=[]; InputMatrixI=[]; InputMatrixIs=[];
else % Poisson inputs
  % E input
  
  % 2 bits:
  % E.Stim(1,2) = [CatA CatB]    % (incongruent stimulus)
  % E.Rule(1) = CatA
  % E.Rule(2) = CatB
  
  % 4 bits:
  % E.Stim(1,2) = [CatA CatD]    % (incongruent stimulus)
  % E.Stim(1,2) = [CatB CatD]    % (congruent stimulus)
  % E.Rule(1) = [CatA CatB]
  % E.Rule(2) = [CatC CatD]

  % I.Stim(1,2) = [Context1 Context2]
  % I.Rule(1) = Context1
  % I.Rule(2) = Context2

  % S1time, S2time, Context1Time, Context2Time
  
  % -------------------------
  % Stimulus 1
  targets=EStimTargets;
  S1=InputGenerator2('NPOP',nE,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
  'AMP',StimAmpE,'PAMP',StimAmpE*bumpscaleE,'TARGETS',targets,'ONSET',S1time(1),'OFFSET',S1time(2),...
  'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
  % Stimulus 2
  targets=EStimTargets;
  S2=InputGenerator2('NPOP',nE,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
  'AMP',StimAmpE,'PAMP',StimAmpE*bumpscaleE,'TARGETS',targets,'ONSET',S2time(1),'OFFSET',S2time(2),...
  'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
  % Context(Rule) Cue
%   if RuleAmpI>0 % add input to restore rate (justification: same axons synapse on E and I cells)
    % Context 1
    targets=ERule1Targets;
    R1=InputGenerator2('NPOP',nE,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
    'AMP',RuleAmpE,'PAMP',RuleAmpE*bumpscaleE,'TARGETS',targets,'ONSET',Context1Time(1),'OFFSET',Context1Time(2),...
    'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
    % Context 2
    targets=ERule2Targets;
    R2=InputGenerator2('NPOP',nE,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
    'AMP',RuleAmpE,'PAMP',RuleAmpE*bumpscaleE,'TARGETS',targets,'ONSET',Context2Time(1),'OFFSET',Context2Time(2),...
    'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
    InputMatrixE{1}=S1+S2+R1+R2;
%   else
%     InputMatrixE{1}=S1+S2; %.15(soma)%1(dend)    
%   end
  % -------------------------
  % Ifast input
  % -------------------------
  % Context(Rule) Cue
  % Context 1
  targets=Context1; % NOTRULE
  IfR1=InputGenerator2('NPOP',nI,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
  'AMP',RuleAmpI,'PAMP',RuleAmpI*bumpscaleI,'TARGETS',targets,'ONSET',Context1Time(1),'OFFSET',Context1Time(2),...
  'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
  % Context 2
  targets=Context2; % NOTRULE
  IfR2=InputGenerator2('NPOP',nI,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
  'AMP',RuleAmpI,'PAMP',RuleAmpI*bumpscaleI,'TARGETS',targets,'ONSET',Context2Time(1),'OFFSET',Context2Time(2),...
  'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
  % Stimulus 1
  targets=unique([Context1 Context2]); % NOTRULE
  IfS1=InputGenerator2('NPOP',nI,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
  'AMP',StimAmpI,'PAMP',StimAmpI*bumpscaleI,'TARGETS',targets,'ONSET',S1time(1),'OFFSET',S1time(2),...
  'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);  
  % Stimulus 2
  targets=unique([Context1 Context2]); % NOTRULE
  IfS2=InputGenerator2('NPOP',nI,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
  'AMP',StimAmpI,'PAMP',StimAmpI*bumpscaleI,'TARGETS',targets,'ONSET',S2time(1),'OFFSET',S2time(2),...
  'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);  
  InputMatrixIf{1}=IfS1+IfS2+IfR1+IfR2;
  % -------------------------
  % Islow input
  % -------------------------
  if nIs>0
    % Context(Rule) Cue
    % Context 1
    targets=Context1; % RULE
    IsR1=InputGenerator2('NPOP',nIs,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
    'AMP',RuleAmpIs,'PAMP',RuleAmpIs*bumpscaleIs,'TARGETS',targets,'ONSET',Context1Time(1),'OFFSET',Context1Time(2),...
    'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
    % Context 2
    targets=Context2; % RULE
    IsR2=InputGenerator2('NPOP',nIs,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
    'AMP',RuleAmpIs,'PAMP',RuleAmpIs*bumpscaleIs,'TARGETS',targets,'ONSET',Context2Time(1),'OFFSET',Context2Time(2),...
    'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);
    % Stimulus 1
    targets=unique([Context1 Context2]);
    IsS1=InputGenerator2('NPOP',nIs,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
    'AMP',StimAmpIs,'PAMP',StimAmpIs*bumpscaleIs,'TARGETS',targets,'ONSET',S1time(1),'OFFSET',S1time(2),...
    'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);  
    % Stimulus 2
    targets=unique([Context1 Context2]);
    IsS2=InputGenerator2('NPOP',nIs,'TIMELIMITS',tspan,'DT',dt,'INPUTTYPE',6,...
    'AMP',StimAmpIs,'PAMP',StimAmpIs*bumpscaleIs,'TARGETS',targets,'ONSET',S2time(1),'OFFSET',S2time(2),...
    'LAMBDA',lambda_modulation,'LAMBDA0',lambda_baseline,'WAVETYPE','sin','F0',modulation_frequency,'SHAREDFRACTION',sharedfraction,'NINPUTS',ninputs,'TAUD',inp_tauD,'TAUR',inp_tauR);  
    InputMatrixIs{1}=IsS1+IsS2+IsR1+IsR2;
  end
  if kickstart
    kickindex=5000:5200; kickamp=10;%.75;
    InputMatrixE{1}(:,kickindex)=kickamp;
    InputMatrixE{2}(:,kickindex)=kickamp;
    InputMatrixIf{1}(:,kickindex)=kickamp;
    InputMatrixIf{2}(:,kickindex)=kickamp;
    InputMatrixIs{1}(:,kickindex)=kickamp;
    InputMatrixIs{2}(:,kickindex)=kickamp;
    kickindex=6000:6500; 
    kickamp=-5;
    InputMatrixE{1}(:,kickindex)=kickamp;
    InputMatrixE{2}(:,kickindex)=kickamp;
    InputMatrixIf{1}(:,kickindex)=kickamp;
    InputMatrixIf{2}(:,kickindex)=kickamp;
    InputMatrixIs{1}(:,kickindex)=kickamp;
    InputMatrixIs{2}(:,kickindex)=kickamp;
  end
  % Plot inputs
  if plotinputs
    if nI2>0 || nE2>0
      nrows=3;
    else
      nrows=2;
    end
    figure('position',[250 380 1550 600],'name','Network inputs'); tt=tspan(1):dt:tspan(2); 
    subplot(nrows,2,1); I=InputMatrixE{1}; imagesc(I); title('Iext->E1'); xlabel('time'); ylabel('E1'); colorbar;
    subplot(nrows,2,2); plot(tt,I(1,:),'b-',tt,I(end,:),'r-'); legend('E(1)','E(end)'); xlabel('time'); ylabel('amp'); %xlim([0 4000]);  
    subplot(nrows,2,3); I=InputMatrixIf{1}; imagesc(I); title('Iext->Ifast'); colorbar; xlabel('time'); ylabel('Ifast'); 
    subplot(nrows,2,4); I=InputMatrixIs{1}; imagesc(I); title('Iext->Islow'); colorbar; xlabel('time'); ylabel('Islow');
    if nE2>0
      subplot(nrows,2,5); I=InputMatrixE{2}; imagesc(I); title('Iext->E2'); xlabel('time'); ylabel('E2'); colorbar;      
    end
    if nI2>0
      subplot(nrows,2,6); I=InputMatrixIf{2}; imagesc(I); title('Iext->Ifast2'); xlabel('time'); ylabel('Ifast2'); colorbar;      
    end
  end
end
