% nBits=2;          % 2 or 4
% nECellsPerBit=5;
% nResponses=2;     % 1 or 2 (eg, Left, Right)
% nECellsPerResponse=2;
% nRules=2;         % 1 or 2 (eg, Color, Orientation)
% nICellsPerRule=2;
% nIsCellsPerRule=2;
% numLayers=2;                % number of layers
% % input layer
% nE=nBits*nECellsPerBit;     % number of E-cells (superficial layer)
% nI=nRules*nICellsPerRule;   % number of fast I-cells (superficial layer)
% nIs=nRules*nIsCellsPerRule; % number of slow I-cells (superficial layer)
% % downstream layers
% nE2=nResponses*nECellsPerResponse; % number of E-cells (downstream layers)
% nI2=1;                         % number of fast I-cells (downstream layers)
% nIs2=0;                        % number of slow I-cells (downstream layers)
% FFfraction=.8; % what random fraction of a Bit's cells connect to each Response cell

% E->E (recurrent) (item maintenance)

% E1->E1: connect E1 cells within Bit
EEmask{1}=zeros(nE,nE);                         % [target x source], E->E connectivity mask. (eg, ones(nE,nE))
bitmat=ones(nECellsPerBit)-eye(nECellsPerBit);  % connectivity within a single bit
for i=1:nBits
  pos=(i-1)*nECellsPerBit; % start point in connectivity matrix for this bit
  ind=pos+(1:nECellsPerBit); % rows/cols for this bit
  EEmask{1}(ind,ind)=bitmat;
end

% E2->E2: connect E2 cells within Response
EEmask{2}=zeros(nE2,nE2); % E->E (L2) (item maintenance)
resmat=ones(nECellsPerResponse)-eye(nECellsPerResponse);
for i=1:nResponses
  pos=(i-1)*nECellsPerResponse; % start point in connectivity matrix for this bit
  ind=pos+(1:nECellsPerResponse); % rows/cols for this bit
  EEmask{2}(ind,ind)=resmat;
end

% E->E (feedforward) (context-dependent S->R)

% ################################################
FFnuminputs=ceil(FFfraction*nECellsPerBit);
FFnumexclude=nECellsPerBit-FFnuminputs;
% note: need to scale gEEampa=(gEEampa/FFfraction) for nECellsPerBit=5
% ################################################

% Probabilistic E1->E2: every other E1 Bit connects to the same E2 Response (works for nResponses <=2)
EEmaskFF{1} = zeros(nE2,nE); % [target x source], E(L1) -> E(L2) (S->R) 
iomap=zeros(nBits,1);
iomap(1:2:end,1)=1; % maps to Response 1
iomap(2:2:end,1)=2; % maps to Response 2
iomat=ones(nECellsPerResponse,nECellsPerBit);
for i=1:nBits
  posi=(i-1)*nECellsPerBit;
  indi=posi+(1:nECellsPerBit);
  poso=(iomap(i)-1)*nECellsPerResponse;
  indo=poso+(1:nECellsPerResponse);
  % limit # inputs to each target based on FFfraction
  thisconn=iomat;
  if i==1 && FFnumexclude>0
    for j=1:nECellsPerResponse
      srcinds=randperm(nECellsPerBit);
      disconnect=srcinds(1:FFnumexclude);
      if j>1 && isequal(disconnect,lastdisconnect)
        while isequal(disconnect,lastdisconnect)
          srcinds=randperm(nECellsPerBit);
          disconnect=srcinds(1:FFnumexclude);
        end
      end
      thisconn(j,disconnect)=0;
      lastdisconnect=disconnect;
    end
    iomat=thisconn;
  end
  EEmaskFF{1}(indo,indi)=thisconn;
end

% Interneuron connections

% E->Ifast
if nBits==1
  EIfmask{1}=ones(nI,nE);
else
  EIfmask{1}=zeros(nI,nE);  
  EIfmask{1}(1:ceil(.5*nI),1:ceil(.5*nE))=1;        % E(0-50%)->PV(0-50%)
  EIfmask{1}(ceil(.5*nI)+1:nI,ceil(.5*nE)+1:nE)=1;  % E(50-100%)->PV(50-100%)
end

% E->Islow
EIsmask{1}=zeros(nIs,nE); % E -> CB (L1)
if nBits==1 && nIs>0
  EIsmask{1}=ones(nIs,nE);
elseif nIs>0
  EIsmask{1}(1:ceil(.5*nIs),1:ceil(.5*nE))=1;         % E(0-50%)->CB(0-50%)
  EIsmask{1}(ceil(.5*nIs)+1:nIs,ceil(.5*nE)+1:nE)=1;  % E(50-100%)->CB(50-100%)
end

% Ifast->Ifast (connect within rule)
% PV -> PV (L1): connect If1 cells within Rule
IfIfmask{1}=zeros(nI,nI);
rulemat=ones(nICellsPerRule);  % connectivity within a single bit
for i=1:nRules
  pos=(i-1)*nICellsPerRule; % start point in connectivity matrix for this bit
  ind=pos+(1:nICellsPerRule); % rows/cols for this bit
  IfIfmask{1}(ind,ind)=rulemat;
end

% Islow->Islow (connect within rule)
% CB -> CB (L1): connect Is1 cells within Rule
IsIsmask{1}=zeros(nIs,nIs);
rulemat=ones(nIsCellsPerRule);  % connectivity within a single bit
for i=1:nRules
  pos=(i-1)*nIsCellsPerRule; % start point in connectivity matrix for this bit
  ind=pos+(1:nIsCellsPerRule); % rows/cols for this bit
  IsIsmask{1}(ind,ind)=rulemat;
end

% Probabilistic E -> Ifast (layer 2, probabilistic WTA)
% ################################################
EInuminputs=ceil(EI2fraction*nECellsPerResponse);
EInumexclude=nECellsPerResponse-EInuminputs;
% note: need to scale gEEampa=(gEEampa/FFfraction) for nECellsPerBit=5
% ################################################
EIfmask{2}=zeros(nI2,nE2); % E -> PV (L2) (semi-all-to-all for WTA)
iomat=ones(nI2,nECellsPerResponse);
for i=1:nResponses
  posi=(i-1)*nECellsPerResponse;
  indi=posi+(1:nECellsPerResponse);
  poso=0; %(iomap(i)-1)*nI2;
  indo=poso+(1:nI2);
  % limit # inputs to each target based on FFfraction
  thisconn=iomat;
  if i==1 && EInumexclude>0
    for j=1:nI2
      srcinds=randperm(nECellsPerResponse);
      disconnect=srcinds(1:EInumexclude);
      if j>1 && isequal(disconnect,lastdisconnect)
        while isequal(disconnect,lastdisconnect)
          srcinds=randperm(nECellsPerResponse);
          disconnect=srcinds(1:EInumexclude);
        end
      end
      thisconn(j,disconnect)=0;
      lastdisconnect=disconnect;
    end
    iomat=thisconn;
  end
  EIfmask{2}(indo,indi)=thisconn;
end


% IfIfmask{1}=eye(nI,nI);       % PV -> PV (L1)
% IsIsmask{1}=eye(nIs,nIs);     % CB -> CB (L1)

IfEmask{1}=EIfmask{1}';       % PV -> E (L1) (reciprocal feedback inhibition)
IfEmask{2}=ones(nE2,nI2);       % PV -> E (L2) (reciprocal feedback inhibition)
IsEmask{1}=EIsmask{1}';       % CB -> E (L1) (reciprocal feedback inhibition)
IfIsmask{1}=zeros(nIs,nI);    % PV -> CB (L1) (no connection)
IsIfmask{1}=zeros(nI,nIs);    % CB -> PV (L1)
IsIfmask{1}(1:ceil(.5*nI),1:ceil(.5*nIs))=1;         % CB(0-50%)->PV(0-50%)
IsIfmask{1}(ceil(.5*nI)+1:nI,ceil(.5*nIs)+1:nIs)=1;  % CB(50-100%)->PV(50-100%)
IfIfmask{2}=ones(nI2,nI2);    % PV -> PV (L2)
IsIsmask{2}=ones(nIs2,nIs2);  % CB -> CB (L2)
% DNE (iff nIs2=0):
EIsmask{2}=EIfmask{2};   % E -> CB (L2)
IsEmask{2}=ones(nE2,nIs2);   % CB -> E (L2)
IfIsmask{2}=zeros(nIs2,nI2);  % PV -> CB (L2)
IsIfmask{2}=ones(nI2,nIs2);  % CB -> PV (L2)

if numLayers>2
  for i=3:numLayers
    EEmaskFF{i-1}=EEmaskFF{1};
    EEmask{i}=EEmask{2};
    EIfmask{i}=EIfmask{2};
    IfEmask{i}=IfEmask{2};
    IsEmask{i}=IsEmask{2};
    IfIsmask{i}=IfIsmask{2};
    IsIfmask{i}=IsIfmask{2};
    IfIfmask{i}=IfIfmask{2};
    IsIsmask{i}=IsIsmask{2};
    EIsmask{i}=EIsmask{2};
    IsEmask{i}=IsEmask{2};
    IfIsmask{i}=IfIsmask{2};
    IsIfmask{i}=IsIfmask{2};    
  end
end

if include_thalamus_flag
  ETCmask = zeros(nTC,nE2); % E(L2) -> TC   [target x source]
  ETCmask(1:ceil(.5*nTC),1:ceil(.5*nE2))=1;            % E2(0-50%)->TC(0-50%)
  if nTC>1
    ETCmask(ceil(.5*nTC)+1:nTC,ceil(.5*nE2)+1:nE2)=1;  % E2(50-100%)->TC(50-100%)
  else
    ETCmask(1,ceil(.5*nE2)+1:nE2)=1;
  end
  TCEmask = zeros(nE,nTC);  % TC -> E(L1)   [target x source]
  TCEmask(1:ceil(.5*nE),1:ceil(.5*nTC))=1;            % TC(0-50%)->E(0-50%)
  if nTC>1
    TCEmask(ceil(.5*nE)+1:nE,ceil(.5*nTC)+1:nTC)=1;  % TC(50-100%)->E(50-100%)
  else
    TCEmask(ceil(.5*nE)+1:nE,1)=1;  % TC(50-100%)->E(50-100%)
  end
  TCREmask = ones(nRE,nTC); % TC -> RE
  RETCmask = ones(nTC,nRE); % RE -> TC
  REREmask = ones(nRE,nRE); % RE -> RE
  if plotconnect
    figure; 
    subplot(1,2,1); imagesc(ETCmask); title('E->TC'); xlabel('E'); ylabel('TC');
    subplot(1,2,2); imagesc(TCEmask); title('TC->E'); xlabel('TC'); ylabel('E');
  end
end

if use_predef_topology
  if nBits~=2 || nECellsPerBit~=5
    error('predefined RNN(1) matrices (EEmaskFF) exist only for nBits=2, nECellsPerBit=5, EI2fraction=1. need to create others first.');
  end
  if EI2fraction<1
    error('predefined matrices (EIfmask,EIsmask) do not exist for EI2fraction<1. need to create them first.');
  end
  switch FFfraction
    case .4
      load('ruletask_ResponseLR_CB_EEmaskFF_40percent.mat','EEmaskFF');
    case .6
      load('ruletask_ResponseLR_CB_EEmaskFF_60percent.mat','EEmaskFF');
    case .8
      load('ruletask_ResponseLR_CB_EEmaskFF_80percent.mat','EEmaskFF');
    case 1
      % do nothing (all-to-all)
    otherwise
      error('no FF matrix defined for this FFfraction.');    
  end
end

% IsIfsame=0;

if plotconnect
  clims=[0 1];
  figure('position',[650 70 600 900],'name','network topology');
  subplot(5,3,1); % layer 1, E<->E
  imagesc(EEmask{1}); title('RNN(1): E<->E'); xlabel('E'); ylabel('E'); caxis(clims);
  subplot(5,3,2); % layer 1->2, E->E
  imagesc(EEmaskFF{1}); title('FF: E->E'); xlabel('E'); ylabel('E'); caxis(clims);
  subplot(5,3,3); % layer 2, E<->E
  imagesc(EEmask{2}); title('RNN(2): E<->E'); xlabel('E'); ylabel('E'); caxis(clims);
  subplot(5,3,4); % layer 1, E->If
  imagesc(EIfmask{1}); title('E->If'); xlabel('E'); ylabel('If'); caxis(clims);
  subplot(5,3,5); % layer 1, If->E
  imagesc(IfEmask{1}); title('If->E'); xlabel('If'); ylabel('E'); caxis(clims);
  subplot(5,3,6); % layer 1, If->If
  imagesc(IfIfmask{1}); title('If->If'); xlabel('If'); ylabel('If'); caxis(clims);
  subplot(5,3,7); % layer 1, E->Is
  imagesc(EIsmask{1}); title('E->Is'); xlabel('E'); ylabel('Is'); caxis(clims);
  subplot(5,3,8); % layer 1, Is->E
  imagesc(IsEmask{1}); title('Is->E'); xlabel('Is'); ylabel('E'); caxis(clims);
  subplot(5,3,9); % layer 1, Is->Is
  imagesc(IsIsmask{1}); title('Is->Is'); xlabel('Is'); ylabel('Is'); caxis(clims);
  subplot(5,3,10); % layer 1, If->Is
  imagesc(IfIsmask{1}); title('If->Is'); xlabel('If'); ylabel('Is'); caxis(clims);
  subplot(5,3,11); % layer 1, Is->If
  imagesc(IsIfmask{1}); title('Is->If'); xlabel('Is'); ylabel('If'); caxis(clims);
  if numLayers>1
    subplot(5,3,12); % layer 2, E->If
    imagesc(EIfmask{2}); title('E->If (layer 2)'); xlabel('E'); ylabel('If'); caxis(clims);
    subplot(5,3,15); % layer 2, E->If
    imagesc(IfEmask{2}); title('If->E (layer 2)'); xlabel('If'); ylabel('E'); caxis(clims);    
    subplot(5,3,13); % layer 2,If->Is
    imagesc(IfIsmask{2}); title('If->Is (layer 2)'); xlabel('If'); ylabel('Is'); caxis(clims);    
    subplot(5,3,14); % layer 2, Is->If
    imagesc(IsIfmask{2}); title('Is->If (layer 2)'); xlabel('Is'); ylabel('If'); caxis(clims);    
  end  
end

% % E->E (recurrent) (item maintenance)
% EEmask{1}=zeros(nE,nE);          % [target x source], E->E connectivity mask. (eg, ones(nE,nE))
% for i=1:nE % two-cell clusters
%   if rem(i,2)==1 % odd
%     EEmask{1}(i+1,i)=1;
%   else
%     EEmask{1}(i-1,i)=1;
%   end
% end
% EEmask{2}=zeros(nE2,nE2); % E->E (L2) (item maintenance)
% % for i=1:nE2
% %   if rem(i,2)==1 % odd
% %     EEmask{2}(i+1,i)=1;
% %   else
% %     EEmask{2}(i-1,i)=1;
% %   end
% % end
% % E->E (feedforward) (context-dependent S->R)
% EEmaskFF{1} = zeros(nE2,nE); % E(L1) -> E(L2) (S->R)
% group=1; sel1=1:ceil(.5*nE2); sel2=ceil(.5*nE2)+1:nE2;
% for i=1:2:nE % alternating E1 clusters to E2(0-50%) vs. E2(50-100%)
%   if group==1 % target 0-50%
%     EEmaskFF{1}(sel1,i:i+1)=1;
%   else        % target 50-100%
%     EEmaskFF{1}(sel2,i:i+1)=1;    
%   end
%   group=-group;
% end

% Interneuron connections
% IsIfsame=0; % INSTRUCTION: set IsIfsame=1 iff IslowType=IfastType (below)
% EIfmask{1}=zeros(nI,nE); % E -> PV (L1)
% if nE<=2
%   EIfmask{1}=ones(nI,nE);
% else
%   if IsIfsame % Is=If (IfastType=IslowType), each to separate half of E-cells
%     EIfmask{1}(1:ceil(.5*nI),1:ceil(.5*nE))=1;        % E(0-50%)->PV(0-50%)
%   else % Is~=If, each to both halves of E-cells
%     EIfmask{1}(1:ceil(.5*nI),1:ceil(.5*nE))=1;        % E(0-50%)->PV(0-50%)
%     EIfmask{1}(ceil(.5*nI)+1:nI,ceil(.5*nE)+1:nE)=1;  % E(50-100%)->PV(50-100%)
%   end
% end
% EIfmask{2}=ones(nI2,nE2); % E -> PV (L2) (all-to-all for WTA)
% EIsmask{1}=zeros(nIs,nE); % E -> CB (L1)
% if nE<=2 && nIs>0
%   EIsmask{1}=ones(nIs,nE);
% elseif nIs>0
%   if IsIfsame % Is=If (IfastType=IslowType), each to separate half of E-cells
%     EIsmask{1}(ceil(.5*nIs)+1:nIs,ceil(.5*nE)+1:nE)=1;  % E(50-100%)->CB(50-100%)
%   else % Is~=If, each to both halves of E-cells
%     EIsmask{1}(1:ceil(.5*nIs),1:ceil(.5*nE))=1;         % E(0-50%)->CB(0-50%)
%     EIsmask{1}(ceil(.5*nIs)+1:nIs,ceil(.5*nE)+1:nE)=1;  % E(50-100%)->CB(50-100%)
%   end
% end
% if 0 % dominant vs. nondominant rule
%   domstim=1; nondomstim=.5;
%   EIsmask{1}(1:ceil(.5*nIs),1:ceil(.5*nE))=domstim;            % DOMINANT rule: E(0-50%)->CB(0-50%)
%   EIsmask{1}(ceil(.5*nIs)+1:nIs,ceil(.5*nE)+1:nE)=nondomstim;  % NONDOMINANT rule: E(50-100%)->CB(50-100%)
% end
% IfEmask{1}=EIfmask{1}';       % PV -> E (L1) (reciprocal feedback inhibition)
% IfEmask{2}=EIfmask{2}';       % PV -> E (L2) (reciprocal feedback inhibition)
% IsEmask{1}=EIsmask{1}';       % CB -> E (L1) (reciprocal feedback inhibition)
% IfIsmask{1}=zeros(nIs,nI);    % PV -> CB (L1) (no connection)
% IsIfmask{1}=zeros(nI,nIs);    % CB -> PV (L1)
% IsIfmask{1}(1:ceil(.5*nI),1:ceil(.5*nIs))=1;         % CB(0-50%)->PV(0-50%)
% IsIfmask{1}(ceil(.5*nI)+1:nI,ceil(.5*nIs)+1:nIs)=1;  % CB(50-100%)->PV(50-100%)
% IfIfmask{1}=eye(nI,nI);       % PV -> PV (L1)
% IfIfmask{2}=ones(nI2,nI2);    % PV -> PV (L2)
% IsIsmask{1}=eye(nIs,nIs);     % CB -> CB (L1)
% IsIsmask{2}=ones(nIs2,nIs2);  % CB -> CB (L2)
% % DNE (iff nIs2=0):
% EIsmask{2}=zeros(nIs2,nE2);   % E -> CB (L2)
% IsEmask{2}=zeros(nE2,nIs2);   % CB -> E (L2)
% IfIsmask{2}=zeros(nIs2,nI2);  % PV -> CB (L2)
% IsIfmask{2}=zeros(nI2,nIs2);  % CB -> PV (L2)
% 
% if plotconnect
%   figure('position',[650 300 600 700],'name','network topology');
%   subplot(4,3,1); % layer 1, E<->E
%   imagesc(EEmask{1}); title('RNN(1): E<->E'); xlabel('E'); ylabel('E');
%   subplot(4,3,2); % layer 1->2, E->E
%   imagesc(EEmaskFF{1}); title('FF: E->E'); xlabel('E'); ylabel('E');
%   subplot(4,3,3); % layer 2, E<->E
%   imagesc(EEmask{2}); title('RNN(2): E<->E'); xlabel('E'); ylabel('E');
%   subplot(4,3,4); % layer 1, E->If
%   imagesc(EIfmask{1}); title('E->If'); xlabel('E'); ylabel('If');
%   subplot(4,3,5); % layer 1, If->E
%   imagesc(IfEmask{1}); title('If->E'); xlabel('If'); ylabel('E');
%   subplot(4,3,6); % layer 1, If->If
%   imagesc(IfIfmask{1}); title('If->If'); xlabel('If'); ylabel('If');  
%   subplot(4,3,7); % layer 1, E->Is
%   imagesc(EIsmask{1}); title('E->Is'); xlabel('E'); ylabel('Is');
%   subplot(4,3,8); % layer 1, Is->E
%   imagesc(IsEmask{1}); title('Is->E'); xlabel('Is'); ylabel('E');
%   subplot(4,3,9); % layer 1, Is->Is
%   imagesc(IsIsmask{1}); title('Is->Is'); xlabel('Is'); ylabel('Is');    
%   subplot(4,3,10); % layer 1, If->Is
%   imagesc(IfIsmask{1}); title('If->Is'); xlabel('If'); ylabel('Is');  
%   subplot(4,3,11); % layer 1, Is->If
%   imagesc(IsIfmask{1}); title('Is->If'); xlabel('Is'); ylabel('If');  
% end


