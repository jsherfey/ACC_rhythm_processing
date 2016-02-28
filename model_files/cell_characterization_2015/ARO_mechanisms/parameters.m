
%{

tauI=5

NE=80:
  noiseless_flag=0;
  gEI=8;
  gIE=(nI/nE)*gEI;
  gII=.5*gIE;
  gEEnmda=0;          % between-layer E->E gNMDA
  gEEampa=.1;%.1*gEEnmda; % between-layer E->E gAMPA
  geeampa=.1; % within-layer E->E gAMPA
  geenmda=.5; % within-layer E->E gNMDA

NE=20:
  noiseless_flag=1;
  gEI=.05;
  gIE=.5;
  gII=.5;
  gEEnmda=0;
  gEEampa=.1;
  geeampa=0;
  geenmda=0;

%}

vals=[0 .25 .5 .75 1];
figure
for i=1:length(vals)
  G = InputGenerator2(...
  'TIMELIMITS',[0 1000],'DT',(0.01),'INPUTTYPE',(6),'NPOP',(20),'AMP',(0.3),...
  'ONSET',(50),'OFFSET',(Inf),'SHAREDFRACTION',vals(i),'NINPUTS',(10),...
  'TAUD',(2),'TAUR',(0.5),'LAMBDA',(5));
  subplot(length(vals),1,i); imagesc(G); title(num2str(vals(i))); colorbar
  mean(G(:))
end

fI = @(t,v) (((1)/(10)).*G(:,max(1,round(t/(0.01))))).*(v - (0));
V = squeeze(data(1).epochs.data(3,:,:));
t = 1000*data(1).epochs.time;
I = zeros(size(V));
for i=1:length(t)
  I(i,:) = fI(t(i),V(i,:)');
end
figure; imagesc(t,1:size(V,2),I'); colorbar




