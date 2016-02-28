function Iinj = getStepProtocolStim(dt,isi,nsteps,steptime,stepsize,membranearea,nsections,tonictime,bltime,timelimits,ramprate,depol,stepreset)
if nargin<1
  dt = .01; % ms
end
if nargin<2
  isi = 1000; % ms
end
if nargin<3
  nsteps = 5; % number of current pulses
end
if nargin<4
  steptime = 400; % ms
end
if nargin<5
  stepsize = 100; % pA. typically: 100-500pA (.1-.5nA)
end
if nargin<6
  membranearea = 1500; % um^2. typically: 1000-2000 um2
end
if nargin<7
  nsections = 5; % number of blocks of current pulses
end
if nargin<8
  tonictime = 60000; % ms
end
if nargin<9
  bltime = 100; % ms, baseline duration. baseline = [0 bltime].
end
if nargin<10
  timelimits = [0 nsections*nsteps*isi+bltime+tonictime];
elseif numel(timelimits)==1
  timelimits = [0 timelimits];
end
if nargin<11
  ramprate=5;
end
if nargin<12
  depol=1;
end
if nargin<13
  stepreset=0;
end
  
%stepsize = 100; % pA. typically: 100-500pA (.1-.5nA)
%membranearea = 1500; % um^2. typically: 1000-2000 um2
CF = (1e-6)/(1e-8); % pA/um^2 => uA/cm^2. note: 1um2=1e-8cm2, 1pA=1e-6uA
Iapp = CF*stepsize/membranearea; % uA/cm^2

% baseline
bl = zeros(size(0:dt:bltime));

% step time vector
tstep=0:dt:isi;

% step template
I0 = zeros(size(tstep));
if depol==1
  I0(tstep<steptime)=Iapp;
end
I0 = repmat(I0,[1 nsteps]);

% hyperpolarizing and depolarizing steps
n0=length(I0);
I1=zeros(1,n0*nsections);
I2=zeros(1,n0*nsections);
k=0;
for a=1:nsections
  I1(k+(1:n0))=-a*I0;
  I2(k+(1:n0))=a*I0;
  k=k+n0;
end

% ramp to threshold (then tonic depolarization maintained by holdstim())
I3 = linspace(0,nsections*Iapp*ramprate,round(tonictime/dt));

% combine current injections
% hyperpolarizing reset
if stepreset
  Ireset=-a*I0;
  Iinj = [bl I1 I2 Ireset Ireset Ireset Ireset I3];
else
  Iinj = [bl I1 I2 I3];
end
allt = timelimits(1):dt:timelimits(2);
if length(Iinj)<length(allt)
  Iinj0 = Iinj;
  Iinj = zeros(1,length(allt));
  Iinj(1:length(Iinj0)) = Iinj0;
end
%t = (0:length(Iinj)-1)*dt;
%figure; plot(t/1000,Iinj);
