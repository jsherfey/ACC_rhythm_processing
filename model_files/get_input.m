function input=get_input(type,ncells,T,f,dc,ac,tau,xc,baseline,phase)
% inputs:
% - type {'sin'}
% - ncells: size of population
% - T: full time vector [ms] [time x 1]
% - f: oscillation frequency [Hz]
% - phase: phase offset [radians]
% 
% example:
% type='rectified_sin';
% ncells=1; % pop size
% T=0:.01:1000; % ms
% f=5; % Hz
% phase=pi; % radians
% I=get_input(type,ncells,T,f,phase);
% figure; plot(T,I);
% 
% s=get_input('poisson',8,T,35,0,1,2,.25,.1,0);

% default inputs
if nargin<10, phase=0; end
if nargin<9, baseline=.1; end % kHz
if nargin<8, xc=.5; end
if nargin<7, tau=2; end % ms
if nargin<6, ac=1; end % kHz (equivalent to 1000 inputs at 1Hz)
if nargin<5, dc=0; end % kHz
if nargin<4, f=10; end % Hz
if nargin<3, T=(0:.01:1000)'; end % ms
if nargin<2, ncells=1; end
if nargin<1 || ~ischar(type), type='sin'; end

% check time dimensions: should be [T] = time x 1
if size(T,1)<size(T,2)
  error('time vector dimensions must be [time x 1]');
  %T=T'; % <-- this is not compatible with coder
end

ntime=length(T); % # time points

% define input
switch type
  case 'sin'
    T=T/1000; % convert time from ms to sec
    input=sin(2*pi*T*f+phase);
    input=repmat(input,[1 ncells]);
  case 'rectified_sin'
    T=T/1000; % convert time from ms to sec
    input=max(0,sin(2*pi*T*f+phase));
    input=repmat(input,[1 ncells]);
  case 'noisy_sin'
    T=T/1000; dt=T(2)-T(1); % time step, sec
    input=nan(ntime,ncells);
    for i=1:ncells
      input(:,i)=10*sin(2*pi*T*f+phase)+.01*randn(ntime,ncells)*sqrt(dt)/dt;
    end
  case 'noisy_rectified_sin'
    T=T/1000; dt=T(2)-T(1); % time step, sec
    input=nan(ntime,ncells);
    for i=1:ncells
      input(:,i)=max(0,sin(2*pi*T*f+phase)+.0002*randn(ntime,1)*sqrt(dt)/dt);
    end
  case 'poisson'
    on=T(1); 
    off=T(end); 
    dt=T(2)-T(1);
    interval=T(end)-T(1);
    latency=.1; 
    kick=1; 
    fspread=.03; 
    consigma=.001;
    s=getGenExtPoissonTotalGating(on,off,latency,f/1000,fspread,phase,consigma,baseline,dc,ac,tau,kick,ncells,interval,dt,xc);
    input=s';
end

%{
kick=.1;
%f=10; % Hz
%phase=0; % degrees?
if f==0
  % generate homogeneous (DC) spikes
  spikes=poissrnd(dc*dt,[ntime ncells]); % [time x cells x inputs]
else
  % conductance decay time
  %tau=2/1000; % ms, AMPA
  % parameters of time-varying poisson rate      
  %dc=.05; % Hz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
  %ac=1; % Hz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
  % define AC component
  wave=sin(2*pi*f*T+2*pi*phase);    
  % convert [-1,1] to [0,1]
  %modulation=(1+wave)/2;
  modulation=max(wave,0);
  % define time-varying poisson rate
  lambda=dc+ac*modulation;
  % copy rate for each source (input) and target (cell)
  lambda=repmat(lambda,[1 ncells]); % [time x cells]
  % generate nonhomogeneous (DC+AC) spikes
  spikes=poissrnd(lambda); % [time x cells x inputs]
end
% generate spike-triggered conductance changes (i.e., exponential response to poisson inputs)
S=zeros(ntime,ncells);
for k=2:ntime
  % conductance decay
  S(k,:)=S(k-1,:)-dt*S(k-1,:)/tau;
  for i=1:ncells
    nspikes=spikes(k,i);
    if nspikes>=1
      % add a kick for each spike
      S(k,i)=S(k,i)+(kick*nspikes)*(1-S(k,i));
    end
  end
end
input=S;
%}
