function input=get_input(type,ncells,T,f,dc,ac,tau,conn,baseline,phase,kick)
% input=get_input(type,ncells,T,f,dc,ac,tau,conn,baseline,phase)
% arguments:
% - type {'poisson','sin','rectified_sin'}
% - ncells: size of population
% - T: full time vector [ms] [time x 1]
% - f: oscillation frequency [Hz]
% - dc
% - ac
% - tau
% - conn [1 x num_cells], input connectivity
% - baseline
% - phase: phase offset [radians]
% - kick
% 
% Example: simple rectified sin (generates "injected current")
% type='rectified_sin';
% ncells=1; % pop size
% T=(0:.01:1000)'; % ms
% f=5; % Hz
% I=get_input(type,ncells,T,f);
% figure; plot(T,I);
% 
% % Poisson examples (generates "poisson-based time-varying conductance s(t)")
% % Example: Baseline activity and AC everywhere
% T=(0:.01:1000)'; % ms
% N=8; f=35; dc=0; ac=1; tau=2; conn=ones(1,N); baseline=.1; phase=0;
% s=get_input('poisson',N,T,f,dc,ac,tau,conn,baseline,phase);
% PlotData(s);
% % Example: Baseline activity everywhere with restricted AC
% N=4; f=35; dc=0; ac=1; tau=2; conn=[1 1 0 0]; baseline=.1; phase=0;
% s=get_input('poisson',N,T,f,dc,ac,tau,conn,baseline,phase);
% PlotData(s);
% % Example: No baseline activity with restricted DC + AC
% N=4; f=35; dc=.1; ac=1; tau=2; conn=[1 1 0 0]; baseline=0; phase=0;
% s=get_input('poisson',N,T,f,dc,ac,tau,conn,baseline,phase);
% PlotData(s);
% 
% See also: getGenExtPoissonTotalGating

% Notes on Salva's Poisson parameters for PFC inputs to MSN network:
% https://mail.google.com/mail/u/0/#inbox/152ec80bc9304667

% default inputs
if nargin<11, kick=1; end
if nargin<10, phase=0; end
if nargin<9, baseline=.1; end % kHz
% if nargin<8, conn=ones(1,ncells); end
if nargin<7, tau=2; end % ms
if nargin<6, ac=1; end % kHz (equivalent to 1000 inputs at 1Hz)
if nargin<5, dc=0; end % kHz
if nargin<4, f=10; end % Hz
if nargin<3, T=(0:.01:1000)'; end % ms
if nargin<2, ncells=1; end
if nargin<1 || ~ischar(type), type='sin'; end
if nargin<8 || isempty(conn), conn=ones(1,ncells); end

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
    fspread=.03; 
    consigma=.001;
    s=getGenExtPoissonTotalGating(on,off,latency,f/1000,fspread,phase,consigma,baseline,dc,ac,tau,kick,ncells,interval,dt,conn');
    input=s';
end

% Note on getGenExtPoissonTotalGating arguments:
% rows: simultaneous inputs
% cols: different epochs
% 
% Npop=5;
% dt=.01;
% tOn_pfcInp = 0;               % in ms
% tOff_pfcInp = 1000;           % in ms
% latency_pfcInp = [.1;.1];           % in ms (rise and decay time constant for the input)
% freq_pfcInp = [.02;.04];              % kHz, sequence of frequencies, e.g. 0 kHz only DC component, 10 Hz modulation (alpha) or 25 Hz (beta) or 35 Hz (low gamma), or first alpha, then beta/gamma, etc.
% normFreqSigma_pfcInp = [0.03;.03];  % normalized Freq sigma.
% phase_pfcInp = [0;0];             % useful to introduce a phase lag between different components
% widthSigma_pfcInp = 0.001;    % 0.001 represents an abrupt connectivity transition (in contrast to 0.1; it only applies to dc+ac not the baseline)
% rate_pfcInp_baseline = .1;     % kHz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
% rate_pfcInp_dc = [1;1];           % kHz, Poisson spike rate (e.g., 1 kHz may correspond to 1000 external neurons firing at 1 Hz)
% rate_pfcInp_ac = [5;5];
% tau_pfcInp = 2;               % ms, exponential decay time constant (AMPA)
% kick_pfcInp = 1;              % conductance increase after a spike
% interval = 2000;              % interval (ms)
% g_pfcInp = 0;                 % mS/cm^2, external conductance (rate normalized)
% E_pfcInp = 0;                 % reversal potential (mV; AMPA)
% x_c=[.25;.75];
% S = getGenExtPoissonTotalGating(tOn_pfcInp,tOff_pfcInp,latency_pfcInp,freq_pfcInp,normFreqSigma_pfcInp,phase_pfcInp,widthSigma_pfcInp,rate_pfcInp_baseline,rate_pfcInp_dc,rate_pfcInp_ac,tau_pfcInp,kick_pfcInp,Npop,interval,dt,x_c);
% PlotData(S','plot_type','waveform','variable','data');
