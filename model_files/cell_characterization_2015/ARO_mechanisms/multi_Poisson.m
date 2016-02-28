function psps = multi_Poisson(no_cells,inputs_per_cell,rate,tau_i,tau_1,tau_d,tau_r,T,dt)

% NE = # E-cells (in input layer)
% ne = # excitatory inputs to each E cell

% 1. Generate NE*ne independent Poisson inputs and ne additional inputs
% 2. Replace c of ne inputs by the shared c inputs (for each E cell)

% iMultiPoissonExp.txt:
% rate=2;     % Hz, Poisson spike rate
% T=2000;     % ms, duration
% g_esyn=1;		% input strength
% g_isyn=1;
% E_esyn=0;
% E_isyn=-85;
% tau_i = 10;
% tau_1 = 1;
% N_einputs = 127;
% N_iinputs = 73;
% Ge = multi_Poisson(Npop,N_einputs,rate,tau_i,tau_1,2,.5,T,dt);
t = 0:dt:T;
% EPSP for spikes at time t = 0.
psp = tau_i*(exp(-max(t-tau_1,0)/tau_d) - exp(-max(t-tau_1,0)/tau_r))/(tau_d-tau_r);
%psp = psp(psp > eps);
psp = [zeros(1,length(psp)) psp];
% input connectivity
no_inputs = inputs_per_cell*no_cells;
C = repmat(eye(no_cells), 1, no_inputs/no_cells);
% input spikes
spikes = rand(no_inputs, ceil(T/dt));
spikes = spikes < rate*dt/1000;
spike_arrivals = C*spikes; % Calculating presynaptic spikes for each cell.
% convolve spike trains with EPSP shape
psps = nan(size(spike_arrivals)); % Calculating EPSP experienced by each cell.
for c = 1:no_cells
  psps(c,:) = conv(spike_arrivals(c,:),psp,'same');
end

% getPoissonExp.m:
% if nargin<8, poiss_id=1; end
% if nargin<7, overwrite_flag=0; end
% if nargin<6, dt=.01; end
% if nargin<5, T=1000; end
% if nargin<4, N=1; end
% if nargin<3, Pmax=1; end
% if nargin<2, tauD=10; end
% if nargin<1, lambda=50; end
% nt=ceil(T/dt);
% spikes=poissrnd(lambda*(1e-4),[N nt]);
% G=zeros(N,nt);
% for t=2:nt
%   G(:,t)=G(:,t-1) - dt*G(:,t-1)/tauD;
%   if independent
%     for i=1:N
%       if spikes(i,t)==1
%         G(i,t) = G(i,t) + Pmax.*(1-G(i,t));
%       end
%     end
%   else
%     if spikes(t)==1
%       G(:,t) = G(:,t) + Pmax.*(1-G(:,t));
%     end
%   end
% end
% G=single(G);

