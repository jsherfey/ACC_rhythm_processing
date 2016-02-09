function [T,E_v,E_iNa_m,E_iNa_h,E_iK_n,If_v,If_iNa_m,If_iNa_h,If_iK_n,Is_v,Is_iNa_m,Is_iNa_h,Is_iK_n,If_E_iAMPA_s,Is_E_iAMPA_s,E_If_iGABAa_s,E_Is_iGABAa_s,E_If_iGABAa_ISYN,E_Is_iGABAa_ISYN,If_E_iAMPA_ISYN,Is_E_iAMPA_ISYN,E_s,If_E_iAMPA_netcon,Is_E_iAMPA_netcon,E_If_iGABAa_netcon,E_Is_iGABAa_netcon]=solve_ode
% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
p=load('params.mat','p'); p=p.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);
% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
E_s = get_input('poisson',p.E_Npop,T,p.E_f,p.E_dc,p.E_ac,p.E_tau,p.E_xc,p.E_baseline,p.E_phase);
If_E_iAMPA_netcon = [+1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000];
Is_E_iAMPA_netcon = [+1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000; +1.00000000   +1.00000000];
E_If_iGABAa_netcon = [+1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000; +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000];
E_Is_iGABAa_netcon = [+1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000; +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000   +1.00000000];
% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng(p.random_seed);
t=0;
% STATE_VARIABLES:
E_v = zeros(nsamp,p.E_Npop);
  E_v(1,:) = -65*ones(1,p.E_Npop);
E_iNa_m = zeros(nsamp,p.E_Npop);
  E_iNa_m(1,:) = p.E_iNa_m_IC+p.E_iNa_IC_noise*rand(1,p.E_Npop);
E_iNa_h = zeros(nsamp,p.E_Npop);
  E_iNa_h(1,:) = p.E_iNa_h_IC+p.E_iNa_IC_noise*rand(1,p.E_Npop);
E_iK_n = zeros(nsamp,p.E_Npop);
  E_iK_n(1,:) = p.E_iK_n_IC+p.E_iK_IC_noise*rand(1,p.E_Npop);
If_v = zeros(nsamp,p.If_Npop);
  If_v(1,:) = -65*ones(1,p.If_Npop);
If_iNa_m = zeros(nsamp,p.If_Npop);
  If_iNa_m(1,:) = p.If_iNa_m_IC+p.If_iNa_IC_noise*rand(1,p.If_Npop);
If_iNa_h = zeros(nsamp,p.If_Npop);
  If_iNa_h(1,:) = p.If_iNa_h_IC+p.If_iNa_IC_noise*rand(1,p.If_Npop);
If_iK_n = zeros(nsamp,p.If_Npop);
  If_iK_n(1,:) = p.If_iK_n_IC+p.If_iK_IC_noise*rand(1,p.If_Npop);
Is_v = zeros(nsamp,p.Is_Npop);
  Is_v(1,:) = -65*ones(1,p.Is_Npop);
Is_iNa_m = zeros(nsamp,p.Is_Npop);
  Is_iNa_m(1,:) = p.Is_iNa_m_IC+p.Is_iNa_IC_noise*rand(1,p.Is_Npop);
Is_iNa_h = zeros(nsamp,p.Is_Npop);
  Is_iNa_h(1,:) = p.Is_iNa_h_IC+p.Is_iNa_IC_noise*rand(1,p.Is_Npop);
Is_iK_n = zeros(nsamp,p.Is_Npop);
  Is_iK_n(1,:) = p.Is_iK_n_IC+p.Is_iK_IC_noise*rand(1,p.Is_Npop);
If_E_iAMPA_s = zeros(nsamp,p.E_Npop);
  If_E_iAMPA_s(1,:) =  p.If_E_iAMPA_IC+p.If_E_iAMPA_IC_noise.*rand(1,p.E_Npop);
Is_E_iAMPA_s = zeros(nsamp,p.E_Npop);
  Is_E_iAMPA_s(1,:) =  p.Is_E_iAMPA_IC+p.Is_E_iAMPA_IC_noise.*rand(1,p.E_Npop);
E_If_iGABAa_s = zeros(nsamp,p.If_Npop);
  E_If_iGABAa_s(1,:) =  p.E_If_iGABAa_IC+p.E_If_iGABAa_IC_noise.*rand(1,p.If_Npop);
E_Is_iGABAa_s = zeros(nsamp,p.Is_Npop);
  E_Is_iGABAa_s(1,:) =  p.E_Is_iGABAa_IC+p.E_Is_iGABAa_IC_noise.*rand(1,p.Is_Npop);
% MONITORS:
E_If_iGABAa_ISYN = zeros(nsamp,p.E_Npop);
  E_If_iGABAa_ISYN(1,:)=(p.E_If_iGABAa_gSYN.*(E_If_iGABAa_s(1,:)*E_If_iGABAa_netcon).*(E_v(1,:)-p.E_If_iGABAa_ESYN));
E_Is_iGABAa_ISYN = zeros(nsamp,p.E_Npop);
  E_Is_iGABAa_ISYN(1,:)=(p.E_Is_iGABAa_gSYN.*(E_Is_iGABAa_s(1,:)*E_Is_iGABAa_netcon).*(E_v(1,:)-p.E_Is_iGABAa_ESYN));
If_E_iAMPA_ISYN = zeros(nsamp,p.If_Npop);
  If_E_iAMPA_ISYN(1,:)=(p.If_E_iAMPA_gSYN.*(If_E_iAMPA_s(1,:)*If_E_iAMPA_netcon).*(If_v(1,:)-p.If_E_iAMPA_ESYN));
Is_E_iAMPA_ISYN = zeros(nsamp,p.Is_Npop);
  Is_E_iAMPA_ISYN(1,:)=(p.Is_E_iAMPA_gSYN.*(Is_E_iAMPA_s(1,:)*Is_E_iAMPA_netcon).*(Is_v(1,:)-p.Is_E_iAMPA_ESYN));
% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  E_v_k1=((-(( p.E_iNa_gNa.*E_iNa_m(n-1,:).^3.*E_iNa_h(n-1,:).*(E_v(n-1,:)-p.E_iNa_ENa))))+((-(( p.E_iK_gK.*E_iK_n(n-1,:).^4.*(E_v(n-1,:)-p.E_iK_EK))))+((-(( (p.E_If_iGABAa_gSYN.*(E_If_iGABAa_s(n-1,:)*E_If_iGABAa_netcon).*(E_v(n-1,:)-p.E_If_iGABAa_ESYN)))))+((-(( (p.E_Is_iGABAa_gSYN.*(E_Is_iGABAa_s(n-1,:)*E_Is_iGABAa_netcon).*(E_v(n-1,:)-p.E_Is_iGABAa_ESYN)))))))))+p.E_input-p.E_gINPUT*E_s(k,:).*(E_v(n-1,:)-0)+0*14;
  E_iNa_m_k1= (( (2.5-.1*(E_v(n-1,:)+65))./(exp(2.5-.1*(E_v(n-1,:)+65))-1))).*(1-E_iNa_m(n-1,:))-(( 4*exp(-(E_v(n-1,:)+65)/18))).*E_iNa_m(n-1,:);
  E_iNa_h_k1= (( .07*exp(-(E_v(n-1,:)+65)/20))).*(1-E_iNa_h(n-1,:))-(( 1./(exp(3-.1*(E_v(n-1,:)+65))+1))).*E_iNa_h(n-1,:);
  E_iK_n_k1= (( (.1-.01*(E_v(n-1,:)+65))./(exp(1-.1*(E_v(n-1,:)+65))-1))).*(1-E_iK_n(n-1,:))-(( .125*exp(-(E_v(n-1,:)+65)/80))).*E_iK_n(n-1,:);
  If_v_k1=((-(( p.If_iNa_gNa.*If_iNa_m(n-1,:).^3.*If_iNa_h(n-1,:).*(If_v(n-1,:)-p.If_iNa_ENa))))+((-(( p.If_iK_gK.*If_iK_n(n-1,:).^4.*(If_v(n-1,:)-p.If_iK_EK))))+((-(( (p.If_E_iAMPA_gSYN.*(If_E_iAMPA_s(n-1,:)*If_E_iAMPA_netcon).*(If_v(n-1,:)-p.If_E_iAMPA_ESYN))))))))+p.If_input;
  If_iNa_m_k1= (( (2.5-.1*(If_v(n-1,:)+65))./(exp(2.5-.1*(If_v(n-1,:)+65))-1))).*(1-If_iNa_m(n-1,:))-(( 4*exp(-(If_v(n-1,:)+65)/18))).*If_iNa_m(n-1,:);
  If_iNa_h_k1= (( .07*exp(-(If_v(n-1,:)+65)/20))).*(1-If_iNa_h(n-1,:))-(( 1./(exp(3-.1*(If_v(n-1,:)+65))+1))).*If_iNa_h(n-1,:);
  If_iK_n_k1= (( (.1-.01*(If_v(n-1,:)+65))./(exp(1-.1*(If_v(n-1,:)+65))-1))).*(1-If_iK_n(n-1,:))-(( .125*exp(-(If_v(n-1,:)+65)/80))).*If_iK_n(n-1,:);
  Is_v_k1=((-(( p.Is_iNa_gNa.*Is_iNa_m(n-1,:).^3.*Is_iNa_h(n-1,:).*(Is_v(n-1,:)-p.Is_iNa_ENa))))+((-(( p.Is_iK_gK.*Is_iK_n(n-1,:).^4.*(Is_v(n-1,:)-p.Is_iK_EK))))+((-(( (p.Is_E_iAMPA_gSYN.*(Is_E_iAMPA_s(n-1,:)*Is_E_iAMPA_netcon).*(Is_v(n-1,:)-p.Is_E_iAMPA_ESYN))))))))+p.Is_input;
  Is_iNa_m_k1= (( (2.5-.1*(Is_v(n-1,:)+65))./(exp(2.5-.1*(Is_v(n-1,:)+65))-1))).*(1-Is_iNa_m(n-1,:))-(( 4*exp(-(Is_v(n-1,:)+65)/18))).*Is_iNa_m(n-1,:);
  Is_iNa_h_k1= (( .07*exp(-(Is_v(n-1,:)+65)/20))).*(1-Is_iNa_h(n-1,:))-(( 1./(exp(3-.1*(Is_v(n-1,:)+65))+1))).*Is_iNa_h(n-1,:);
  Is_iK_n_k1= (( (.1-.01*(Is_v(n-1,:)+65))./(exp(1-.1*(Is_v(n-1,:)+65))-1))).*(1-Is_iK_n(n-1,:))-(( .125*exp(-(Is_v(n-1,:)+65)/80))).*Is_iK_n(n-1,:);
  If_E_iAMPA_s_k1= -If_E_iAMPA_s(n-1,:)./p.If_E_iAMPA_tauD + ((1-If_E_iAMPA_s(n-1,:))/p.If_E_iAMPA_tauR).*(1+tanh(E_v(n-1,:)/10));
  Is_E_iAMPA_s_k1= -Is_E_iAMPA_s(n-1,:)./p.Is_E_iAMPA_tauD + ((1-Is_E_iAMPA_s(n-1,:))/p.Is_E_iAMPA_tauR).*(1+tanh(E_v(n-1,:)/10));
  E_If_iGABAa_s_k1= -E_If_iGABAa_s(n-1,:)./p.E_If_iGABAa_tauD + ((1-E_If_iGABAa_s(n-1,:))/p.E_If_iGABAa_tauR).*(1+tanh(If_v(n-1,:)/10));
  E_Is_iGABAa_s_k1= -E_Is_iGABAa_s(n-1,:)./p.E_Is_iGABAa_tauD + ((1-E_Is_iGABAa_s(n-1,:))/p.E_Is_iGABAa_tauR).*(1+tanh(Is_v(n-1,:)/10));
  t=t+.5*dt;
  E_v_k2=((-(( p.E_iNa_gNa.*E_iNa_m(n-1,:).^3.*E_iNa_h(n-1,:).*(E_v(n-1,:)-p.E_iNa_ENa))))+((-(( p.E_iK_gK.*E_iK_n(n-1,:).^4.*(E_v(n-1,:)-p.E_iK_EK))))+((-(( (p.E_If_iGABAa_gSYN.*(E_If_iGABAa_s(n-1,:)*E_If_iGABAa_netcon).*(E_v(n-1,:)-p.E_If_iGABAa_ESYN)))))+((-(( (p.E_Is_iGABAa_gSYN.*(E_Is_iGABAa_s(n-1,:)*E_Is_iGABAa_netcon).*(E_v(n-1,:)-p.E_Is_iGABAa_ESYN)))))))))+p.E_input-p.E_gINPUT*E_s(k,:).*(E_v(n-1,:)-0)+0*14;
  E_iNa_m_k2= (( (2.5-.1*(E_v(n-1,:)+65))./(exp(2.5-.1*(E_v(n-1,:)+65))-1))).*(1-E_iNa_m(n-1,:))-(( 4*exp(-(E_v(n-1,:)+65)/18))).*E_iNa_m(n-1,:);
  E_iNa_h_k2= (( .07*exp(-(E_v(n-1,:)+65)/20))).*(1-E_iNa_h(n-1,:))-(( 1./(exp(3-.1*(E_v(n-1,:)+65))+1))).*E_iNa_h(n-1,:);
  E_iK_n_k2= (( (.1-.01*(E_v(n-1,:)+65))./(exp(1-.1*(E_v(n-1,:)+65))-1))).*(1-E_iK_n(n-1,:))-(( .125*exp(-(E_v(n-1,:)+65)/80))).*E_iK_n(n-1,:);
  If_v_k2=((-(( p.If_iNa_gNa.*If_iNa_m(n-1,:).^3.*If_iNa_h(n-1,:).*(If_v(n-1,:)-p.If_iNa_ENa))))+((-(( p.If_iK_gK.*If_iK_n(n-1,:).^4.*(If_v(n-1,:)-p.If_iK_EK))))+((-(( (p.If_E_iAMPA_gSYN.*(If_E_iAMPA_s(n-1,:)*If_E_iAMPA_netcon).*(If_v(n-1,:)-p.If_E_iAMPA_ESYN))))))))+p.If_input;
  If_iNa_m_k2= (( (2.5-.1*(If_v(n-1,:)+65))./(exp(2.5-.1*(If_v(n-1,:)+65))-1))).*(1-If_iNa_m(n-1,:))-(( 4*exp(-(If_v(n-1,:)+65)/18))).*If_iNa_m(n-1,:);
  If_iNa_h_k2= (( .07*exp(-(If_v(n-1,:)+65)/20))).*(1-If_iNa_h(n-1,:))-(( 1./(exp(3-.1*(If_v(n-1,:)+65))+1))).*If_iNa_h(n-1,:);
  If_iK_n_k2= (( (.1-.01*(If_v(n-1,:)+65))./(exp(1-.1*(If_v(n-1,:)+65))-1))).*(1-If_iK_n(n-1,:))-(( .125*exp(-(If_v(n-1,:)+65)/80))).*If_iK_n(n-1,:);
  Is_v_k2=((-(( p.Is_iNa_gNa.*Is_iNa_m(n-1,:).^3.*Is_iNa_h(n-1,:).*(Is_v(n-1,:)-p.Is_iNa_ENa))))+((-(( p.Is_iK_gK.*Is_iK_n(n-1,:).^4.*(Is_v(n-1,:)-p.Is_iK_EK))))+((-(( (p.Is_E_iAMPA_gSYN.*(Is_E_iAMPA_s(n-1,:)*Is_E_iAMPA_netcon).*(Is_v(n-1,:)-p.Is_E_iAMPA_ESYN))))))))+p.Is_input;
  Is_iNa_m_k2= (( (2.5-.1*(Is_v(n-1,:)+65))./(exp(2.5-.1*(Is_v(n-1,:)+65))-1))).*(1-Is_iNa_m(n-1,:))-(( 4*exp(-(Is_v(n-1,:)+65)/18))).*Is_iNa_m(n-1,:);
  Is_iNa_h_k2= (( .07*exp(-(Is_v(n-1,:)+65)/20))).*(1-Is_iNa_h(n-1,:))-(( 1./(exp(3-.1*(Is_v(n-1,:)+65))+1))).*Is_iNa_h(n-1,:);
  Is_iK_n_k2= (( (.1-.01*(Is_v(n-1,:)+65))./(exp(1-.1*(Is_v(n-1,:)+65))-1))).*(1-Is_iK_n(n-1,:))-(( .125*exp(-(Is_v(n-1,:)+65)/80))).*Is_iK_n(n-1,:);
  If_E_iAMPA_s_k2= -If_E_iAMPA_s(n-1,:)./p.If_E_iAMPA_tauD + ((1-If_E_iAMPA_s(n-1,:))/p.If_E_iAMPA_tauR).*(1+tanh(E_v(n-1,:)/10));
  Is_E_iAMPA_s_k2= -Is_E_iAMPA_s(n-1,:)./p.Is_E_iAMPA_tauD + ((1-Is_E_iAMPA_s(n-1,:))/p.Is_E_iAMPA_tauR).*(1+tanh(E_v(n-1,:)/10));
  E_If_iGABAa_s_k2= -E_If_iGABAa_s(n-1,:)./p.E_If_iGABAa_tauD + ((1-E_If_iGABAa_s(n-1,:))/p.E_If_iGABAa_tauR).*(1+tanh(If_v(n-1,:)/10));
  E_Is_iGABAa_s_k2= -E_Is_iGABAa_s(n-1,:)./p.E_Is_iGABAa_tauD + ((1-E_Is_iGABAa_s(n-1,:))/p.E_Is_iGABAa_tauR).*(1+tanh(Is_v(n-1,:)/10));
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  E_v(n,:)=E_v(n-1,:)+dt*E_v_k2;
  E_iNa_m(n,:)=E_iNa_m(n-1,:)+dt*E_iNa_m_k2;
  E_iNa_h(n,:)=E_iNa_h(n-1,:)+dt*E_iNa_h_k2;
  E_iK_n(n,:)=E_iK_n(n-1,:)+dt*E_iK_n_k2;
  If_v(n,:)=If_v(n-1,:)+dt*If_v_k2;
  If_iNa_m(n,:)=If_iNa_m(n-1,:)+dt*If_iNa_m_k2;
  If_iNa_h(n,:)=If_iNa_h(n-1,:)+dt*If_iNa_h_k2;
  If_iK_n(n,:)=If_iK_n(n-1,:)+dt*If_iK_n_k2;
  Is_v(n,:)=Is_v(n-1,:)+dt*Is_v_k2;
  Is_iNa_m(n,:)=Is_iNa_m(n-1,:)+dt*Is_iNa_m_k2;
  Is_iNa_h(n,:)=Is_iNa_h(n-1,:)+dt*Is_iNa_h_k2;
  Is_iK_n(n,:)=Is_iK_n(n-1,:)+dt*Is_iK_n_k2;
  If_E_iAMPA_s(n,:)=If_E_iAMPA_s(n-1,:)+dt*If_E_iAMPA_s_k2;
  Is_E_iAMPA_s(n,:)=Is_E_iAMPA_s(n-1,:)+dt*Is_E_iAMPA_s_k2;
  E_If_iGABAa_s(n,:)=E_If_iGABAa_s(n-1,:)+dt*E_If_iGABAa_s_k2;
  E_Is_iGABAa_s(n,:)=E_Is_iGABAa_s(n-1,:)+dt*E_Is_iGABAa_s_k2;
  % ------------------------------------------------------------
  % Update monitors:
  % ------------------------------------------------------------
  E_If_iGABAa_ISYN(n,:)=(p.E_If_iGABAa_gSYN.*(E_If_iGABAa_s(n,:)*E_If_iGABAa_netcon).*(E_v(n,:)-p.E_If_iGABAa_ESYN));
  E_Is_iGABAa_ISYN(n,:)=(p.E_Is_iGABAa_gSYN.*(E_Is_iGABAa_s(n,:)*E_Is_iGABAa_netcon).*(E_v(n,:)-p.E_Is_iGABAa_ESYN));
  If_E_iAMPA_ISYN(n,:)=(p.If_E_iAMPA_gSYN.*(If_E_iAMPA_s(n,:)*If_E_iAMPA_netcon).*(If_v(n,:)-p.If_E_iAMPA_ESYN));
  Is_E_iAMPA_ISYN(n,:)=(p.Is_E_iAMPA_gSYN.*(Is_E_iAMPA_s(n,:)*Is_E_iAMPA_netcon).*(Is_v(n,:)-p.Is_E_iAMPA_ESYN));
  n=n+1;
end
