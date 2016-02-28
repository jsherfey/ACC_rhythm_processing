function [P,f,stats]=js_pwelch(X,t,Fmin,Fmax,PxxSmooth,powthreshprc,Fwin)
% usage: [P,f,stats]=js_pwelch(X,t,2,150,5,95,5);
% Fmin=2;
% Fmax=150;
% PxxSmooth=5;
% powthreshprc=95;
% Fwin=5;
% [P,f,stats]=js_pwelch(X,t,Fmin,Fmax,PxxSmooth,powthreshprc,Fwin)
% stats.OscFreq
% stats.AreaPower

if nargin<3, Fmin=2; end
if nargin<4, Fmax=inf; end
if nargin<5, PxxSmooth=1; end
if nargin<6, powthreshprc=95; end
if nargin<7, Fwin=5; end

dt = t(2)-t(1);
Fs = fix(1/dt);
NFFT=2^(nextpow2(length(X)-1)-2);
WINDOW=2^(nextpow2(NFFT-1)-3);
FreqRange=[max(Fmin,2/t(end)) Fmax]; % frequencies to consider
X=smooth(X,ceil(.1/dt));
[P,f] = pwelch(detrend(X),NFFT,[],NFFT,Fs);
sel = find(FreqRange(1)<=f & f<=FreqRange(end));
tmpPxx=smooth(P,PxxSmooth);
ht=prctile(log10(tmpPxx(sel)),powthreshprc);
[PeakPower,PPind]=findpeaks(log10(tmpPxx(sel)),'MinPeakHeight',ht,'NPeaks',3);
if ~isempty(PPind)
  %PPind=PPind(1);
  PPind=PPind(PeakPower==max(PeakPower));
  OscFreq = f(sel(PPind));
  flo=OscFreq-Fwin/2;
  fhi=OscFreq+Fwin/2;
  sel2=find(flo<=f & f<=fhi);
  AreaPower = sum(P(sel2))*(f(2)-f(1));
else
  OscFreq=nan;
  AreaPower=nan;
end
stats.OscFreq=OscFreq;
stats.AreaPower=AreaPower;
