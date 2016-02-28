function[Psd, freq, PsdManSmooth, freqbin] = powerspectrum(t,v)

% Computation of the power spectrum (Psd) for a signal v(t)
% [t] = msec
% [f] = Hz

vm = mean(v);
y = v - vm;


T = t(end);
L = length(y);
Nfft = 2^nextpow2(L);
Y = fft(y,Nfft)/(L/2);
Ys = abs(Y(1:Nfft/2));
%Ys = abs(Y(1:Nfft/2).^2);
% if you want to smooth the PSD uncomment the next two lines
% Ysmooth = smooth(Ys,'lowess');
% Ys = Ysmooth;
Fs = L/(T/1000);
freq = Fs/2*linspace(0,1,Nfft/2);
Psd = Ys;

% Manual Smooth (PsdManSmooth) of the "row" power spectrum Psd

fmax = 40;
fgraphcut = find(freq>fmax,1);
binsize = 0.25; %in Hz
freq(2)-freq(1);
binsize = floor(binsize./(freq(2)-freq(1)));
if binsize<=0
    binsize=1;
end


freqbin = zeros(1,ceil(fgraphcut/binsize));
PsdManSmooth = zeros(1,ceil(fgraphcut/binsize));

for k=1:ceil(fgraphcut/binsize)
    freqbin(k) = freq(k*binsize);
    PsdManSmooth(k) = max(Ys(k*binsize:(k+1)*binsize-1));
end


figure
hold on
%plot(freq(1:fgraphcut),Psd(1:fgraphcut),'o','linewidth',1);
plot(freqbin,PsdManSmooth,'o');
axis([0 fmax,-0.01,max(Ys)*1.2]);
set(gca,'fontsize',20);
xlabel('Freq.  [Hz]');
ylabel('|Y(f)|');





