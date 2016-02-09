function S = getGenExtPoissonTotalGating(tOn,tOff,latency_,freq_,normFreqSigma_,phase_,widthSigma,rate_baseline,rate_dc_,rate_ac_,tau,kick,N,interval,dt,x_c_)

% 04-Feb-2016: JSS added x_c as input argument

  if nargin<1,  tOn = 0;                end
  if nargin<2,  tOff = 0;               end
  if nargin<3,  latency_ = 0;           end
  if nargin<4,  freq_ = 0;              end
  if nargin<5,  normFreqSigma_ = 0.03;  end
  if nargin<6,  phase_ = 0;             end
  if nargin<7,  widthSigma = 0.001;     end % 0.001 represents an abrupt connectivity transition (in contrast to 0.1; it only applies to dc+ac not the baseline)
  if nargin<8,  rate_baseline = 0;      end
  if nargin<9,  rate_dc_ = 0;           end
  if nargin<10, rate_ac_ = 0;           end
  if nargin<11, tau = 2;                end
  if nargin<12, kick = 1;               end
  if nargin<13, N = 100;                end
  if nargin<14, interval = 1000;        end
  if nargin<15, dt = 0.05;              end
  if nargin<16, x_c_=.25;                end

  % broadcasting sizes with respect to tOn if they don't match
  if length(tOn) > size(freq_,2), freq = freq_*ones(1,length(tOn)); else, freq = freq_; end

  if length(tOn) > size(normFreqSigma_,2), normFreqSigma = normFreqSigma_*ones(1,length(tOn)); else, normFreqSigma = normFreqSigma_; end

  if length(tOn) > size(phase_,2), phase = phase_*ones(1,length(tOn)); else, phase = phase_; end

  if length(tOn) > size(rate_dc_,2), rate_dc = rate_dc_*ones(1,length(tOn)); else, rate_dc = rate_dc_; end

  if length(tOn) > size(rate_ac_,2), rate_ac = rate_ac_*ones(1,length(tOn)); else, rate_ac = rate_ac_; end

  if length(tOn) > size(latency_,2), latency = latency_*ones(1,length(tOn)); else, latency = latency_; end

  if length(tOn) > size(x_c_,2), x_c = x_c_*ones(1,length(tOn)); else, x_c = x_c_; end
  
  time = 0:dt:interval;
  S_ini = zeros(N,1);
  S = zeros(N,length(time));
  dynrate = rate_baseline*ones(N,length(time));
  ratecomp = zeros(length(tOn),length(time));
  for i = 1:length(tOn)
    timeWindow = zeros(1,length(time));
    if latency(i) ~= 0
      timeWindow(time >= tOn(i)) = 1;
    else
      timeWindow(time >= tOn(i) & time <= tOff(i)) = 1;
    end
    ratecomp(i,:) = rate_dc(i);
    if ~all(freq(:,i)==0)
      % freq modulation
      numFreqs = 1000;
      step = 2*pi/numFreqs;
      %Ph = load('sinPhases'); % from test.m in this directory
      phase_ic = 0:step:2*pi*(1-1/numFreqs);
      shuffle = randperm(length(phase_ic));
      phase_ic = phase_ic(shuffle);
      Ph.phase_ic=phase_ic;
      % save('sinPhases','phase_ic')
      for j = 1:length(freq(:,i))
        if freq(j,i) ~= 0
          freqSigma = normFreqSigma(i)*freq(j,i); % 0.03*freq(j,i);
          freqSet = -freq(j,i)/5:2*freq(j,i)/5/(numFreqs-1):freq(j,i)/5;
          freqVar = exp(-freqSet.^2/(2*freqSigma^2));
          m = sum(freqVar)/numFreqs;
          sumCos = zeros(size(time));
          for k = 1:numFreqs
            sumCos = sumCos + freqVar(k)*cos(2*pi*freqSet(k)*time + Ph.phase_ic(k));
          end
          ratecomp(i,:) = ratecomp(i,:) + rate_ac(j,i)*sin(2*pi*freq(j,i)*time+m*sumCos+phase(j,i));
        end
      end
    end
    ratecomp(i,:) = ratecomp(i,:).*timeWindow;
    if latency(i) ~= 0
      ratecomp(i,time>=tOn(i)) = ratecomp(i,time>=tOn(i)).*(1-exp(-(time(time>=tOn(i))-tOn(i))/latency(i)));
      ratecomp(i,time>=tOff(i)) = ratecomp(i,time>=tOff(i)).*exp(-(time(time>=tOff(i))-tOff(i))/latency(i));
    end
  end
  x_N = ((1:N)-0.5)/N;
  %x_c = 0.25; % center of first category
  conn = 1./(1+exp(-cos(2*pi*(x_N-x_c(i)))/widthSigma));
  dynrate = dynrate + conn'*sum(ratecomp,1);
  S = nonhomPoissonGeneratorSpikeTimes(S_ini,dynrate,tau,kick,N,interval,dt);
end
