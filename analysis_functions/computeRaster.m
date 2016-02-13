function raster  = computeRaster(t,V)
  raster = [];
  [indTimes,neuronSpikes] = find (V > 0);
  if ~isempty(neuronSpikes)
    tSpikes = t(indTimes); % in s
    raster(:,1) = tSpikes;
    raster(:,2) = neuronSpikes;
    raster(diff(raster(:,2))==0 & diff(raster(:,1))<=1.05*dsfact*dt,:) = []; % removing artificial spikes that come from two consecutive voltages above 0 mV
    [~,indSort] = sort(raster(:,1));
    raster = raster(indSort,:);
  end
