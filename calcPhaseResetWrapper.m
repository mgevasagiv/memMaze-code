%We examined the phenomenon of phase reset for slow theta and fast theta.
% This analysis was performed at the electrode level incorporating subject as a random effect (Rizzuto et al., 2003, 2006; Fell et al., 2004).
% EEG signals were filtered using Butterworth bandpass filters of 2–5 and 5–9 Hz.
% Phase values were then calculated during the study phase with the Hilbert transform at every time point (2 ms separation)
% the mean phase difference between study and test was compared across electrodes using the Rayleigh test.
% The comparison was conducted separately for the anterior hippocampus and posterior hippocampus.
function calcPhaseResetWrapper()

lowLimit = 4;
highLimit = 10;
samplingRate = 1000;

dataFilteredSW = bandpass(probeData, [lowLimit, highLimit], samplingRate);
nanInds = isnan(dataFilteredSW);
dataFilteredSW(nanInds) = 0;
dataFilteredSW(~isSleep) = 0;
dataFilteredSW(ISIpeakTimes - obj.minT_spike:ISIpeakTimes + obj.minT_spike) = 0;
currentPhases = angle(hilbert(dataFilteredSW-mean(dataFilteredSW))); % valuse range: [-pi:pi]


pacc = PACCalculator;
pacc.findPreferredPhase = false;
pacc.couplingIndexType = 'MI';
pacc.removeOutliers = false;

%SW-Sp
pacc.lowLimitPhase = obj.lowLimitSW;
pacc.highLimitPhase = obj.highLimitSW;
pacc.lowLimitAmp = obj.lowLimitSpindle;
pacc.highLimitAmp = obj.highLimitSpindle;
MISWSp = pacc.calculatePAC(dataSW,[],IISTimes{2});

end

% Helper func
function BP = bandpass(obj, timecourse, lowLimit, highLimit, filterOrder)

%bandpass code - from Maya
if (nargin < 5)
    filterOrder = obj.defaultFilterOrder;
end

% Maya GS - handle NAN values
indices = find(isnan(timecourse));
if length(indices) > obj.nanWarning*length(timecourse)
    warning('many NaN values in filtered signal')
end
timecourse(indices) = 0;
%

[b, a] = butter(filterOrder, [(lowLimit/obj.samplingRate)*2 (highLimit/obj.samplingRate)*2]);
BP = filtfilt(b, a, timecourse );
BP(indices) = NaN;
end