function y = vhtBeamformingRemoveCSD(x, chanBW, numSTSTotal)
% vhtBeamformingRemoveCSD Remove effect of transmitter cyclic shifts
%
%   Y = vhtBeamformingRemoveCSD(X, CHANBW, NUMSTSTOTAL) returns a channel
%   estimate matrix without the effect of cyclic shifts applied at the
%   transmitter.
%
%   Y is a complex Nst-by-Nsts-by-Nr array containing the estimated channel
%   at data and pilot subcarriers without the effects of transmitter cyclic
%   shifts. Nst is the number of occupied subcarriers, Nsts is the total
%   number of space-time streams and Nr is the number of receive antennas.
%
%   X is a complex Nst-by-Nsts-by-Nr array containing the estimated channel
%   at data and pilot subcarriers with the effects of transmitter cyclic
%   shifts.

%   Copyright 2015-2017 The MathWorks, Inc.

%#codegen

csd = wlan.internal.getCyclicShiftVal('VHT', numSTSTotal, ...
    wlan.internal.cbwStr2Num(chanBW));
cfgOFDM = wlan.internal.wlanGetOFDMConfig(chanBW, 'Long', 'VHT', numSTSTotal);
Nfft = cfgOFDM.FFTLength;
% Active frequency indices
k = sort([cfgOFDM.DataIndices; cfgOFDM.PilotIndices])-cfgOFDM.FFTLength/2-1;
% Negate CSD to remove cyclic shift
y = wlan.internal.cyclicShiftChannelEstimate(x, -csd, Nfft, k);

end
