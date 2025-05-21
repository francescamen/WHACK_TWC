function est = helperNoiseEstimate(rxSym,varargin)
%helperNoiseEstimate Estimate noise power using L-LTF (non-HT, HT, and VHT) 
%and LTF1 in S1G
%
%   EST = helperNoiseEstimate(RXSYM) estimates the mean noise power in
%   watts using the demodulated L-LTF symbols (non-HT, HT, and VHT) or the
%   demodulated S1G-LTF1 symbols in S1G, assuming 1ohm resistance. The
%   estimated noise power in non-HT packet is averaged over the number of
%   of receive antennas.
%
%   RXSYM is the frequency-domain signal corresponding to the L-LTF or
%   S1G-LTF1. It is a complex matrix or 3-D array of size Nst-by-2-by-Nr,
%   where Nst represents the number of used subcarriers in the L-LTF or
%   S1G-LTF1, and Nr represents the number of receive antennas. Two OFDM
%   symbols in the L-LTF or S1G-LTF1 fields are used to estimate the noise
%   power. Noise estimate from the S1G-LTF1 field for S1G 1MHz format is
%   not supported.
%
%   EST = helperNoiseEstimate(RXSYM,CHANBW,NUMSTS) returns the estimated
%   noise power in beamformed fields using the specified channel bandwidth,
%   CHANBW, and total number of space-time streams, NUMSTS. The number of
%   subcarriers used within each field and the scaling applied during
%   demodulation differs between the non-HT and HT/VHT fields for VHT, HT
%   and non-HT. Therefore the estimated noise power after demodulation in
%   HT/VHT fields is calculated by scaling the estimated noise power in the
%   L-LTF. The number of space-time streams is required for scaling if the
%   demodulated RXSYM is not scaled according to the number of space-time
%   streams.
%
%   EST = helperNoiseEstimate(...,'Per Antenna') specifies the option of
%   estimating the noise for each receive antenna. When this option is
%   specified EST is a row vector of length Nr.
%
%   Example:
%   %  Estimate the noise variance of an HT packet.
%
%      cfgHT = wlanHTConfig; % Create packet configuration
%      chanBW = cfgHT.ChannelBandwidth;
%      numSTS = cfgHT.NumSpaceTimeStreams;
%      noisePower = -20;
%      awgnChannel = comm.AWGNChannel;
%      awgnChannel.NoiseMethod = 'Variance';
%      awgnChannel.Variance = 10^(noisePower/10);
%
%      Nst = 56;  % Data and pilot OFDM subcarriers in 20MHz, HT format
%      Nfft = 64; % FFT size for 20MHz bandwidth
%      nVarHT = 10^(noisePower/10)*(Nst/Nfft); % non-HT noise variance
%
%      NumRxAnts = 1;
%      % Average noise estimate over 100 independent noise realization
%      for n=1:100
%         % Generate LLTF and add noise
%         rxSym = awgnChannel(wlanLLTF(cfgHT));
%         y = wlanLLTFDemodulate(rxSym,cfgHT);
%         noiseEst(n) = helperNoiseEstimate(y,chanBW,numSTS);
%      end
%
%      % Check noise variance estimates without Channel
%      noiseEstError = 10*log10(mean(noiseEst))-10*log10(nVarHT);
%      disp(['Error between noise variance and mean estimated noise ', ...
%      'power(dB): ' num2str(noiseEstError,'%2.2f ')]);
%
%   See also wlanLLTF, wlanLLTFDemodulate.

%   Copyright 2015-2017 The MathWorks, Inc.

%#codegen

narginchk(1,4);

% Validate symbol type
validateattributes(rxSym,{'double'},{'3d','finite'},mfilename,'OFDM symbol(s)');

% Two L-LTF symbols (non-HT, HT, and VHT) or two S1G-LTF1 symbols (S1G) are
% required to estimate the noise
coder.internal.errorIf(size(rxSym,2)~=2,'wlan:helperNoiseEstimate:IncorrectNumSyms');

numSC = size(rxSym,1);

% Minimal optional parameter checks
if nargin == 2
    % (rxSym,'Per Antenna')
    validateNType(varargin{1});
    average = false;

    scalingFactor = 1; % Noise scaling factor
elseif nargin == 3
    % (rxSym, chanBW, numSTSTotal)
    chanBW = varargin{1};        % chanBW: CBW2/4/8/16/20/40/80/160
    coder.internal.errorIf(strcmp(chanBW,'CBW1'),'wlan:helperNoiseEstimate:InvalidS1G1M');
    validateInput(chanBW,numSC);

    numSTSTotal = varargin{2};   % numSTSTotal: 1,...,8 
    average = true;
    
    % Noise scaling factor
    scalingFactor = noiseScaling(chanBW,numSC,numSTSTotal);
elseif nargin == 4
    % (rxSym, chanBW, numSTSTotal, 'Per Antenna')
    chanBW = varargin{1};       % chanBW: CBW2/4/8/16/20/40/80/160
    coder.internal.errorIf(strcmp(chanBW,'CBW1'),'wlan:helperNoiseEstimate:InvalidS1G1M');
    validateInput(chanBW,numSC);

    numSTSTotal = varargin{2};  % numSTSTotal: 1,...,8 

    validateNType(varargin{3});
    average = false;

    % Noise scaling factor
    scalingFactor = noiseScaling(chanBW,numSC,numSTSTotal);    
else
    % (rxSym)
    average = true;
    scalingFactor = 1;
end

% Noise estimate
noiseEst = sum(abs(rxSym(:,1,:)-rxSym(:,2,:)).^2,1)/(2*numSC);
if average
    noise = mean(noiseEst);
else
    noise = squeeze(noiseEst).';
end

% Scale
est = noise*scalingFactor;

end

%-------------------------------------------------------------------------
function out = noiseScaling(chanBW,numSC,numSTSTotal)

if any(strcmp(chanBW,{'CBW2','CBW4','CBW8','CBW16'})) 
    % In S1G, Data and LTF1 fields have the same number of occupied
    % subcarriers. Only apply scaling by the number of space-time streams.
    Nst = numSC;
else
    % Get the number of occupied subcarriers in HT and VHT fields.
    %   The number of used subcarriers for HT and VHT are same therefore
    %   fix the character vector input of the following helper function to 
    %   VHT. The guard type is not relevant for numbers alone.
    [~,vhtData,vhtPilots] = wlan.internal.wlanGetOFDMConfig(chanBW,'Long','VHT');
    Nst = numel(vhtData)+numel(vhtPilots);
end
out = (Nst/numSC)*numSTSTotal;

end

%-------------------------------------------------------------------------
function validateInput(chanBW,numSC)

if any(strcmp(chanBW,{'CBW2','CBW4','CBW8','CBW16'})) % S1G
    [~,s1gData,s1gPilots] = wlan.internal.s1gOFDMConfig(chanBW,'Long','LTF1');
    Nst = numel(s1gData)+numel(s1gPilots);
else % nonHT, HT, and VHT
    % Get number of used subcarriers in NonHT format
    [~,nonhtData,nonhtPilots] = wlan.internal.wlanGetOFDMConfig(chanBW,'Long','Legacy');
    Nst = numel(nonhtData)+numel(nonhtPilots);
end

% Validate number of subcarriers in input
coder.internal.errorIf(numSC~=Nst,'wlan:helperNoiseEstimate:IncorrectNumSC',Nst,numSC);

end

%-------------------------------------------------------------------------
function validateNType(nType)

coder.internal.errorIf(~(strcmpi(nType,'Per Antenna')),'wlan:helperNoiseEstimate:InvalidNoiseEstType');

end
