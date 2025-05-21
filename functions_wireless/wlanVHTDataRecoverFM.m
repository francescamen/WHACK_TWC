function [bits, crcBits, eqDataSym, GMatrix, varargout] = wlanVHTDataRecoverFM( ...
    rxVHTData, chanEst, noiseVarEst, cfgVHT, varargin)
%wlanVHTDataRecover Recover bits from VHT Data field signal
%
%   [BITS, CRCBITS] = wlanVHTDataRecover(RXVHTDATA, CHANEST, NOISEVAREST,
%   CFGVHTSU) recovers the bits in the VHT-Data field for a VHT format
%   single-user transmission.
%
%   BITS is an int8 column vector of length 8*CFGVHT.PSDULength containing
%   the recovered information bits.
%
%   CRCBITS is an int8 column vector of length 8 containing the VHT-Data
%   field checksum bits.
%
%   RXVHTDATA is the received time-domain VHT Data field signal, specified
%   as an Ns-by-Nr matrix of real or complex values. Ns represents the
%   number of time-domain samples in the VHT Data field and Nr represents
%   the number of receive antennas. Ns can be greater than the VHT Data
%   field length; in this case additional samples at the end of RXVHTDATA
%   are not used.
%
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the VHT-LTF. It is an array of size Nst-by-Nsts-by-Nr, where Nst
%   represents the total number of occupied subcarriers, Nsts represents
%   the total number of space-time streams used for the transmission and Nr
%   is the number of receive antennas.
%
%   NOISEVAREST is the noise variance estimate. It is a nonnegative scalar.
%
%   CFGVHTSU is the format configuration object of type <a 
%   href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, which
%   specifies the parameters for the single-user VHT format.
%
%   [BITS, CRCBITS] = wlanVHTDataRecover(RXVHTDATA, CHANEST, NOISEVAREST,
%   CFGVHTMU, USERNUMBER) recovers the bits in the VHT-Data field of a VHT
%   format multi-user transmission for an individual user of interest.
%
%   CFGVHTMU is the VHT format configuration for a multi-user transmission,
%   specified as a <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> object.
%
%   USERNUMBER is the user of interest, specified as an integer between 1
%   and NumUsers, where NumUsers is the number of users in the
%   transmission.
%
%   [BITS, CRCBITS] = wlanVHTDataRecover(RXVHTDATA, CHANEST, NOISEVAREST,
%   CFGVHTSU, USERNUMBER, NUMSTS) recovers the bits in the VHT-Data field
%   of a VHT format multi-user transmission for an individual user of
%   interest.
%
%   CFGVHTSU is the VHT format configuration for the user of interest,
%   specified as a <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> object.
%
%   NUMSTS is the number of space-time streams, specified as a
%   1-by-NumUsers vector. Element values specify the number of space-time
%   streams per user.
%
%   [BITS, CRCBITS] = wlanVHTDataRecover(..., NAME, VALUE) specifies
%   additional name-value pair arguments described below. When a name-value
%   pair is not specified, its default value is used.
%
%   'OFDMSymbolOffset'          OFDM symbol sampling offset. Specify the
%                               OFDMSymbolOffset as a fraction of the
%                               cyclic prefix (CP) length for every OFDM
%                               symbol, as a double precision, real scalar
%                               between 0 and 1, inclusive. The OFDM
%                               demodulation is performed based on Nfft
%                               samples following the offset position,
%                               where Nfft denotes the FFT length. The
%                               default value of this property is 0.75,
%                               which means the offset is three quarters of
%                               the CP length.
%
%   'EqualizationMethod'        Specify the equalization method as one of
%                               'MMSE' | 'ZF'. 'MMSE' indicates that the
%                               receiver uses a minimum mean square error
%                               equalizer. 'ZF' indicates that the receiver
%                               uses a zero-forcing equalizer. The default
%                               value of this property is 'MMSE'.
%
%   'PilotPhaseTracking'        Specify the pilot phase tracking performed
%                               as one of 'PreEQ' | 'None'. 'PreEQ' pilot
%                               phase tracking estimates and corrects a
%                               common phase offset across all subcarriers
%                               and receive antennas for each received OFDM
%                               symbol before equalization. 'None'
%                               indicates that pilot phase tracking does
%                               not occur. The default is 'PreEQ'.
%
%   'LDPCDecodingMethod'        Specify the LDPC decoding algorithm as one
%                               of these values:
%                               - 'bp'            : Belief propagation (BP)
%                               - 'layered-bp'    : Layered BP
%                               - 'norm-min-sum'  : Normalized min-sum
%                               - 'offset-min-sum': Offset min-sum
%                               The default is 'bp'.
%
%   'MinSumScalingFactor'       Specify the scaling factor for normalized
%                               min-sum LDPC decoding algorithm as a scalar
%                               in the interval (0,1]. This argument
%                               applies only when you set
%                               LDPCDecodingMethod to 'norm-min-sum'. The
%                               default is 0.75.
%
%   'MinSumOffset'              Specify the offset for offset min-sum LDPC
%                               decoding algorithm as a finite real scalar
%                               greater than or equal to 0. This argument
%                               applies only when you set
%                               LDPCDecodingMethod to 'offset-min-sum'. The
%                               default is 0.5.
%
%   'MaximumLDPCIterationCount' Specify the maximum number of iterations in
%                               LDPC decoding as a positive scalar integer.
%                               This applies when you set the channel
%                               coding property of format configuration
%                               object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> to 'LDPC'. 
%                               The default is 12.
%
%   'EarlyTermination'          To enable early termination of LDPC
%                               decoding, set this property to true. Early
%                               termination applies if all parity-checks
%                               are satisfied before reaching the number of
%                               iterations specified in the
%                               'MaximumLDPCIterationCount' input. To let
%                               the decoding process iterate for the number
%                               of iterations specified in the
%                               'MaximumLDPCIterationCount' input, set this
%                               argument to false. This applies when you
%                               set the channel coding property of format
%                               configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>
%                               to 'LDPC'.The default is false.
%
%   [..., EQDATASYM, CPE] = wlanVHTDataRecover(...) also returns the
%   equalized subcarriers and common phase error.
%
%   EQDATASYM is a complex Nsd-by-Nsym-by-Nss array containing the
%   equalized symbols at data carrying subcarriers. Nsd represents the
%   number of data subcarriers, Nsym represents the number of OFDM symbols
%   in the VHT-Data field, and Nss represents the number of spatial
%   streams assigned to the user.
%
%   CPE is a column vector of length Nsym containing the common phase error
%   between each received and expected OFDM symbol.
%
%   Example: 
%   %  Recover bits in VHT Data field via channel estimation on VHT-LTF 
%   %  over a 2 x 2 quasi-static fading channel
%
%     % Configure a VHT configuration object 
%     chanBW = 'CBW160';
%     cfgVHT = wlanVHTConfig('ChannelBandwidth',    chanBW, ...
%         'NumTransmitAntennas', 2, 'NumSpaceTimeStreams', 2, ...
%         'APEPLength',          512); 
%
%     % Generate VHT-LTF and VHT Data field signals
%     txDataBits = randi([0 1], 8*cfgVHT.PSDULength, 1);
%     txVHTLTF  = wlanVHTLTF(cfgVHT); 
%     txVHTData = wlanVHTData(txDataBits, cfgVHT);
%
%     % Pass through a 2 x 2 quasi-static fading channel with AWGN 
%     H = 1/sqrt(2)*complex(randn(2, 2), randn(2, 2));
%     rxVHTLTF  = awgn(txVHTLTF  * H, 10);
%     rxVHTData = awgn(txVHTData * H, 10);
%
%     % Perform channel estimation based on VHT-LTF
%     demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF, cfgVHT, 1);
%     chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF, cfgVHT);
%
%     % Recover information bits in VHT Data
%     rxDataBits = wlanVHTDataRecover(rxVHTData, chanEst, 0.1, ...
%         cfgVHT, 'EqualizationMethod', 'ZF');
%
%     % Compare against original information bits
%     disp(isequal(txDataBits, rxDataBits));
%
%   See also wlanVHTConfig, wlanVHTData, wlanVHTLTF, wlanVHTLTFDemodulate,
%   wlanVHTLTFChannelEstimate.

%   Copyright 2015-2019 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

narginchk(4, 22);
nargoutchk(0, 4);

% Calculate CPE if requested
if nargout>3
    calculateCPE = true;
else
    calculateCPE = false;
end

% VHT configuration input self-validation
validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, 'VHT format configuration object');

% Validate and parse optional inputs
[muSpec, userNum, numSTSVec, recParams] = wlan.internal.parseVHTOptionalInputs(mfilename, cfgVHT.NumSpaceTimeStreams, varargin{:});

% Check MU
if muSpec==2 % SU CFGVHT
    % Have a SU cfgVHT object as input
    propIdx = 1;
    numSTSu = numSTSVec(userNum);
elseif muSpec==1 % MU CFGVHT
    validateattributes(varargin{1}, {'numeric'}, {'real','integer','scalar','>=',1,'<=',cfgVHT.NumUsers}, mfilename, 'USERNUMBER');
    % Have a MU cfgVHT object as input
    propIdx = userNum;
    numSTSu = numSTSVec(propIdx);
else % not specified, set defaults
    % Single-user case
    propIdx = 1;
    numSTSu = numSTSVec(propIdx);
end

cfgInfo = validateConfig(cfgVHT, 'MCS');
mcsTable = wlan.internal.getRateTable(cfgVHT);
chanBW = cfgVHT.ChannelBandwidth;

% NDP only for SU, so idx is (1)
if cfgVHT.APEPLength(1) == 0 
    bits     = zeros(0, 1, 'int8');
    crcBits  = zeros(0, 1, 'int8');
    eqDataSym = zeros(mcsTable.NSD(1), 0, mcsTable.Nss(1));
    if calculateCPE==true
        varargout{1} = []; % CPE
    end
    return;
end

% All optional params: parsed and validated
numSTSTotal = sum(numSTSVec);

% Signal input self-validation
validateattributes(rxVHTData, {'double'}, {'2d','finite'}, 'rxVHTData', 'VHT-Data field signal'); 
validateattributes(chanEst, {'double'}, {'3d','finite'}, 'chanEst', 'channel estimation'); 
validateattributes(noiseVarEst, {'double'}, {'real','scalar','nonnegative','finite'}, 'noiseVarEst', 'noise variance estimation'); 

% Set up some implicit configuration parameters
numBPSCS   = mcsTable.NBPSCS(propIdx);    % Number of coded bits per single carrier
numCBPS    = mcsTable.NCBPS(propIdx);     % Number of coded bits per OFDM symbol
numDBPS    = mcsTable.NDBPS(propIdx);
rate       = mcsTable.Rate(propIdx);
numES      = mcsTable.NES(propIdx);       % Number of encoded streams
numSS      = mcsTable.Nss(propIdx);       % Number of spatial streams
numSeg     = strcmp(chanBW, 'CBW160') + 1;
% Number of coded bits per OFDM symbol, per spatial stream, per segment
numCBPSSI  = numCBPS/numSS/numSeg;  
numRx      = size(rxVHTData, 2);

% Get OFDM configuration
[cfgOFDM, dataInd, pilotInd] = wlan.internal.wlanGetOFDMConfig(chanBW, cfgVHT.GuardInterval, 'VHT', numSTSTotal);
    
% Set channel coding
coder.varsize('channelCoding',[1,4]);
channelCoding = getChannelCoding(cfgVHT);

% Cross-validation between inputs
if muSpec==2
    coder.internal.errorIf(cfgVHT.NumSpaceTimeStreams(1) ~= numSTSu, ...
    'wlan:wlanVHTDataRecover:InvalidNumSTS', numSTSu, cfgVHT.NumSpaceTimeStreams(1));
end
coder.internal.errorIf(cfgVHT.STBC && muSpec > 0, 'wlan:wlanVHTDataRecover:InvalidSTBCMU');

numST = numel([dataInd; pilotInd]); % Total number of occupied subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST, 'wlan:wlanVHTDataRecover:InvalidChanEst1D', numST);
coder.internal.errorIf(size(chanEst, 2) ~= numSTSTotal, 'wlan:wlanVHTDataRecover:InvalidChanEst2D', numSTSTotal);
coder.internal.errorIf(size(chanEst, 3) ~= numRx, 'wlan:wlanVHTDataRecover:InvalidChanEst3D');

% Extract data and pilot subcarriers from channel estimate
chanEstData = chanEst(dataInd,:,:);
chanEstPilots = chanEst(pilotInd,:,:);

% Cross-validation between inputs
numOFDMSym = cfgInfo.NumDataSymbols;
symLen = cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength;
minInputLen = numOFDMSym*(symLen);
coder.internal.errorIf(size(rxVHTData, 1) < minInputLen, ...
    'wlan:wlanVHTDataRecover:ShortDataInput', minInputLen);

% Calculate the number of whole OFDM symbols in the input signal. This
% accounts for the MU-padding and LDPC extra symbol, if needed.
if muSpec==2 && any(strcmp(channelCoding,'LDPC'))
    numOFDMSym = floor(length(rxVHTData)/(symLen));
    % Recalculate the minimum input signal length required to process LDPC
    % encoded data with an extra LDPC encoded symbol.
    minInputLen = numOFDMSym*(symLen);
end

% OFDM demodulation
[ofdmDemodData, ofdmDemodPilots] = wlan.internal.wlanOFDMDemodulate(rxVHTData(1:minInputLen, :), cfgOFDM, recParams.OFDMSymbolOffset);

% Index into streams for the user of interest
stsIdx = sum(numSTSVec(1:(userNum-1)))+(1:numSTSu); 

% Pilot phase tracking
if calculateCPE==true || strcmp(recParams.PilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from Eqn 22-95, IEEE Std 802.11ac-2013
    % Offset by 4 to allow for L-SIG, VHT-SIG-A, VHT-SIG-B pilot symbols
    n = (0:numOFDMSym-1).';
    z = 4;
    refPilots = wlan.internal.vhtPilots(n, z, chanBW, sum(numSTSVec));
    
    % Estimate CPE and phase correct symbols
    cpe = wlan.internal.commonPhaseErrorEstimate(ofdmDemodPilots, chanEstPilots(:,stsIdx,:), refPilots(:,:,stsIdx));
    if strcmp(recParams.PilotPhaseTracking, 'PreEQ')
        ofdmDemodData = wlan.internal.commonPhaseErrorCorrect(ofdmDemodData, cpe);
    end
    if calculateCPE==true
        varargout{1} = cpe.'; % Permute to Nsym-by-1
    end
end

% Equalization
if cfgVHT.STBC  % Only SU
    [eqDataSym, dataCSI] = wlan.internal.wlanSTBCCombine(ofdmDemodData, chanEstData, numSS, recParams.EqualizationMethod, noiseVarEst);
else    % Both SU and MU
    [eqDataSym, dataCSI, GMatrix] = wlanEqualizeFM(ofdmDemodData, chanEstData, recParams.EqualizationMethod, noiseVarEst, stsIdx);
end

% Segment parsing of symbols 
parserOut  = wlanSegmentParseSymbols(eqDataSym, chanBW); % [Nsd/Nseg Nsym Nss Nseg]
csiParserOut = wlanSegmentParseSymbols(reshape(dataCSI, [], 1, numSS), chanBW); % [Nsd/Nseg 1 Nss Nseg]

% LDPC Tone demapping
if strcmp(channelCoding{propIdx},'LDPC')
    mappingIndicesLDPC = wlan.internal.getToneMappingIndices(chanBW);
    parserOut = parserOut(mappingIndicesLDPC,:,:,:);
    csiParserOut = csiParserOut(mappingIndicesLDPC,:,:,:);
end
 
% Constellation demapping
qamDemodOut = wlanConstellationDemap(parserOut, noiseVarEst, numBPSCS); % [Ncbpssi,Nsym,Nss,Nseg]

% Apply bit-wise CSI and concatenate OFDM symbols in the first dimension
qamDemodOut = bsxfun(@times, ...
        reshape(qamDemodOut, numBPSCS, [], numOFDMSym, numSS, numSeg), ...
        reshape(csiParserOut, 1, [], 1, numSS, numSeg));
qamDemodOut = reshape(qamDemodOut, [], numSS, numSeg); % [(Ncbpssi*Nsym),Nss,Nseg]

% BCC Deinterleaving
if strcmp(channelCoding{propIdx}, 'BCC')
    deintlvrOut = wlanBCCDeinterleave(qamDemodOut, 'VHT', numCBPSSI, chanBW); % [(Ncbpssi*Nsym),Nss,Nseg]
else
    % Deinterleaving is not required for LDPC
    deintlvrOut = qamDemodOut;
end

% Segment deparsing of bits
segDeparserOut = wlanSegmentDeparseBits(deintlvrOut, chanBW, numES, numCBPS, numBPSCS); % [(Ncbpss*Nsym),Nss]

% Stream deparsing
streamDeparserOut = wlanStreamDeparse(segDeparserOut(:,:), numES, numCBPS, numBPSCS); % [(Ncbps*Nsym/Nes),Nes]
% Indexing for codegen

if strcmp(channelCoding{propIdx}, 'BCC')
    % Channel decoding for BCC
    numTailBits = 6;
    chanDecOutPreDeparse = wlanBCCDecode(streamDeparserOut, mcsTable.Rate(propIdx));
    % BCC decoder deparser
    chanDecOut = reshape(chanDecOutPreDeparse(1:end-numTailBits,:)', [], 1);
else
   % Channel decoding for LDPC
   % Calculate numSymMaxInit as specified in IEEE Std 802.11ac-2013,
   % Section 22.3.10.5.4, Eq 22-65 and Section 22.3.21, Eq 22-107
   numSym = cfgInfo.NumDataSymbols(1);
   
   % Estimate number of OFDM symbols as specified in IEEE Std
   % 802.11ac-2013, Section 22.3.21, Eq 22-107.
   mSTBC = (cfgVHT.NumUsers == 1)*(cfgVHT.STBC ~= 0) + 1;
   numSymMaxInit = numSym - mSTBC*cfgInfo.ExtraLDPCSymbol;
   
   % Compute the number of payload bits as specified in IEEE Std
   % 802.11ac-2013, Section 22.3.10.5.4, Eq 22-61 and Eq 22-66.
   numPLD = numSymMaxInit*numDBPS;
   
   % LDPC decoding parameters as per IEEE Std 802.11-2012, Section
   % 20.3.11.17.4 and IEEE Std 802.11ac-2013, Section 22.3.10.5.4.
   cfg = wlan.internal.getLDPCparameters(numDBPS, rate, mSTBC, numPLD, numOFDMSym);
   chanDecOut = wlan.internal.wlanLDPCDecode(streamDeparserOut, cfg, ...
        recParams.algChoice, recParams.alphaBeta, recParams.MaximumLDPCIterationCount, recParams.EarlyTermination);
end

% Derive initial state of the scrambler 
scramInit = wlan.internal.scramblerInitialState(chanDecOut(1:7));

% Remove pad and tail bits, and descramble
if all(scramInit==0)
    % Scrambler initialization invalid (0), therefore do not descramble
    descramBits = chanDecOut(1:16+8*cfgVHT.PSDULength(propIdx));
else
    descramBits = wlanScramble(chanDecOut(1:16+8*cfgVHT.PSDULength(propIdx)), scramInit);
end

% Outputs
crcBits = descramBits(9:16);
bits = descramBits(17:end);

end