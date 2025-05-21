function [y, CSI, GMatrix] = wlanEqualizeFM(x, chanEst, eqMethod, varargin)
%wlanEqualize Perform MIMO channel equalization. 
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   [Y, CSI] = wlanEqualize(X, CHANEST, 'ZF') performs equalization using
%   the signal input X and the channel estimation input CHANEST, and
%   returns the estimation of transmitted signal in Y and the soft channel
%   state information in CSI. The zero-forcing (ZF) method is used. The
%   inputs X and CHANEST can be double precision 2-D matrices or 3-D arrays
%   with real or complex values. X is of size Nsd x Nsym x Nr, where Nsd
%   represents the number of data subcarriers (frequency domain), Nsym
%   represents the number of OFDM symbols (time domain), and Nr represents
%   the number of receive antennas (spatial domain). CHANEST is of size Nsd
%   x Nsts x Nr, where Nsts represents the number of space-time streams.
%   The double precision output Y is of size Nsd x Nsym x Nsts. Y is
%   complex when either X or CHANEST is complex and is real otherwise. The
%   double precision, real output CSI is of size Nsd x Nsts.
%
%   [Y, CSI] = wlanEqualize(X, CHANEST, 'MMSE', NOISEVAR) performs the
%   equalization using the minimum-mean-square-error (MMSE) method. The
%   noise variance input NOISEVAR is a double precision, real, nonnegative
%   scalar.
%
%   See also wlanSTBCCombine.

%   Copyright 2015-2020 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% Input validation
narginchk(3, 5);

validateattributes(x, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanEqualize:InSignal', 'signal input');
validateattributes(chanEst, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanEqualize:ChanEst', 'channel estimation input');   
coder.internal.errorIf(~strcmp(eqMethod, 'ZF') && ~strcmp(eqMethod, 'MMSE'), ...
    'wlan:wlanEqualize:InvalidEqMethod');
coder.internal.errorIf(size(x, 1) ~= size(chanEst, 1), ...
    'wlan:wlanEqualize:UnequalFreqCarriers');
coder.internal.errorIf(size(x, 3) ~= size(chanEst, 3), ...
    'wlan:wlanEqualize:UnequalNumRx');

if strcmp(eqMethod, 'MMSE')
    narginchk(5,5);
    validateattributes(varargin{1}, {'double'}, {'real','scalar','nonnegative','finite','nonempty'}, ...
        'wlanEqualizer:noiseVarEst', 'noise variance estimation input'); 
    noiseVarEst = varargin{1};
    stsIdx = varargin{2};
else % ZF
    noiseVarEst = 0;
end

% Perform equalization
[numSc, numTx, numRx] = size(chanEst);

CSI = zeros(size(x, 1), numRx); % Pre-allocation here for code generation
numSym = size(x, 2);
% MU-MIMO
xTmp = permute(x, [3 2 1]); % numRx-by-numSym-by-numSc
chanEstTmp = permute(chanEst,[3 2 1]); % numRx-by-numTx(Nss tot)-by-numSc
y = coder.nullcopy(complex(zeros(numSc, numSym, size(stsIdx, 2))));
GMatrix = complex(zeros(numSc, size(stsIdx, 2), numRx));
for idx = 1:numSc
    H = chanEstTmp(:,:,idx);
    invH = inv(H*H'+noiseVarEst*eye(numRx));
    H_usr = H(:, stsIdx);
    G = H_usr' * invH;
    GMatrix(idx, :, :) = G;
    y(idx, 1:numSym, 1:size(stsIdx, 2)) = (G * xTmp(:, :, idx)).'; %#ok<MINV>
    
    invH_usr = inv(H_usr*H_usr'+noiseVarEst*eye(numRx));
    CSI(idx, :)  = 1./real(diag(invH_usr));
end

end