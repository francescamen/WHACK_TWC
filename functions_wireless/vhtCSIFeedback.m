function [mat, chanEstMinusCSD] = vhtCSIFeedback(rxSig, cfgVHT, userNum, numSTSVec)
% Computes the channel state information (CSI) feedback matrix per user, 
% for each subcarrier.

%   Copyright 2016 The MathWorks, Inc.

chanBW = cfgVHT.ChannelBandwidth;
numSTSTotal = sum(cfgVHT.NumSpaceTimeStreams); % == numTx for NDP
ind = wlanFieldIndices(cfgVHT);

% Estimate channel
rxVHTLTF  = rxSig(ind.VHTLTF(1):ind.VHTLTF(2),:);
demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF, chanBW, numSTSTotal);
chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF, chanBW, numSTSTotal); % Nst-by-Nsts-by-Nr

% Remove cyclic shift effect from the channel estimate
chanEstMinusCSD = vhtBeamformingRemoveCSD(chanEst, chanBW, numSTSTotal);  % Nst-by-Nsts-by-Nr
chanEstPerm = permute(chanEstMinusCSD, [3 2 1]); % Nsts-by-Nr-by-Nst

% Compute the feedback matrix using singular value decomposition
% for the streams allocated to the user
V = complex(zeros(length(chanEst), numSTSTotal, numSTSVec(userNum))); % Nst-by-Nr-by-Nsts
for i = 1:length(chanEst)
    [~, ~, V(i,:,:)] = svd(chanEstPerm(:,:,i), 'econ');
end
mat = V;

end

% [EOF]
