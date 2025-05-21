clear all
close all

addpath('../')

%% WHACK attack emulation
% This is a modification of the WINNERVHTMUMIMOExample.m provided by MathWorks, Inc.

% Copyright 2022 Francesca Meneghello

%% 802.11ac Multiuser MIMO Precoding with WINNER II Channel Model
% This example shows the transmit and receive processing for a 802.11ac(TM)
% multiuser downlink transmission over experimentally emulated channels

%% Simulation Parameters and Configuration
% For 802.11ac, a maximum of eight spatial streams is allowed. 
% Different rate parameters and payload sizes for each user are specified
% as vector parameters for the transmission configuration.

% s = rng(20);                           % Set RNG seed for repeatability

% Transmission parameters
chanBW      = 'CBW80';               % Channel bandwidth
numUsers    = 2;                     % Number of users
numSTSVec   = [1 1];               % Number of streams per user
userPos     = [0 1];               % User positions (inside the Group ID Management frame together with the Membership Status, i.e., GroupID inside the wlanVHTConfig)
mcsVec      = [4 4];               % MCS per user
apepVec     = [520 520];         % Payload per user, in bytes
chCodingVec = {'LDPC','LDPC'}; % Channel coding per user

% Precoding and equalization parameters
precodingType = 'ZF';                % Precoding type; ZF or MMSE
snr           = 25;                  % SNR in dB
eqMethod      = 'MMSE';                 % Equalization method

numTx = sum(numSTSVec);
% The number of transmit antennas is set to be the sum total of all the
% used space-time streams. This implies no space-time block coding (STBC)
% or spatial expansion is employed for the transmission.

% Create the multiuser VHT format configuration object
cfgVHTMU = wlanVHTConfig('ChannelBandwidth',chanBW, ...
    'NumUsers',numUsers, ...
    'NumTransmitAntennas',numTx, ...
    'GroupID',2, ...
    'NumSpaceTimeStreams',numSTSVec,...
    'UserPositions',userPos, ...
    'MCS',mcsVec, ...
    'APEPLength',apepVec, ...
    'ChannelCoding',chCodingVec);

%% Sounding (NDP) Configuration
% For precoding, channel sounding is first used to determine the channel
% experienced by the users (receivers). This channel state information is
% sent back to the transmitter, for it to be used for subsequent data
% transmission. It is assumed that the channel varies slowly over the two
% transmissions. For multiuser transmissions, the same NDP (Null Data
% Packet) is transmitted to each of the scheduled users [ <#18 2> ].

% VHT sounding (NDP) configuration, for same number of streams
cfgVHTNDP = wlanVHTConfig('ChannelBandwidth',chanBW, ...
    'NumUsers',1, ...
    'NumTransmitAntennas',numTx, ...
    'GroupID',0, ...
    'NumSpaceTimeStreams',sum(numSTSVec),...
    'MCS',0, ...
    'APEPLength',0);

% The number of streams specified is the sum total of all space-time
% streams used. This allows the complete channel to be sounded.

% Generate the null data packet, with no data
txNDPSig = wlanWaveformGenerator([],cfgVHTNDP);
NPDSigLen = size(txNDPSig, 1);
        
%% Define common variables for the loop
tx_ant_vect = linspace(0, numTx-1, numTx).';
numPadZeros = 10;

% Files of channel parameters of the users
base_name = '../mat_files/020822/020822_pos';

attacker_idx = 12;
victim_idx = 14;

attacker_idx = num2str(attacker_idx);

end_time = 2000;  
end_iter = 1;

%% Create the channel for each user and sound the channel
victim_name = num2str(victim_idx);

file_attacker = strcat(base_name, attacker_idx, '_nlos');
file_victim = strcat(base_name, victim_name, '_nlos');
folders_users = {file_victim, file_attacker};

power_limit = 1e5;
numSTSelected = 150;

name_save_folder = strcat('../emulation_results/020822/change_subcarriers_nlos_NCG/no_subc_', num2str(numSTSelected));

if isfolder(name_save_folder)
    name_save_max = strcat(name_save_folder, '/results_time', ...
        num2str(end_time), '_iteration', ...
        num2str(end_iter), '.mat');
    if isfile(name_save_max)
        return
    end
else
    mkdir(name_save_folder);
end

amplitudes_users = cell(numUsers, 1);
aoas_users = cell(numUsers, 1);
toas_users = cell(numUsers, 1);

for us=1:numUsers
    % Upload the experimental channel parameters        
    amplitudes_users{us} = importdata(strcat(folders_users{us}, '/paths_amplitude_list.m'));
    aoas_users{us} = importdata(strcat(folders_users{us}, '/paths_aoa_list.m'));
    toas_user = importdata(strcat(folders_users{us}, '/paths_toa_list.m'));
    for ii=1:size(toas_user, 2)
        toas_user{ii} = toas_user{ii} - toas_user{ii}(:, 1);
    end
    toas_users{us} = toas_user;
end

%% NCG compute the symbolic expressions
syms v1_i [numTx numSTSVec(2)] matrix
syms v2_i [numTx numSTSVec(2)] matrix

[functpowerST, gradPower] = power_constraint_gradient_symbolic(numSTSVec, v1_i, v2_i);
[functobjectiveST, gradObjective] = objective_function_gradient_symbolic_newobjdet(numSTSVec, v1_i, v2_i);

for time_idx=1:end_time
    disp(['time_idx ' num2str(time_idx)]);
    time_idx_expanded = sprintf('%05d', time_idx);
    for num_iter=1:end_iter

        num_iter_expanded = sprintf('%02d', num_iter);
        name_save_file = strcat(name_save_folder, '/results_time', ...
            num2str(time_idx_expanded), '_iteration', ...
            num2str(num_iter_expanded), '.mat');

        if isfile(name_save_file)
             continue
        end

        % disp(['num_iter ' num2str(num_iter)]);
        chanOutNDP = cell(numUsers, 1);
        chanDelay = zeros(numUsers, 1);
        chanFiltUsers = cell(numUsers, 1);
        amplitudesChannelUsers = cell(numUsers, 1);
    
        for us=1:numUsers
            % Upload the experimental channel parameters
                
            amplitudes = amplitudes_users{us}{time_idx};
            aoas = aoas_users{us}{time_idx};
            toas = toas_users{us}{time_idx};
            
            chanFilt = comm.ChannelFilter('SampleRate', 8e7, ...
                                          'PathDelays', toas);
                
            amplitudes_channel = zeros(1, size(toas, 2), numTx, numSTSVec(us));  %1-by-Nr-by-Nt-by-Np or Ns-by-Nr-by-Nt-by-Np
            amplitudes_channel(1, 1:size(amplitudes, 2), :, 1)  = amplitudes.' .* exp(1i*pi*cos(aoas).' *tx_ant_vect.');
        
            chanFiltUsers{us} = chanFilt;
            amplitudesChannelUsers{us} = amplitudes_channel;
            
            % Sound the emulated channel
            chanOutNDP{us} = chanFilt([txNDPSig;zeros(numPadZeros,numTx)], amplitudes_channel); 
            chanDelay(us) = 7;  % as in the simulation
        end
        % Add AWGN
        rxNDPSig = cellfun(@awgn,chanOutNDP, ...
            num2cell(snr*ones(numUsers,1)),'UniformOutput',false);
            
        %% Channel State Information Feedback
        % Each user estimates its own channel using the received NDP signal and
        % computes the channel state information that it can send back to the
        % transmitter. This example uses the singular value decomposition of the
        % channel seen by each user to compute the CSI feedback.
        
        mat = cell(numUsers,1);  % V matrix
        csi_feed = cell(numUsers,1); % chanEstMinusCSD: channel estimate without cyclic prefix effect
        for uIdx = 1:numUsers
            % Compute the feedback matrix based on received signal per user
            [mat{uIdx}, csi_feed{uIdx}] = vhtCSIFeedback(rxNDPSig{uIdx}(chanDelay(uIdx)+1:end,:), ...
                cfgVHTNDP,uIdx,numSTSVec);  % Nst-by-Nr-by-Nsts
        end
        numST = length(mat{1});         % Number of subcarriers
        
        %% Beamforming feedback quantization and reconstruction
        mat_rec = cell(numUsers,1);  % V tilde matrix reconstructed
        for uIdx = 1:numUsers
            V_mat = mat{uIdx};
            V_mat_tilde_rec = zeros(numST, numTx, numSTSVec(uIdx));
            
            for s_i = 1:numST
                V_subcarr = squeeze(V_mat(s_i, :, :));
                if length(size(V_mat)) == 2
                    V_subcarr = permute(V_subcarr, [2, 1]);
                end
                [Vtilde_rec, Vtildecompl] = beamformingFeedbackQuantization(V_subcarr, numTx, numSTSVec(uIdx));
                if length(size(V_mat)) == 2
                    Vtilde_rec = permute(Vtilde_rec, [2, 1]);
                end
                V_mat_tilde_rec(s_i, :, :) = Vtilde_rec;
            end
            mat_rec{uIdx} = V_mat_tilde_rec;
        end

        mat_rec_no_attack = mat_rec;
                
        %% Attack vector random and NCG to enforce constraints
        
        % Check function before the attack
        sumPowerOrig = power_constraint(mat_rec{2}, mat_rec{1}, numST, numSTSVec);
        [~, minObjOrig] = min_det_prob(mat_rec{2}, mat_rec{1}, numST);
        
        size_x = size(mat_rec{2});        

        done = false;
        rand_factor = 0.1;
        mat_rec_attacker_init = rand_factor*randn(size_x, "like", 1i);
        mat_rec_attacker_init(:, size_x(2)) = abs(real(mat_rec_attacker_init(:, size_x(2))));

        max_fails = 10;
        num_fails = 0;
        while ~ done
            try
                [mat_rec_attacker, num_iterations] = applyNCG(mat_rec_attacker_init, ...
                    size_x, mat_rec{2}, numST, numSTSVec, numSTSelected, mat_rec{1}, ...
                    power_limit, functpowerST, gradPower, functobjectiveST, gradObjective, v1_i, v2_i);
                done = true;
            catch ME
                disp(ME.identifier)
                num_fails = num_fails + 1;
                if num_fails >= max_fails
                    done = true;
                end
                mat_rec_attacker_init = rand_factor*randn(size_x, "like", 1i);
                mat_rec_attacker_init(:, size_x(2), :) = abs(real(mat_rec_attacker_init(:, size_x(2), :)));
                disp('notdone')
                continue
            end
        end
        if num_fails >= max_fails
            disp('discarding problematic channel')
            continue
        end
        
        % Apply attack
        mat_rec{2} = 1*squeeze(mat_rec_attacker);
        
        sumPowerAttack = power_constraint(mat_rec_attacker, mat_rec{1}, numST, numSTSVec);
        [~, minObjAttack] = min_det_prob(mat_rec_attacker, mat_rec{1}, numST); zeros(numST, 1);
                
        %% Precoding
        % The transmitter computes the steering matrix for the data
        % transmission using either Zero-Forcing or Minimum-Mean-Square-Error
        % (MMSE) based precoding techniques. Both methods attempt to cancel out the
        % intra-stream interference for the user of interest and interference due
        % to other users. The MMSE-based approach avoids the noise enhancement
        % inherent in the zero-forcing technique. As a result, it performs better
        % at low SNRs.
        
        % Pack the per user CSI into a matrix
        steeringMatrix = zeros(numST,sum(numSTSVec),sum(numSTSVec));
        %   Nst-by-Nt-by-Nsts
        for uIdx = 1:numUsers
            stsIdx = sum(numSTSVec(1:uIdx-1))+(1:numSTSVec(uIdx));
            steeringMatrix(:,:,stsIdx) = mat_rec{uIdx};     % Nst-by-Nt-by-Nsts
        end
        
        % Zero-forcing or MMSE precoding solution
        if strcmp(precodingType, 'ZF')
            delta = 0; % Zero-forcing
        else
            delta = (numTx/(10^(snr/10))) * eye(numTx); % MMSE
        end
        for i = 1:numST
            % Channel inversion precoding
            v = squeeze(steeringMatrix(i,:,:));
            steeringMatrix(i,:,:) = v/(v'*v + delta);
        end
        
        % Set the spatial mapping based on the steering matrix
        cfgVHTMU.SpatialMapping = 'Custom';
        cfgVHTMU.SpatialMappingMatrix = permute(steeringMatrix,[1 3 2]); % Nst-by-Nsts-by-Nt
        
        %% Data Transmission
        % Random bits are used as the payload for the individual users. A cell
        % array is used to hold the data bits for each user, |txDataBits|. For a
        % multiuser transmission the individual user payloads are padded such that
        % the transmission duration is the same for all users. This padding process
        % is described in Section 9.12.6 of [ <#18 1> ]. In this example for
        % simplicity the payload is padded with zeros to create a PSDU for each
        % user.
        
        % Create data sequences, one for each user
        txDataBits = cell(numUsers,1);
        psduDataBits = cell(numUsers,1);
        for uIdx = 1:numUsers
            % Generate payload for each user
            txDataBits{uIdx} = randi([0 1],cfgVHTMU.APEPLength(uIdx)*8,1,'int8');
            
            % Pad payload with zeros to form a PSDU
            psduDataBits{uIdx} = [txDataBits{uIdx}; ...
                zeros((cfgVHTMU.PSDULength(uIdx)-cfgVHTMU.APEPLength(uIdx))*8,1,'int8')];
        end
        
        %%
        % Using the format configuration, |cfgVHTMU|, with the steering matrix, to
        % generate the multiuser VHT waveform.
        txSig = wlanWaveformGenerator(psduDataBits,cfgVHTMU);
        
        %%
        % As we restart the channel, we want
        % it to re-process the NDP before the waveform so as to accurately mimic
        % the channel continuity. Only the waveform portion of the channel's output
        % is extracted for the subsequent processing of each user.
        
        % Transmit through the emulated channel for all users, with 10 all-zero
        % samples appended to account for channel filter delay
        
        chanOut = cell(numUsers, 1);
        for us=1:numUsers
            chanOut{us} = chanFiltUsers{us}([txNDPSig; zeros(numPadZeros,numTx); ...
                txSig; zeros(numPadZeros,numTx)], amplitudesChannelUsers{us}); 
        end
        
        % Extract the waveform output for each user
        chanOut = cellfun(@(x) x(NPDSigLen+numPadZeros+1:end,:),chanOut,'UniformOutput',false);
        
        % Add AWGN
        rxSig = cellfun(@awgn,chanOut, ...
            num2cell(snr*ones(numUsers,1)),'UniformOutput',false);
        
        %% Data Recovery Per User
        % The receive signals for each user are processed individually. The example
        % assumes that there are no front-end impairments and that the transmit
        % configuration is known by the receiver for simplicity.
        %
        % A user number specifies the user of interest being decoded for the
        % transmission. This is also used to index into the vector properties of
        % the configuration object that are user-specific.
        
        % Get field indices from configuration, assumed known at receiver
        ind = wlanFieldIndices(cfgVHTMU);
        
        % Single-user receivers recover payload bits
        rxDataBits = cell(numUsers,1);
        GMatrices = cell(numUsers,1);
        for uIdx = 1:numUsers
            rxNSig = rxSig{uIdx}(chanDelay(uIdx)+1:end, :);
            
            % User space-time streams
            stsU = numSTSVec(uIdx);
            
            % Estimate noise power in VHT fields
            lltf = rxNSig(ind.LLTF(1):ind.LLTF(2),:);
            demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
            nVar = helperNoiseEstimate(demodLLTF,chanBW,sum(numSTSVec));
            % fprintf('noise %.5f\n', nVar);
            
            % Perform channel estimation based on VHT-LTF
            rxVHTLTF  = rxNSig(ind.VHTLTF(1):ind.VHTLTF(2),:);
            demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF,chanBW,numSTSVec);
            chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF,chanBW,numSTSVec); % H tilde Nst-by-Nsts-by-Nr
        
            % Recover information bits in VHT Data field
            rxVHTData = rxNSig(ind.VHTData(1):ind.VHTData(2),:);
            [rxDataBits{uIdx},~,eqsym,GMatrices{uIdx}] = wlanVHTDataRecoverFM(rxVHTData, ...
                chanEst,nVar,cfgVHTMU,uIdx, ...
                'EqualizationMethod',eqMethod,'PilotPhaseTracking','None', ...
                'LDPCDecodingMethod','layered-bp','MaximumLDPCIterationCount',6);
        end
        
        %%
        % Per-stream equalized symbol constellation plots validate the simulation
        % parameters and convey the effectiveness of the technique. Note the
        % discernible 16QAM, 64QAM and QPSK constellations per user as specified on
        % the transmit end. Also observe the EVM degradation over the different
        % streams for an individual user. This is a representative characteristic
        % of the channel inversion technique.
        %
        % The recovered data bits are compared with the transmitted payload bits to
        % determine the bit error rate.
        
        % Compare recovered bits against per-user APEPLength information bits
        ber = inf(1, numUsers);
        for uIdx = 1:numUsers
            idx = (1:cfgVHTMU.APEPLength(uIdx)*8).';
            [~,ber(uIdx)] = biterr(txDataBits{uIdx}(idx),rxDataBits{uIdx}(idx));
            % disp(['Bit Error Rate for User ' num2str(uIdx) ': ' num2str(ber(uIdx))]);
        end
        
        % rng(s); % Restore RNG state
        
        %% Save data
        dict.snr = snr;
        dict.mat_rec_no_attack = mat_rec_no_attack;
        dict.mat_rec = mat_rec;
        dict.sumPowerOrig = sumPowerOrig;
        dict.minObjOrig = minObjOrig;
        dict.sumPowerAttack = sumPowerAttack;
        dict.minObjAttack = minObjAttack;
        dict.ber = ber;
        dict.num_iterations = num_iterations;
        
        save(name_save_file, '-struct', 'dict');
    end
end

% for i=1:numTx
%     figure(); 
%     for ii=1:numUsers
%         plot(abs(mat_rec_no_attack{ii}(:, i, 1)), 'DisplayName', strcat('user', num2str(ii)))
%         hold on;
%     end
%     legend();
%     title(['V matrix attack - TX ', num2str(i)])
% end
% 
% for i=1:numTx
%     figure(); 
%     for ii=1:numUsers
%         plot(imag(mat_rec{ii}(:, i, 1)), 'DisplayName', strcat('user', num2str(ii)))
%         hold on;
%     end
%     legend();
%     title(['V matrix attack - TX ', num2str(i)])
% end
