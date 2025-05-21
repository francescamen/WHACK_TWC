function [Vtilde_rec, Vtildecompl] = beamformingFeedbackQuantization_complete(V, Nr, Nc)

%% parameters 

psi_bit = 7;
phi_bit = psi_bit + 2;

const1_phi = 2^(phi_bit-1);
const2_phi = 2^(phi_bit);

const1_psi = 2^(psi_bit+1);
const2_psi = 2^(psi_bit+2);

%% computation

D_i_matrices = {};
G_li_matrices = {};

vm_angles_vector = angle(V(Nr, :)); 
Dtilde = diag(exp(1i*vm_angles_vector));
Omega = V * Dtilde';
for i = 1:min(Nc, Nr-1)
    %% compute phi and build D matrix
    D_i = eye(Nr);
    for l = i:Nr-1
        phi_li = angle(Omega(l, i));

        quantized_phi_li = quantizeVmat(phi_li, const1_phi, const2_phi);
        phi_li_rad = inverseQuantizeVmat(quantized_phi_li, const1_phi, const2_phi);

        D_i(l, l) = exp(1i*phi_li_rad);
    end
    D_i_matrices{i} = D_i;
    Omega = D_i' * Omega;

    %% compute psi and build G matrices
    for l = i+1:Nr
        G_li = eye(Nr);
        psi_li = acos(Omega(i, i)/ sqrt(abs(Omega(i, i))^2 + abs(Omega(l, i))^2));

        quantized_psi_li = quantizeVmat(psi_li, const1_psi, const2_psi);
        psi_li_rad = inverseQuantizeVmat(quantized_psi_li, const1_psi, const2_psi);

        G_li(i, i) = cos(psi_li_rad);
        G_li(l, l) = cos(psi_li_rad);
        G_li(i, l) = sin(psi_li_rad);
        G_li(l, i) = -sin(psi_li_rad);
        G_li_matrices{l, i} = G_li;
        Omega = G_li * Omega;
    end
end   

%% reconstruct V tilde matrix
I_matrix = eye(Nr, Nc);
Vtildecompl = eye(Nr, Nr);
for i = 1:min(Nc, Nr-1)
    Vtildecompl = Vtildecompl * cell2mat(D_i_matrices(i));
    for l = i+1:Nr
        Vtildecompl = Vtildecompl * cell2mat(G_li_matrices(l, i)).';
    end
end
Vtilde_rec = Vtildecompl * I_matrix;

% Vreconstruct = Vtilde * Dtilde;

end