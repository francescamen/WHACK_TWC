function [xNew, num_iterations] = applyNCG_smart(x, size_x, unperturbed_mat, numST, numSTSVec, numTx, ...
    numSTSelected, mat_rec_victim, power_limit, functpowerST, gradPower, ...
    functobjectiveST, gradObjective, varargin)

[functObstST, gradObstST] = assign_value_to_gradient(unperturbed_mat, mat_rec_victim, ...
    numST, numSTSVec, numTx, functobjectiveST, gradObjective, varargin);
% figure(); plot(abs(gradObstST(1, :)));

xNew = x;
[perturbedSTs, xNew] = selectST_smart(xNew, unperturbed_mat, numST, numSTSelected, functObstST, gradObstST);
% figure(); plot(abs(x(:, 1))); hold on; plot(abs(xNew(:, 1))); figure(); plot(abs(mat_rec_victim(:, 1))); hold on; plot(abs(unperturbed_mat(:, 1))); 

[powerST, gradpowerST] = assign_value_to_gradient(xNew, mat_rec_victim, ...
    numST, numSTSVec, numTx, functpowerST, gradPower, varargin);
% figure(); plot(abs(powerST(1, :)));
% figure(); plot(abs(gradpowerST(1, :)));

deltas = - gradpowerST.';

num_iterations = 0;
sumPower = sum(powerST);
mu_update = 0.5;

while sumPower > power_limit || sumPower < 0

    num_iterations = num_iterations + 1;
    % disp(num_iterations)

    delats_norm = deltas./abs(sum(deltas, 1));
    
    if num_iterations == 2 % avoids to remain stuck
         error('MyError:MaxNumIterationsReached', 'Maximum number of iterations reached')
    end

    for stream = 1:numTx
        perturbedST = perturbedSTs{stream}; 
        
        xNew_real = real(xNew(perturbedST, stream)) + real(mu_update/mod(num_iterations, 5)*(delats_norm(perturbedST, stream)));
        xNew_imag = imag(xNew(perturbedST, stream)) + imag(mu_update/mod(num_iterations, 5)*(delats_norm(perturbedST, stream)));
        xNew(perturbedST, stream) = xNew_real + 1i*xNew_imag;
    end

    [perturbedSTs, xNew] = selectST_smart(xNew, unperturbed_mat, numST, numSTSelected, functObstST, gradObstST);
    
    % Last antenna real and positive as for the IEEE 802.11 standard
    xNew(:, size_x(2)) = abs(real(xNew(:, size_x(2))));
    
    %constraints
    bound = 1.1;
    xNewBounded_real = max(real(xNew), -bound);
    xNewBounded_real = min(xNewBounded_real, bound);
    xNewBounded_imag = max(imag(xNew), -bound);
    xNewBounded_imag = min(xNewBounded_imag, bound);
    xNew = xNewBounded_real + 1i*xNewBounded_imag;

    [powerST, gradpowerST] = assign_value_to_gradient(xNew, ...
        mat_rec_victim, numST, numSTSVec, numTx, functpowerST, gradPower, varargin);

    sumPower = sum(powerST);
    deltas_new = - gradpowerST.';

    chi = diag((deltas_new - deltas)*deltas_new' ./ (deltas*deltas'));
    chi = max(chi, 0);
    deltas = deltas_new + chi.*deltas;

end

end