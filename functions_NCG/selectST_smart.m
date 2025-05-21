function [perturbedSTs, childSelectST] = selectST_smart(child, unperturbed_mat, ...
    numST, numSTSelected, functObstST, gradObjST)

% unperturbed_idxs = find(powerST(1, :) > power_limit);
% unperturbed_idxs = sort(unperturbed_idxs, 2, 'descend');
% function_sel = gradObjST;
% function_sel(:, unperturbed_idxs) = [];
% % figure(); plot(abs(powerST(1, :)));
% 
% unperturbed_idxs = find(powerST_unperturbed(1, :) > 10*mean(powerST_unperturbed));
% unperturbed_idxs = sort(unperturbed_idxs, 2, 'descend');
% function_sel(:, unperturbed_idxs) = [];
% % figure(); plot(abs(powerST_unperturbed(1, :)));

childSelectST = child;
function_sel = gradObjST; % figure(); plot(abs(function_sel(1, :)));
perturbedSTs = [];

for tx_ant = 1:size(function_sel, 1)
    [~, sorted_idxs] = sort(function_sel(tx_ant, :), 2, 'descend'); % perturb the ones with highest values, keep the others unchanges
    
    zero_idxs = find(functObstST == 0);
    num_zeros = size(zero_idxs, 2);  % to enforce th constraint that the determinant should not be zero

    perturbedST = [zero_idxs, sorted_idxs(1:numSTSelected-num_zeros)];
    perturbedST = sort(perturbedST, 1, 'descend');
    
    unperturbedST = linspace(1, numST, numST);
    unperturbedST(perturbedST) = [];
    
    childSelectST(unperturbedST, tx_ant) = unperturbed_mat(unperturbedST, tx_ant);
    perturbedSTs{tx_ant} = perturbedST;
end

end

% figure(); plot(gradObjST(1, :)); hold on; scatter(sorted_idxs(1, 1:numSTSelected), gradObjST(1, sorted_idxs(1, 1:numSTSelected)));
% figure(); plot(gradObjST(2, :)); hold on; scatter(sorted_idxs(1, 1:numSTSelected), gradObjST(2, sorted_idxs(1, 1:numSTSelected)));
% figure(); plot(abs(child(:, 1))); hold on; plot(abs(childSelectST(:, 1)));
