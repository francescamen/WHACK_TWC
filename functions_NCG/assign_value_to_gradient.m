function [functionST, gradST] = assign_value_to_gradient(x, ...
    mat_rec_victim, numST, numSTSVec, numTx, funct_sym, grad, varargin)

varargin = varargin{1};
functionST = zeros(1, numST, 1);
gradST = zeros(numTx, numST);  % gradient with respect to the feedback of the adversary (control variables)
num_users = length(varargin);

v_i_cell = cell(num_users*numTx, 1);
for ind_user=1:num_users
    for ind_tx=1:numTx
        v_i_cell{ind_tx + numTx*(ind_user-1)} = varargin{ind_user}(ind_tx);
    end
end

if ~ iscell(mat_rec_victim)
    mat_rec_victim = {mat_rec_victim};
end

victims_attacker_matrix = cell(numST, num_users*numTx, 1);
for ind_user=1:num_users - 1  % all the victims, remove the attacker
    for ind_user_ST = 1:numST
        for ind_tx=1:numTx
            victims_attacker_matrix{ind_user_ST, ind_tx + numTx*(ind_user-1)} = mat_rec_victim{ind_user}(ind_user_ST, ind_tx);
        end
    end
end
for ind_user_ST = 1:numST
    for ind_tx=1:numTx
        victims_attacker_matrix{ind_user_ST, ind_tx + numTx*(num_users-1)} = x(ind_user_ST, ind_tx);
    end
end

%% loop - assign the values to all subcarriers
for i = 1:numST  % can be parfor
    functionST(i) = double(subs(funct_sym, ...
        v_i_cell, ...
        victims_attacker_matrix(i, :).'));

    if ~isempty(grad)
        gradST(:, i) = double(subs(grad(numTx*(num_users-1)+1:end), v_i_cell, ...
            victims_attacker_matrix(i, :).')); 
    end
end

end