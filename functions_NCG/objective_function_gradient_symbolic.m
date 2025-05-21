function [f_obj_matrix, grad] = objective_function_gradient_symbolic(numSTSVec, ...
    numTx, varargin)

num_users = length(varargin);
mat_rec_i = varargin.';

syms v_i [numTx, sum(numSTSVec)] matrix
for uIdx = 1:num_users
    stsIdx = sum(numSTSVec(1:uIdx-1))+(1:numSTSVec(uIdx));
    v_i(:,stsIdx) = mat_rec_i{uIdx};     % Nst-by-Nt-by-Nsts
end

v_i_square = v_i'*v_i; 

% syms f_constr(v1_i, v2_i) [1 1] matrix keepargs
f_obj = det(v_i_square)*det(v_i_square)'; 
f_obj_matrix = symmatrix2sym(f_obj);

% fsurf(f_constr_matrix,[-200 200],'ShowContours','on');

% grad_v2_i1 = diff(f_obj_matrix, v2_i(1));  
% grad_v2_i2 = diff(f_obj_matrix, v2_i(2));

% syms grad(v1_i, v2_i) [sum(numSTSVec(2))*sum(numSTSVec(2))] matrix keepargs
grad = gradient(f_obj_matrix);

% [X, Y] = meshgrid(-1:.1:1,-1:.1:1);
% X = X(1, :).';
% Y = Y(1, :).';
% G1 = subs(grad(1),[symmatrix2sym(v2_i(1)) symmatrix2sym(v2_i(2))],{X,Y});
% G2 = subs(grad(2),[symmatrix2sym(v2_i(1)) symmatrix2sym(v2_i(2))],{X,Y});
% quiver(X,Y,G1,G2)

end