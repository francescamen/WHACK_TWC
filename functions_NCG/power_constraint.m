function sumPower = power_constraint(x, mat_rec_victim, numST, numSTSVec)

mat_rec = cell(2, 1);
mat_rec{1} = mat_rec_victim;
mat_rec{2} = x;
v = zeros(numST,sum(numSTSVec),sum(numSTSVec));
%   Nst-by-Nt-by-Nsts
for uIdx = 1:2
    stsIdx = sum(numSTSVec(1:uIdx-1))+(1:numSTSVec(uIdx));
    v(:,:,stsIdx) = mat_rec{uIdx};     % Nst-by-Nt-by-Nsts
end
% steeringMatrixAttack = cat(3, mat_rec_victim_complex, mat_rec_attacker);

steeringMatrixAttack = zeros(numST,sum(numSTSVec),sum(numSTSVec));
for i = 1:numST
    % Channel inversion precoding
    v_i = squeeze(v(i,:,:));
    if det(v_i'*v_i) == 0
        sumPower = -100;
        return
    end   
    steeringMatrixAttack(i,:,:) = v_i/(v_i'*v_i);
end
powerST = zeros(numST, 1);
for st=1:numST
    powerST(st) = trace(squeeze(steeringMatrixAttack(st, :, :))*squeeze(steeringMatrixAttack(st, :, :))');
end
sumPower = sum(powerST);

end