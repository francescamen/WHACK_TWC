function [out_obj_sum, out_obj] = min_det_prob(x, mat_rec_victim, numST)

for st=1:numST      
    det_expr = squeeze(mat_rec_victim(st, :, :))*squeeze(mat_rec_victim(st, :, :))'*squeeze(x(st, :, :))*squeeze(x(st, :, :))' ...
        - squeeze(mat_rec_victim(st, :, :))*squeeze(x(st, :, :))'*squeeze(x(st, :, :))*squeeze(mat_rec_victim(st, :, :))';
    
    out_obj(st) = det_expr * det_expr'; % det * det'
end
out_obj_sum = sum(out_obj);

end
