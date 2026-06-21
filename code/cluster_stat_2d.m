function [max_stat, mask_out] = cluster_stat_2d(sig_mask, t_map)
% Suma-t del cluster positivo más grande en una máscara 2D
CC = bwconncomp(sig_mask);
max_stat = 0;  mask_out = false(size(sig_mask));  best_cc = [];
for k = 1:CC.NumObjects
    cs = sum(t_map(CC.PixelIdxList{k}));
    if cs > max_stat, max_stat = cs; best_cc = CC.PixelIdxList{k}; end
end
if ~isempty(best_cc), mask_out(best_cc) = true; end
end
