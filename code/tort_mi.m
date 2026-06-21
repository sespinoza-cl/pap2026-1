function mi = tort_mi(phase, amp, n_bins)
% tort_mi  Modulation Index de Tort et al. 2010
%   phase: vector columna de ángulos de fase [-pi, pi]
%   amp:   vector columna de amplitud (misma longitud)
%   n_bins: número de bins (típico: 18)
%   mi:    escalar MI normalizado [0, 1/log(n_bins)]
    edges = linspace(-pi, pi, n_bins+1);
    amp_m = zeros(n_bins,1);
    for k = 1:n_bins
        ik = phase >= edges(k) & phase < edges(k+1);
        if any(ik), amp_m(k) = mean(amp(ik)); end
    end
    amp_m = amp_m / (sum(amp_m) + eps);
    H     = -sum(amp_m .* log(amp_m + eps));
    H_max = log(n_bins);
    mi    = (H_max - H) / H_max;
end
