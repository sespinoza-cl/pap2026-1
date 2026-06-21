function MI_vec = tort_mi_batch(phase, amp_matrix, n_bins)
% tort_mi_batch  Tort MI vectorizado para múltiples señales de amplitud.
%
%   phase      : [1 × T] ángulos de fase EMG [-pi, pi]
%   amp_matrix : [N_A × T] envolventes de amplitud EEG (una fila por banda/freq)
%   n_bins     : número de bins de fase (típico: 18)
%   MI_vec     : [N_A × 1] índice de modulación Tort para cada señal de amplitud
%
% Implementación idéntica a tort_mi.m pero opera sobre N_A señales en paralelo.
% La fórmula es: MI = (log(N) - H) / log(N)
% donde H = -sum(p_j * log(p_j)) y p_j = amplitud media normalizada en bin j.

    edges  = linspace(-pi, pi, n_bins + 1);
    N_A    = size(amp_matrix, 1);
    amp_m  = zeros(N_A, n_bins);   % amplitud media por bin [N_A × n_bins]

    for k = 1:n_bins
        ik = phase >= edges(k) & phase < edges(k+1);
        if any(ik)
            amp_m(:, k) = mean(amp_matrix(:, ik), 2);   % [N_A × 1]
        end
    end

    row_sum = sum(amp_m, 2) + eps;        % [N_A × 1]
    p       = amp_m ./ row_sum;           % [N_A × n_bins] normalizado
    H       = -sum(p .* log(p + eps), 2); % [N_A × 1] entropía
    MI_vec  = (log(n_bins) - H) / log(n_bins);  % [N_A × 1]
end
