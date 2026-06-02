%% get_fooof_offset.m — extrae offset b desde AP fit + exponent
% b = AP(fmin) + χ·log10(fmin)   [modelo FOOOF fixed: L(f) = b - χ·log10(f)]
clear; clc;
P  = S0_paths();
W  = load(P.file_fooof_ws);
GR = W.GR;

f = GR.Cases.f(:);
[f_min, i_min] = min(f);   % primer bin de frecuencia del fit (e.g. 3 Hz)

% AP_Ch / AP_Nc = curva aperiódica ajustada [nFreq × nSubj]
AP_ch = real(GR.Cases.AP_Ch);   % [nFreq × n]
AP_nc = real(GR.Cases.AP_Nc);

exp_ch = double(GR.Cases.exp_Ch(:));   % exponent χ, Chew
exp_nc = double(GR.Cases.exp_Nc(:));   % exponent χ, NoChew

% offset b = AP(fmin) + χ · log10(fmin)
b_ch = AP_ch(i_min,:)' + exp_ch .* log10(f_min);
b_nc = AP_nc(i_min,:)' + exp_nc .* log10(f_min);

n = sum(~isnan(b_ch));
[p_b, ~] = signrank(b_ch, b_nc);
d_b = (mean(b_ch,'omitnan') - mean(b_nc,'omitnan')) / ...
      sqrt((std(b_ch,'omitnan')^2 + std(b_nc,'omitnan')^2)/2);

fprintf('\n=== FOOOF Aperiodic OFFSET (b)  —  fmin=%.1f Hz ===\n', f_min);
fprintf('  Cases No-Chew : %.3f ± %.3f\n', mean(b_nc,'omitnan'), std(b_nc,'omitnan'));
fprintf('  Cases Chew    : %.3f ± %.3f\n', mean(b_ch,'omitnan'), std(b_ch,'omitnan'));
fprintf('  Wilcoxon Ch vs Nc : p = %.4f\n', p_b);
fprintf('  Cohen d           = %.3f\n\n', d_b);

% También R² del ajuste (fit quality, responde al comentario M7)
r2_ch = 1 - sum((real(GR.Cases.PSD_Ch) - AP_ch - real(GR.Cases.Res_Ch)).^2, 1) ./ ...
            sum((real(GR.Cases.PSD_Ch)  - mean(real(GR.Cases.PSD_Ch),1)).^2,  1);
r2_nc = 1 - sum((real(GR.Cases.PSD_Nc) - AP_nc - real(GR.Cases.Res_Nc)).^2, 1) ./ ...
            sum((real(GR.Cases.PSD_Nc)  - mean(real(GR.Cases.PSD_Nc),1)).^2,  1);

fprintf('=== FOOOF Fit Quality (R²) ===\n');
fprintf('  Cases Chew    : R² = %.3f ± %.3f\n', mean(r2_ch,'omitnan'), std(r2_ch,'omitnan'));
fprintf('  Cases No-Chew : R² = %.3f ± %.3f\n', mean(r2_nc,'omitnan'), std(r2_nc,'omitnan'));
fprintf('\nPásame todos estos números para insertarlos en Results §3.\n');
