%% SGap1_chi_emg.m — Gap #1 (alta prioridad): Δχ × EMG masetero RMS
% Controla si el aplanamiento del exponente aperiódico (χ 1.110→0.735)
% es artefacto muscular: si Δχ correlaciona con la amplitud del EMG masetero,
% el efecto puede ser espurio. Si n.s., el exponent flattening es cortical.
%
% Δχ = χ_Nc - χ_Ch (positivo = aplanamiento en Chew)
% EMG RMS ya está computado en S6b_A3_muscle.mat (rectified masseter, 20-200 Hz)
%
% Análisis adicional: el aplanamiento Δχ versus Δtheta(ROI) — consistencia
% mecanismo (ambos deberían ir juntos si son corticales).
%
% Input:  data/computed/v1_S3_specparam_results.csv
%         outputs/stats/S6b_A3_muscle.mat
% Output: outputs/reviewer/SGap1_chi_emg_result.txt
%
% Correr desde Analysis_V1_Final/:
%   cd('.../Analysis_V1_Final'); run S0_config; run code/SGap1_chi_emg

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

fprintf('\n=== SGap1_chi_emg.m — Gap #1: Δχ × EMG RMS ===\n');

%% ── Cargar Δχ desde CSV ──────────────────────────────────────────────────
fprintf('\nCargando FOOOF CSV: %s\n', FILE_CSV);
T = readtable(FILE_CSV);

% Filtrar solo Cases
T_cas = T(strcmp(T.group, 'Cases'), :);
subj_list = unique(T_cas.subject, 'stable');
n_subj = numel(subj_list);
fprintf('Sujetos Cases en CSV: %d (esperado %d)\n', n_subj, N_CASES);

% Extraer exponente por condición y calcular Δχ = χ_Nc - χ_Ch
chi_nc = nan(n_subj, 1);
chi_ch = nan(n_subj, 1);
for s = 1:n_subj
    row_nc = T_cas(strcmp(T_cas.subject, subj_list{s}) & ...
                   strcmp(T_cas.condition, 'Nc'), :);
    row_ch = T_cas(strcmp(T_cas.subject, subj_list{s}) & ...
                   strcmp(T_cas.condition, 'Ch'), :);
    if ~isempty(row_nc), chi_nc(s) = row_nc.ap_exponent(1); end
    if ~isempty(row_ch), chi_ch(s) = row_ch.ap_exponent(1); end
end
delta_chi = chi_nc - chi_ch;   % positivo = aplanamiento (chi cae con chewing)

n_valid = sum(~isnan(delta_chi));
fprintf('Δχ (χ_Nc - χ_Ch): media=%.3f ± %.3f | mediana=%.3f (n=%d)\n', ...
    mean(delta_chi,'omitnan'), std(delta_chi,'omitnan'), ...
    median(delta_chi,'omitnan'), n_valid);
fprintf('Esperado: ~0.375 (= 1.110 - 0.735 del canónico)\n');

if abs(mean(delta_chi,'omitnan') - 0.375) > 0.05
    warning('Δχ medio (%.3f) difiere del canónico (0.375). Verificar CSV.', ...
        mean(delta_chi,'omitnan'));
end

%% ── Cargar EMG RMS desde S6b ─────────────────────────────────────────────
S6b_file = fullfile(OUT_STATS, 'S6b_A3_muscle.mat');
if ~exist(S6b_file, 'file')
    error(['S6b_A3_muscle.mat no encontrado.\n' ...
           'Correr primero: run code/S6b_A3_muscle\n']);
end
S6b = load(S6b_file, 'emg_rms', 'dtheta');
emg_rms = S6b.emg_rms(:);    % [31×1] RMS masetero, 20-200 Hz, bloque chewing
dtheta  = S6b.dtheta(:);     % [31×1] Δtheta ROI (Ch-Nc, dB)

fprintf('\nEMG RMS cargado de S6b: n=%d, media=%.4g ± %.4g\n', ...
    sum(~isnan(emg_rms)), mean(emg_rms,'omitnan'), std(emg_rms,'omitnan'));

% Alinear sujetos: S6b usa orden CASES cell array, CSV puede diferir
% Verificar que n_subj == N_CASES y usar mismo orden
if n_subj ~= N_CASES
    warning('n_subj CSV (%d) != N_CASES (%d). Verificar exclusiones.', n_subj, N_CASES);
end
if numel(emg_rms) ~= N_CASES
    warning('numel(emg_rms) (%d) != N_CASES (%d).', numel(emg_rms), N_CASES);
end

% Usar la longitud mínima para la correlación
n_corr = min([n_valid, numel(emg_rms), numel(dtheta)]);
dc  = delta_chi(1:n_corr);
er  = emg_rms(1:n_corr);
dth = dtheta(1:n_corr);

%% ── Correlaciones principales ────────────────────────────────────────────
fprintf('\n===== GAP #1: Δχ × EMG RMS (Spearman) =====\n');
[rho_emg, p_emg] = corr(dc, er, 'Type', 'Spearman', 'Rows', 'complete');
fprintf('Δχ × EMG masetero RMS: rho=%+.3f  p=%.4f  → %s\n', rho_emg, p_emg, ...
    verdict(p_emg, 'Δχ covaria con EMG — posible artefacto', ...
                    'n.s. (el aplanamiento NO escala con amplitud muscular — OK)'));

fprintf('\n===== Consistencia mecanismo: Δχ × Δtheta(ROI) =====\n');
[rho_th, p_th] = corr(dc, dth, 'Type', 'Spearman', 'Rows', 'complete');
fprintf('Δχ × Δtheta ROI:       rho=%+.3f  p=%.4f  → %s\n', rho_th, p_th, ...
    verdict(p_th, 'Δχ correlaciona con Δtheta — consistente con mecanismo común', ...
                   'n.s. (disociación Δχ y Δtheta — interpretar con cautela)'));

fprintf('\n===== Referencia: Δtheta × EMG RMS (canónico A3) =====\n');
[rho_dth_emg, p_dth_emg] = corr(dth, er, 'Type', 'Spearman', 'Rows', 'complete');
fprintf('Δtheta × EMG RMS:      rho=%+.3f  p=%.4f  (esperado ρ≈-0.22 n.s.)\n', ...
    rho_dth_emg, p_dth_emg);

%% ── Guardar resultados ───────────────────────────────────────────────────
out_file = fullfile(OUT_REVIEWER, 'SGap1_chi_emg_result.txt');
fid = fopen(out_file, 'w');
fprintf(fid, '=== GAP #1: Δχ × EMG RMS — P2V1 ===\n');
fprintf(fid, 'Fecha: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf(fid, 'N_CASES=%d | seed=%d\n\n', N_CASES, RNG_SEED);
fprintf(fid, 'Δχ = χ_Nc - χ_Ch (>0 = aplanamiento durante chewing)\n');
fprintf(fid, 'Δχ media=%.3f ± %.3f (mediana=%.3f) | n=%d\n\n', ...
    mean(dc,'omitnan'), std(dc,'omitnan'), median(dc,'omitnan'), n_corr);
fprintf(fid, '--- RESULTADO PRINCIPAL ---\n');
fprintf(fid, 'Δχ × EMG masetero RMS (Spearman): rho=%+.3f  p=%.4f\n', rho_emg, p_emg);
fprintf(fid, 'Veredicto: %s\n\n', verdict(p_emg, ...
    'SIG: Δχ covaria con EMG — revisar interpretación cortical del exponente', ...
    'n.s. (OK): aplanamiento aperiodico no escala con amplitud muscular'));
fprintf(fid, '--- CONSISTENCIA MECANISMO ---\n');
fprintf(fid, 'Δχ × Δtheta ROI (Spearman): rho=%+.3f  p=%.4f\n', rho_th, p_th);
fprintf(fid, 'Δtheta × EMG RMS (ref A3):  rho=%+.3f  p=%.4f (canonical ≈-0.22 n.s.)\n', ...
    rho_dth_emg, p_dth_emg);
fprintf(fid, '\n--- INTERPRETACIÓN ---\n');
if p_emg >= 0.05
    fprintf(fid, 'El aplanamiento de χ no escala con la amplitud del masetero.\n');
    fprintf(fid, 'Junto con A3 (Δtheta n.s. con EMG), ambos efectos se disocian\n');
    fprintf(fid, 'del artefacto muscular → interpretación cortical sostenida.\n');
else
    fprintf(fid, 'ATENCIÓN: Δχ correlaciona con EMG RMS (rho=%+.3f p=%.4f).\n', rho_emg, p_emg);
    fprintf(fid, 'Considerar añadir esta limitación al texto (Section Limitations).\n');
end
fclose(fid);
fprintf('\nGuardado: %s\n', out_file);

prov = struct('script', mfilename('fullpath'), 'date', datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
    'fooof_csv', FILE_CSV, 's6b_mat', S6b_file);
save(fullfile(OUT_STATS,'SGap1_chi_emg.mat'), 'delta_chi','emg_rms','dtheta', ...
    'rho_emg','p_emg','rho_th','p_th','rho_dth_emg','p_dth_emg','prov','-v7.3');

function s = verdict(p, msg_sig, msg_ns)
    if p < 0.05, s = msg_sig; else, s = msg_ns; end
end
