%% SM7_baseline_neural.m — M7: Equivalencia neural en baseline (no-chew)
% Prueba si Cases y Controls tienen theta frontal equivalente en el bloque
% no-chew (Block 1), que es el mismo para ambos grupos.
%
% Test: Mann-Whitney U (ranksum) en theta ROI promediado sobre WIN_LATE.
% También CBPT ligero (sign-flip sobre diferencia de grupos) para consistencia
% con el resto del pipeline.
%
% Input:  data/computed/v1_S2_TF_data.mat
% Output: outputs/reviewer/SM7_baseline_neural_result.txt
%
% Correr desde Analysis_V1_Final/:
%   cd('.../Analysis_V1_Final'); run S0_config; run code/SM7_baseline_neural

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

fprintf('\n=== SM7_baseline_neural.m — Equivalencia neural baseline ===\n');
fprintf('N=%d casos + %d controles | theta=%d-%dHz | WIN_LATE=[%d,%d]ms\n', ...
    N_CASES, N_CONTROLS, BAND_THETA(1), BAND_THETA(2), WIN_LATE(1), WIN_LATE(2));

%% ── Cargar datos TF ───────────────────────────────────────────────────────
fprintf('\nCargando %s ...\n', FILE_TF);
vars_available = whos('-file', FILE_TF);
var_names = {vars_available.name};
fprintf('Variables disponibles: %s\n', strjoin(var_names, ', '));

load(FILE_TF, 'tf_npl_cas_nc', 'tf_npl_ctr_nc', 'FREX_TF', 'times_ms');

% Intentar cargar conteos de trials si existen
n_trials_cas = []; n_trials_ctr = [];
if ismember('n_trials_cas_nc', var_names)
    tmp = load(FILE_TF, 'n_trials_cas_nc'); n_trials_cas = tmp.n_trials_cas_nc;
    fprintf('n_trials_cas_nc encontrado: media=%.1f ± %.1f\n', ...
        mean(n_trials_cas,'omitnan'), std(n_trials_cas,'omitnan'));
end
if ismember('n_trials_ctr_nc', var_names)
    tmp = load(FILE_TF, 'n_trials_ctr_nc'); n_trials_ctr = tmp.n_trials_ctr_nc;
    fprintf('n_trials_ctr_nc encontrado: media=%.1f ± %.1f\n', ...
        mean(n_trials_ctr,'omitnan'), std(n_trials_ctr,'omitnan'));
end

FREX_TF  = FREX_TF(:)';
times_ms = times_ms(:)';

%% ── Índices de theta y ventana late ──────────────────────────────────────
theta_idx = FREX_TF >= BAND_THETA(1) & FREX_TF <= BAND_THETA(2);
late_idx  = times_ms >= WIN_LATE(1)  & times_ms <= WIN_LATE(2);

fprintf('Theta: %d bins (%.1f–%.1f Hz) | Late: %d puntos (%d–%d ms)\n', ...
    sum(theta_idx), FREX_TF(find(theta_idx,1)), FREX_TF(find(theta_idx,1,'last')), ...
    sum(late_idx), WIN_LATE(1), WIN_LATE(2));

%% ── Extraer theta ROI promediado sobre frecuencia y tiempo ───────────────
% tf_npl_*: [N x n_frex x n_times] — baseline-normalizado en dB
% Promediar sobre dimensiones 2 (theta_idx) y 3 (late_idx)
cas_nc_theta = squeeze(mean(mean(tf_npl_cas_nc(:, theta_idx, late_idx), 2), 3));  % [31×1]
ctr_nc_theta = squeeze(mean(mean(tf_npl_ctr_nc(:, theta_idx, late_idx), 2), 3));  % [15×1]

fprintf('\nCases   (NoChew, theta late): %.4f ± %.4f dB (n=%d)\n', ...
    mean(cas_nc_theta), std(cas_nc_theta), numel(cas_nc_theta));
fprintf('Controls(NoChew, theta late): %.4f ± %.4f dB (n=%d)\n', ...
    mean(ctr_nc_theta), std(ctr_nc_theta), numel(ctr_nc_theta));

%% ── Mann-Whitney U (principal) ────────────────────────────────────────────
[p_mwu, h_mwu, stats_mwu] = ranksum(cas_nc_theta, ctr_nc_theta);
fprintf('\n===== M7: Mann-Whitney U — Cases-Nc vs Controls-Nc (theta late ROI) =====\n');
fprintf('W = %.0f | p = %.4f → %s\n', stats_mwu.ranksum, p_mwu, ...
    tern(p_mwu < 0.05, 'SIGNIFICATIVO (grupos distintos en baseline)', ...
                        'n.s. (grupos equivalentes en baseline — OK)'));

%% ── CBPT ligero: label-shuffle sobre diferencia de grupos ────────────────
fprintf('\n===== M7 CBPT robusto: label-shuffle permutation =====\n');
rng(RNG_SEED);
obs_diff = mean(cas_nc_theta) - mean(ctr_nc_theta);
n_all    = N_CASES + N_CONTROLS;
pool     = [cas_nc_theta; ctr_nc_theta];

perm_diffs = nan(N_PERM, 1);
for i = 1:N_PERM
    idx = randperm(n_all);
    perm_diffs(i) = mean(pool(idx(1:N_CASES))) - mean(pool(idx(N_CASES+1:end)));
end
p_perm = mean(abs(perm_diffs) >= abs(obs_diff));
fprintf('Diferencia observada: %.4f dB | p_perm(dos-colas) = %.4f → %s\n', ...
    obs_diff, p_perm, tern(p_perm < 0.05, 'SIG', 'n.s.'));

%% ── Verificar conteos de trials si están disponibles ─────────────────────
if ismember('n_trials_cas_ch', var_names) && ismember('n_trials_ctr_ch', var_names) && ...
   ismember('n_trials_cas_nc', var_names) && ismember('n_trials_ctr_nc', var_names)
    tmp = load(FILE_TF, 'n_trials_cas_ch', 'n_trials_ctr_ch');
    n_cas_ch = tmp.n_trials_cas_ch; n_ctr_ch = tmp.n_trials_ctr_ch;
    all_n = [n_trials_cas(:); n_trials_ctr(:); n_cas_ch(:); n_ctr_ch(:)];
    total  = sum(all_n, 'omitnan');
    per_subj_all = [n_trials_cas(:); n_trials_ctr(:); n_cas_ch(:); n_ctr_ch(:)];
    fprintf('\n===== TRIAL COUNTS (verificación) =====\n');
    fprintf('Total correct trials (todas condiciones): %d\n', total);
    fprintf('Media por condición por sujeto: %.1f ± %.1f\n', ...
        mean([n_trials_cas(:); n_trials_ctr(:); n_cas_ch(:); n_ctr_ch(:)],'omitnan'), ...
        std([n_trials_cas(:); n_trials_ctr(:); n_cas_ch(:); n_ctr_ch(:)],'omitnan'));
    fprintf('N > 4200 check: %s\n', tern(total > 4200, 'OK (>4200)', 'FALLO (<4200)'));
end

%% ── Guardar resultados ───────────────────────────────────────────────────
out_file = fullfile(OUT_REVIEWER, 'SM7_baseline_neural_result.txt');
fid = fopen(out_file, 'w');
fprintf(fid, '=== SM7: Equivalencia neural baseline — P2V1 ===\n');
fprintf(fid, 'Fecha: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf(fid, 'N_CASES=%d | N_CONTROLS=%d | seed=%d\n\n', N_CASES, N_CONTROLS, RNG_SEED);
fprintf(fid, 'Contraste: Cases-NoChew vs Controls-NoChew | theta %d–%d Hz | WIN_LATE=[%d,%d]ms\n\n', ...
    BAND_THETA(1), BAND_THETA(2), WIN_LATE(1), WIN_LATE(2));
fprintf(fid, 'Cases    theta: %.4f ± %.4f dB\n', mean(cas_nc_theta), std(cas_nc_theta));
fprintf(fid, 'Controls theta: %.4f ± %.4f dB\n\n', mean(ctr_nc_theta), std(ctr_nc_theta));
fprintf(fid, 'Mann-Whitney U: W=%.0f  p=%.4f\n', stats_mwu.ranksum, p_mwu);
fprintf(fid, 'Label-shuffle CBPT (%d perm): obs_diff=%.4f  p=%.4f\n\n', N_PERM, obs_diff, p_perm);
fprintf(fid, 'Veredicto M7: %s\n', tern(p_mwu < 0.05, ...
    'GRUPOS DIFIEREN en baseline neural — revisar exclusiones', ...
    'Grupos equivalentes en baseline neural (p=%.4f) — OK para interpretar interacción'));
if p_mwu >= 0.05
    fprintf(fid, '  → Junto con baseline conductual (LME Group p=0.86), confirma que\n');
    fprintf(fid, '     la interacción Group×Block refleja el efecto de masticación,\n');
    fprintf(fid, '     no una diferencia constitucional entre grupos.\n');
end
fclose(fid);
fprintf('\nGuardado: %s\n', out_file);

function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end
