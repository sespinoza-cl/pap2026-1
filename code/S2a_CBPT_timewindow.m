%% S2a_CBPT_timewindow.m  —  FASE 1 Paso 1: Ventana temporal significativa
% CBPT 2D (frecuencia x tiempo) dentro de cada banda, restringido a
% WIN_ANALYSIS=[-200,1500]ms. Metodología consistente con S2_TF_cluster.m:
%   - ALPHA_THRESH=0.05 para voxel (igual que Analysis_paper/S0_config)
%   - cluster_stat_2d.m para detectar clusters 2D (freq×time)
%   - N_PERM=5000, sign-flip permutation
%
% Input:  data/computed/v1_S2_TF_data.mat (N=31, freq=200, ROI FC promediado)
% Output: outputs/stats/S2a_CBPT_timewindow.mat
%         outputs/reviewer/S2a_CBPT_timewindow_result.txt
%
% La ventana del trial real es estimulo 200ms + ISI 2000ms.
% El TF fue computado sobre [-2000,2496]ms para margen de borde wavelet.
% Se restringe a WIN_ANALYSIS=[-200,1500]ms para el analisis estadistico.

% Correr siempre desde Analysis_V1_Final/ (no desde code/)
%   cd('.../Analysis_V1_Final'); run S0_config; run code/S2a_CBPT_timewindow
if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file')
        run('S0_config.m');
    else
        error('S0_config.m no encontrado. Correr desde Analysis_V1_Final/');
    end
end

% Umbral de voxel para CBPT 2D: 0.05 (igual que S2_TF_cluster.m original)
% NOTA: ALPHA_VOXEL=0.01 se reserva para el fullspace 3D (S2_TF_fullspace)
ALPHA_THRESH = 0.05;

% cluster_stat_2d.m esta en code/ (copiado ahi para disponibilidad en parfor)
% pctRunOnAll garantiza que los workers lo tengan incluso si el pool ya estaba abierto

fprintf('\n=== S2a_CBPT_timewindow.m — FASE 1 Paso 1 (2D CBPT) ===\n');
fprintf('N=%d casos + %d controles | N_PERM=%d | alpha_voxel=%.2f | seed=%d\n', ...
    N_CASES, N_CONTROLS, N_PERM, ALPHA_THRESH, RNG_SEED);
fprintf('WIN_ANALYSIS = [%.0f, %.0f] ms\n', WIN_ANALYSIS(1), WIN_ANALYSIS(2));

%% ── Cargar datos ──────────────────────────────────────────────────────────
fprintf('\nCargando %s ...\n', FILE_TF);
load(FILE_TF, 'tf_npl_cas_ch','tf_npl_cas_nc', ...
              'tf_npl_ctr_ch','tf_npl_ctr_nc', ...
              'FREX_TF','times_ms');

FREX_TF  = FREX_TF(:)';
times_ms = times_ms(:)';

assert(size(tf_npl_cas_ch,1) == N_CASES, ...
    'N_CASES mismatch: datos=%d, config=%d', size(tf_npl_cas_ch,1), N_CASES);
assert(size(tf_npl_ctr_ch,1) == N_CONTROLS, ...
    'N_CONTROLS mismatch: datos=%d, config=%d', size(tf_npl_ctr_ch,1), N_CONTROLS);

fprintf('Datos: [%d subj x %d frex x %d times] | epoca=[%.0f,%.0f]ms @ %.0fHz\n', ...
    N_CASES, size(tf_npl_cas_ch,2), size(tf_npl_cas_ch,3), ...
    times_ms(1), times_ms(end), 1000/(times_ms(2)-times_ms(1)));

%% ── Trim a ventana de análisis del trial ─────────────────────────────────
tidx = times_ms >= WIN_ANALYSIS(1) & times_ms <= WIN_ANALYSIS(2);
times_anal = times_ms(tidx);
n_times = numel(times_anal);
fprintf('Ventana analisis: [%.0f,%.0f]ms → %d puntos\n', ...
    WIN_ANALYSIS(1), WIN_ANALYSIS(2), n_times);

% Trim matrices de tiempo (dim 3)
cas_ch = tf_npl_cas_ch(:,:,tidx);   % [N x 200 x n_times]
cas_nc = tf_npl_cas_nc(:,:,tidx);
ctr_ch = tf_npl_ctr_ch(:,:,tidx);
ctr_nc = tf_npl_ctr_nc(:,:,tidx);

%% ── Indices de bandas ─────────────────────────────────────────────────────
theta_idx = FREX_TF >= BAND_THETA(1) & FREX_TF <= BAND_THETA(2);
alpha_idx = FREX_TF >= BAND_ALPHA(1) & FREX_TF <= BAND_ALPHA(2);
beta_idx  = FREX_TF >= BAND_BETA(1)  & FREX_TF <= BAND_BETA(2);

fprintf('Theta: %d frec bins | Alpha: %d bins | Beta: %d bins\n', ...
    sum(theta_idx), sum(alpha_idx), sum(beta_idx));

%% ── Funcion: 2D CBPT para una banda ──────────────────────────────────────
function [p_clust, stat_obs, clust_mask, t_obs_out, p_obs_out] = ...
        run_cbpt_2d(data_NxFxT, N_subj, n_perm, alpha_thresh, rng_seed, label)
    % data_NxFxT: [N x n_freq x n_times]
    % Retorna: p_cluster, cluster_stat, mascara 2D [n_freq x n_times], t_obs, p_obs
    [~, n_freq, n_times] = size(data_NxFxT);
    rng(rng_seed, 'twister');

    % Estadistico observado
    [~, p_obs] = ttest(data_NxFxT);           % [1 x n_freq x n_times]
    p_obs = squeeze(p_obs);                    % [n_freq x n_times]
    mn    = squeeze(mean(data_NxFxT));         % [n_freq x n_times]
    sd    = squeeze(std(data_NxFxT));
    t_obs = mn ./ (sd / sqrt(N_subj));

    sig_obs = p_obs < alpha_thresh & t_obs > 0;
    [stat_obs, clust_mask] = cluster_stat_2d(sig_obs, t_obs);

    % Distribucion nula
    signs_all = (rand(n_perm, N_subj) > 0.5)*2 - 1;   % [n_perm x N]
    null_dist = zeros(n_perm, 1);

    parfor pp = 1:n_perm
        sp       = signs_all(pp,:)';
        data_p   = bsxfun(@times, data_NxFxT, reshape(sp, [N_subj,1,1]));
        mn_p     = squeeze(mean(data_p));
        sd_p     = squeeze(std(data_p));
        t_p      = mn_p ./ (sd_p / sqrt(N_subj));
        p_p      = 2*(1-tcdf(abs(t_p), N_subj-1));
        [cs, ~]  = cluster_stat_2d(p_p < alpha_thresh & t_p > 0, t_p);
        null_dist(pp) = cs;
    end

    p_clust = mean(null_dist >= stat_obs);
    p_obs_out = p_obs;
    t_obs_out = t_obs;

    fprintf('  %s: stat=%.2f, p=%.4f\n', label, stat_obs, p_clust);
end

%% ── Pool paralelo ─────────────────────────────────────────────────────────
if isempty(gcp('nocreate'))
    parpool('local', min(12, feature('numcores')));
end
pctRunOnAll addpath(fullfile(ROOT_FINAL, 'code'));

%% ── CBPT 2D por banda — CASOS (Chew > NoChew) ────────────────────────────
fprintf('\n--- CASOS (Chew > NoChew) ---\n');

diff_cas = cas_ch - cas_nc;   % [N x 200 x n_times]

[p_cas_theta, stat_cas_theta, mask_cas_theta, t_cas_theta, pv_cas_theta] = ...
    run_cbpt_2d(diff_cas(:,theta_idx,:), N_CASES, N_PERM, ALPHA_THRESH, RNG_SEED,   'Theta  4-7Hz  Casos');
[p_cas_alpha, stat_cas_alpha, mask_cas_alpha, t_cas_alpha, pv_cas_alpha] = ...
    run_cbpt_2d(diff_cas(:,alpha_idx,:), N_CASES, N_PERM, ALPHA_THRESH, RNG_SEED+1, 'Alpha 8-13Hz  Casos');
[p_cas_beta,  stat_cas_beta,  mask_cas_beta,  t_cas_beta,  pv_cas_beta]  = ...
    run_cbpt_2d(diff_cas(:,beta_idx, :), N_CASES, N_PERM, ALPHA_THRESH, RNG_SEED+2, 'Beta 13-30Hz  Casos');

%% ── CBPT 2D — CONTROLES ───────────────────────────────────────────────────
fprintf('\n--- CONTROLES (Chew > NoChew) ---\n');
diff_ctr = ctr_ch - ctr_nc;

[p_ctr_theta, stat_ctr_theta, mask_ctr_theta, t_ctr_theta, ~] = ...
    run_cbpt_2d(diff_ctr(:,theta_idx,:), N_CONTROLS, N_PERM, ALPHA_THRESH, RNG_SEED+3, 'Theta Controles');

%% ── Extraer ventana temporal del cluster ─────────────────────────────────
% Colapsar mascara de cluster sobre frecuencias -> mascara de tiempo
function win = cluster_to_timewin(clust_mask, times_vec, p_val, alpha_clust)
    if p_val < alpha_clust && any(clust_mask(:))
        time_mask = any(clust_mask, 1);     % [1 x n_times]
        t_sig = times_vec(time_mask);
        win = [min(t_sig), max(t_sig)];
    else
        win = [];
    end
end

win_cas_theta = cluster_to_timewin(mask_cas_theta, times_anal, p_cas_theta, ALPHA_CLUST);
win_cas_alpha = cluster_to_timewin(mask_cas_alpha, times_anal, p_cas_alpha, ALPHA_CLUST);
win_cas_beta  = cluster_to_timewin(mask_cas_beta,  times_anal, p_cas_beta,  ALPHA_CLUST);
win_ctr_theta = cluster_to_timewin(mask_ctr_theta, times_anal, p_ctr_theta, ALPHA_CLUST);

%% ── WIN_CBPT canónico ─────────────────────────────────────────────────────
if ~isempty(win_cas_theta)
    WIN_CBPT = win_cas_theta;
    fprintf('\n✓ WIN_CBPT = [%.0f, %.0f] ms (cluster theta casos, p=%.4f)\n', ...
        WIN_CBPT(1), WIN_CBPT(2), p_cas_theta);
elseif ~isempty(win_cas_alpha)
    WIN_CBPT = win_cas_alpha;
    fprintf('\nWARN: theta n.s. WIN_CBPT desde alpha: [%.0f, %.0f] ms\n', WIN_CBPT(1), WIN_CBPT(2));
else
    WIN_CBPT = WIN_LATE;
    fprintf('\nWARN: ninguna banda sig. WIN_CBPT = WIN_LATE fallback.\n');
end

%% ── Guardar ───────────────────────────────────────────────────────────────
out_file = fullfile(OUT_STATS, 'S2a_CBPT_timewindow.mat');
save(out_file, ...
    'p_cas_theta','stat_cas_theta','mask_cas_theta','t_cas_theta','pv_cas_theta','win_cas_theta', ...
    'p_cas_alpha','stat_cas_alpha','mask_cas_alpha','t_cas_alpha','pv_cas_alpha','win_cas_alpha', ...
    'p_cas_beta', 'stat_cas_beta', 'mask_cas_beta', 't_cas_beta', 'pv_cas_beta', 'win_cas_beta', ...
    'p_ctr_theta','stat_ctr_theta','mask_ctr_theta','t_ctr_theta','win_ctr_theta', ...
    'WIN_CBPT','times_anal','FREX_TF','ALPHA_THRESH', ...
    'theta_idx','alpha_idx','beta_idx', ...
    'N_CASES','N_CONTROLS','N_PERM','RNG_SEED','WIN_ANALYSIS');
fprintf('\nGuardado: %s\n', out_file);

%% ── Figura ───────────────────────────────────────────────────────────────
FREX_theta = FREX_TF(theta_idx);
FREX_alpha = FREX_TF(alpha_idx);
FREX_beta  = FREX_TF(beta_idx);

fig = figure('Name','S2a CBPT 2D Ventana Temporal','Units','normalized',...
             'Position',[0.01 0.02 0.97 0.88]);

configs = {
    t_cas_theta, mask_cas_theta, p_cas_theta, win_cas_theta, FREX_theta, sprintf('Theta (4-7 Hz) Casos | p=%.4f',p_cas_theta);
    t_cas_alpha, mask_cas_alpha, p_cas_alpha, win_cas_alpha, FREX_alpha, sprintf('Alpha (8-13 Hz) Casos | p=%.4f',p_cas_alpha);
    t_cas_beta,  mask_cas_beta,  p_cas_beta,  win_cas_beta,  FREX_beta,  sprintf('Beta (13-30 Hz) Casos | p=%.4f',p_cas_beta);
    t_ctr_theta, mask_ctr_theta, p_ctr_theta, win_ctr_theta, FREX_theta, sprintf('Theta Controles | p=%.4f',p_ctr_theta);
};

for b = 1:4
    subplot(2,2,b);
    t_map = configs{b,1};   % [n_freq x n_times]
    mask  = configs{b,2};
    p_val = configs{b,3};
    win   = configs{b,4};
    frex  = configs{b,5};
    ttl   = configs{b,6};

    % TF map
    imagesc(times_anal, frex, t_map);
    axis xy; colormap(jet); colorbar;
    set(gca,'CLim',[-3 5]);
    hold on;

    % Contorno del cluster
    if ~isempty(mask) && any(mask(:)) && p_val < ALPHA_CLUST
        contour(times_anal, frex, double(mask), 1, 'w-', 'LineWidth', 2);
    end

    % Líneas de referencia
    xline(0, 'k--', 'LineWidth', 1.5);
    xline(WIN_EARLY(1), 'w:', 'LineWidth', 1);
    xline(WIN_LATE(1),  'w:', 'LineWidth', 1);
    xline(WIN_LATE(2),  'w:', 'LineWidth', 1);

    % Cuadro de WIN_CBPT si es este el resultado canónico
    if b == 1 && ~isempty(WIN_CBPT)
        rectangle('Position',[WIN_CBPT(1),frex(1),diff(WIN_CBPT),frex(end)-frex(1)], ...
            'EdgeColor','k','LineWidth',2,'LineStyle','--');
    end

    xlabel('Time (ms)'); ylabel('Freq (Hz)');
    title(ttl, 'FontSize', FIG_FS);
    xlim([times_anal(1) times_anal(end)]);
    set(gca, 'FontName', FIG_FONT, 'FontSize', FIG_FS);
end

sgtitle(sprintf('CBPT 2D (freq×time) — N=%d casos | alpha\\_voxel=%.2f | N\\_PERM=%d', ...
    N_CASES, ALPHA_THRESH, N_PERM), 'FontWeight','bold');
saveas(fig, fullfile(OUT_FIGS,'S2a_CBPT_timewindow.png'));
fprintf('Figura guardada.\n');

%% ── Reporte para el revisor ───────────────────────────────────────────────
fid = fopen(fullfile(OUT_REVIEWER,'S2a_CBPT_timewindow_result.txt'),'w');
fprintf(fid,'=== S2a CBPT 2D Ventana Temporal — P2V1 Revision ===\n');
fprintf(fid,'Fecha: %s\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf(fid,'N_CASES=%d | N_CONTROLS=%d | N_PERM=%d | alpha_voxel=%.2f | seed=%d\n\n', ...
    N_CASES, N_CONTROLS, N_PERM, ALPHA_THRESH, RNG_SEED);
fprintf(fid,'Metodo: CBPT 2D en [freq x tiempo] por banda, dentro de WIN_ANALYSIS=[%.0f,%.0f]ms\n', ...
    WIN_ANALYSIS(1), WIN_ANALYSIS(2));
fprintf(fid,'Consistente con S2_TF_cluster.m (Analysis_paper) — misma alpha_voxel=0.05\n\n');
fprintf(fid,'--- CASOS (Chew > NoChew) ---\n');
for bname = {{'Theta','4-7',p_cas_theta,stat_cas_theta,win_cas_theta}, ...
             {'Alpha','8-13',p_cas_alpha,stat_cas_alpha,win_cas_alpha}, ...
             {'Beta','13-30',p_cas_beta,stat_cas_beta,win_cas_beta}}
    b = bname{1};
    fprintf(fid,'%s (%s Hz): stat=%.2f, p=%.4f', b{1},b{2},b{4},b{3});
    if ~isempty(b{5})
        fprintf(fid,' | WIN=[%.0f, %.0f] ms', b{5}(1), b{5}(2));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n--- CONTROLES ---\n');
fprintf(fid,'Theta (4-7 Hz): stat=%.2f, p=%.4f', stat_ctr_theta, p_ctr_theta);
if ~isempty(win_ctr_theta)
    fprintf(fid,' | WIN=[%.0f, %.0f] ms', win_ctr_theta(1), win_ctr_theta(2));
end
fprintf(fid,'\n\n--- WIN_CBPT CANONICO ---\n');
if ~isempty(WIN_CBPT)
    fprintf(fid,'WIN_CBPT = [%.0f, %.0f] ms\n', WIN_CBPT(1), WIN_CBPT(2));
else
    fprintf(fid,'WIN_CBPT = [] (ninguna sig.)\n');
end
fclose(fid);

fprintf('\n✓ S2a completo.\n');
fprintf('Siguiente: S2b_CBPT_electrodes.m (usa WIN_CBPT=[%.0f,%.0f]ms)\n\n', ...
    WIN_CBPT(1), WIN_CBPT(2));
