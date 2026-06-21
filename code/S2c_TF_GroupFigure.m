%% S2c_TF_GroupFigure.m  —  Figura principal TF + Topo (paper figure)
% Layout 2x3:
%   Fila 1 (topos):   [Cases Ch-Nc]  [Controls Ch-Nc]  [Interaction]
%   Fila 2 (TF maps): [Cases Ch-Nc]  [Controls Ch-Nc]  [Interaction]
%
% CBPT 2D sobre TODAS las frecuencias (1-40Hz) juntas — el cluster puede
% cruzar bandas. ALPHA_THRESH=0.05 voxel-level (consistente con V2).
%
% Contrastes:
%   Cases:       Chew - NoChew (sign-flip, N=31)
%   Controls:    Chew - NoChew (sign-flip, N=15)
%   Interaction: Cases_diff - Controls_diff (permutacion de etiquetas de grupo)
%
% Topoplots: theta Late (WIN_LATE=[900,1300]ms) por electrodo, mismos contrastes.
% Electrodos significativos: marcados con punto negro (t-test FDR o FWER).

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

ALPHA_THRESH = 0.05;   % voxel-level (consistente con S2_TF_cluster.m original)
addpath(fullfile(ROOT_P2V1, 'Analysis_paper'));   % cluster_stat_2d.m

fprintf('\n=== S2c_TF_GroupFigure.m ===\n');
fprintf('CBPT 2D fullband (1-40Hz) | alpha_voxel=%.2f | N_PERM=%d\n', ...
    ALPHA_THRESH, N_PERM);

%% ── Pool paralelo ─────────────────────────────────────────────────────────
if isempty(gcp('nocreate'))
    parpool('local', min(12, feature('numcores')));
end
% Asegurar que cluster_stat_2d.m este disponible en todos los workers,
% incluso si el pool estaba abierto antes del addpath de esta sesion.
code_dir = fullfile(ROOT_FINAL, 'code');
pctRunOnAll addpath(code_dir);

%% ── Cargar datos ──────────────────────────────────────────────────────────
load(FILE_TF, 'tf_npl_cas_ch','tf_npl_cas_nc', ...
              'tf_npl_ctr_ch','tf_npl_ctr_nc', ...
              'theta_topo_ch','theta_topo_nc', ...
              'theta_topo_ch_ctr','theta_topo_nc_ctr', ...
              'FREX_TF','times_ms');

FREX_TF  = FREX_TF(:)';
times_ms = times_ms(:)';

% Verificar dimensiones — MATLAB debe ver [N_subj x N_frex x N_times]
fprintf('Dimensiones MATLAB (sujetos x frec x tiempo):\n');
fprintf('  cas_ch: [%d x %d x %d]  — esperado [%d x 200 x 1152]\n', ...
    size(tf_npl_cas_ch,1), size(tf_npl_cas_ch,2), size(tf_npl_cas_ch,3), N_CASES);
fprintf('  ctr_ch: [%d x %d x %d]  — esperado [%d x 200 x 1152]\n', ...
    size(tf_npl_ctr_ch,1), size(tf_npl_ctr_ch,2), size(tf_npl_ctr_ch,3), N_CONTROLS);
if size(tf_npl_cas_ch,1) ~= N_CASES
    error('ERROR: dim 1 (%d) != N_CASES (%d). Verificar carga del archivo.', ...
        size(tf_npl_cas_ch,1), N_CASES);
end

%% ── Trim a ventana de análisis ────────────────────────────────────────────
tidx       = times_ms >= WIN_ANALYSIS(1) & times_ms <= WIN_ANALYSIS(2);
times_anal = times_ms(tidx);
n_times    = numel(times_anal);
n_frex     = numel(FREX_TF);

cas_ch = tf_npl_cas_ch(:, :, tidx);   % [31 x 200 x n_times]
cas_nc = tf_npl_cas_nc(:, :, tidx);
ctr_ch = tf_npl_ctr_ch(:, :, tidx);   % [15 x 200 x n_times]
ctr_nc = tf_npl_ctr_nc(:, :, tidx);

diff_cas = cas_ch - cas_nc;   % [31 x 200 x n_times]
diff_ctr = ctr_ch - ctr_nc;   % [15 x 200 x n_times]

%% ── Mapas medios (dB) ─────────────────────────────────────────────────────
map_cas = squeeze(mean(diff_cas, 1));   % [200 x n_times]
map_ctr = squeeze(mean(diff_ctr, 1));
map_int = map_cas - map_ctr;            % Interaction: diferencia de diferencias

%% ── CBPT 2D fullband — helper function ───────────────────────────────────
function [p_c, stat_obs, mask_obs, t_obs] = cbpt_fullband_signflip(data_NxFxT, N_subj, n_perm, alpha_t, rng_seed)
    rng(rng_seed, 'twister');
    [~, p_obs] = ttest(data_NxFxT);
    p_obs    = squeeze(p_obs);
    mn       = squeeze(mean(data_NxFxT));
    sd_      = squeeze(std(data_NxFxT));
    t_obs    = mn ./ (sd_ / sqrt(N_subj));
    [stat_obs, mask_obs] = cluster_stat_2d(p_obs < alpha_t & t_obs > 0, t_obs);
    signs_all = (rand(n_perm, N_subj) > 0.5)*2 - 1;
    null_dist = zeros(n_perm, 1);
    parfor pp = 1:n_perm
        d_p     = bsxfun(@times, data_NxFxT, reshape(signs_all(pp,:)', [N_subj,1,1]));
        mn_p    = squeeze(mean(d_p));
        sd_p    = squeeze(std(d_p));
        t_p     = mn_p ./ (sd_p / sqrt(N_subj));
        p_p     = 2*(1-tcdf(abs(t_p), N_subj-1));
        [cs, ~] = cluster_stat_2d(p_p < alpha_t & t_p > 0, t_p);
        null_dist(pp) = cs;
    end
    p_c = mean(null_dist >= stat_obs);
end

function [p_c, stat_obs, mask_obs, t_obs] = cbpt_fullband_interaction(diff_cas_NxFxT, diff_ctr_NxFxT, n_perm, alpha_t, rng_seed)
    % Permutacion de etiquetas de grupo (casos vs controles)
    n1 = size(diff_cas_NxFxT, 1);
    n2 = size(diff_ctr_NxFxT, 1);
    N  = n1 + n2;
    all_diffs = cat(1, diff_cas_NxFxT, diff_ctr_NxFxT);  % [N x F x T]
    rng(rng_seed, 'twister');

    % t-stat observado (Welch)
    m1 = squeeze(mean(all_diffs(1:n1,:,:)));
    m2 = squeeze(mean(all_diffs(n1+1:end,:,:)));
    s1 = squeeze(std(all_diffs(1:n1,:,:)));
    s2 = squeeze(std(all_diffs(n1+1:end,:,:)));
    se = sqrt(s1.^2/n1 + s2.^2/n2);
    t_obs  = (m1 - m2) ./ se;
    p_obs  = 2*(1 - tcdf(abs(t_obs), n1+n2-2));

    [stat_obs, mask_obs] = cluster_stat_2d(p_obs < alpha_t & t_obs > 0, t_obs);

    null_dist = zeros(n_perm, 1);
    parfor pp = 1:n_perm
        perm_idx = randperm(N);
        g1 = all_diffs(perm_idx(1:n1),:,:);
        g2 = all_diffs(perm_idx(n1+1:end),:,:);
        m1p = squeeze(mean(g1));  m2p = squeeze(mean(g2));
        s1p = squeeze(std(g1));   s2p = squeeze(std(g2));
        t_p = (m1p - m2p) ./ sqrt(s1p.^2/n1 + s2p.^2/n2);
        p_p = 2*(1 - tcdf(abs(t_p), N-2));
        [cs, ~] = cluster_stat_2d(p_p < alpha_t & t_p > 0, t_p);
        null_dist(pp) = cs;
    end
    p_c = mean(null_dist >= stat_obs);
end

%% ── Correr CBPT ───────────────────────────────────────────────────────────
fprintf('\nCBPT Cases (sign-flip, N=%d)...\n', N_CASES);
[p_cas, stat_cas, mask_cas, t_cas] = ...
    cbpt_fullband_signflip(diff_cas, N_CASES, N_PERM, ALPHA_THRESH, RNG_SEED);
fprintf('  Cases: stat=%.2f, p=%.4f\n', stat_cas, p_cas);

fprintf('CBPT Controls (sign-flip, N=%d)...\n', N_CONTROLS);
[p_ctr, stat_ctr, mask_ctr, t_ctr] = ...
    cbpt_fullband_signflip(diff_ctr, N_CONTROLS, N_PERM, ALPHA_THRESH, RNG_SEED+1);
fprintf('  Controls: stat=%.2f, p=%.4f\n', stat_ctr, p_ctr);

fprintf('CBPT Interaction (group-label perm, N=%d+%d)...\n', N_CASES, N_CONTROLS);
[p_int, stat_int, mask_int, t_int] = ...
    cbpt_fullband_interaction(diff_cas, diff_ctr, N_PERM, ALPHA_THRESH, RNG_SEED+2);
fprintf('  Interaction: stat=%.2f, p=%.4f\n', stat_int, p_int);

%% ── Topografías theta Late ────────────────────────────────────────────────
% theta_topo_ch/nc: [N x 64] — theta power en WIN_LATE por electrodo
topo_cas = median(theta_topo_ch    - theta_topo_nc,    1);   % [1 x 64]
topo_ctr = median(theta_topo_ch_ctr - theta_topo_nc_ctr, 1);
topo_int = topo_cas - topo_ctr;

% Significancia por electrodo en topoplots
% Cases: t-test por electrodo (ya en roi_cbpt.mat, pero recalculamos aquí)
[~, p_topo_cas, ~, st_cas] = ttest(theta_topo_ch - theta_topo_nc);
t_topo_cas = st_cas.tstat;

[~, p_topo_ctr, ~, st_ctr] = ttest(theta_topo_ch_ctr - theta_topo_nc_ctr);
t_topo_ctr = st_ctr.tstat;

% Interaction topo: two-sample t-test
diff_topo_cas = theta_topo_ch    - theta_topo_nc;     % [31 x 64]
diff_topo_ctr = theta_topo_ch_ctr - theta_topo_nc_ctr; % [15 x 64]
[~, p_topo_int, ~, st_int] = ttest2(diff_topo_cas, diff_topo_ctr);
t_topo_int = st_int.tstat;

% FDR correccion para topos (Benjamin-Hochberg)
fdr_cas = fdr_bh(p_topo_cas(:));   sig_cas_fdr = fdr_cas < 0.05;
fdr_ctr = fdr_bh(p_topo_ctr(:));   sig_ctr_fdr = fdr_ctr < 0.05;
fdr_int = fdr_bh(p_topo_int(:));   sig_int_fdr = fdr_int < 0.05;

fprintf('\nTopos — sig FDR: Cases=%d, Controls=%d, Interaction=%d electrodos\n', ...
    sum(sig_cas_fdr), sum(sig_ctr_fdr), sum(sig_int_fdr));

%% ── Cargar chanlocs ───────────────────────────────────────────────────────
ep_ref = fullfile(DIR_EPOCHS, [CASES{1} '_Nc_ep.set']);
if exist(ep_ref, 'file')
    EEG_ref  = pop_loadset(ep_ref);
    chanlocs = EEG_ref.chanlocs(1:EEG_N);
    all_labels = {chanlocs.labels};
else
    chanlocs = []; all_labels = {};
    fprintf('WARN: chanlocs no disponible — topos sin labels\n');
end

%% ── Colorescala compartida ────────────────────────────────────────────────
clim_tf   = 1.5;    % ±dB para TF maps
clim_topo = max(abs([topo_cas topo_ctr topo_int])) * 1.1;

%% ── Figura 2x3 ───────────────────────────────────────────────────────────
fig = figure('Name','TF+Topo Group Figure','Units','normalized',...
             'Position',[0.01 0.02 0.97 0.90], 'Color','w');

titles_topo = {sprintf('Cases Ch>Nc\n(N=%d)', N_CASES), ...
               sprintf('Controls Ch>Nc\n(N=%d)', N_CONTROLS), ...
               'Interaction\n(Cases - Controls)'};
titles_tf   = {sprintf('Cases  p=%.4f', p_cas), ...
               sprintf('Controls  p=%.4f', p_ctr), ...
               sprintf('Interaction  p=%.4f', p_int)};

topo_data  = {topo_cas, topo_ctr, topo_int};
topo_sig   = {sig_cas_fdr, sig_ctr_fdr, sig_int_fdr};
tf_maps    = {map_cas, map_ctr, map_int};
tf_masks   = {mask_cas, mask_ctr, mask_int};
tf_pvals   = {p_cas, p_ctr, p_int};

for col = 1:3
    %% --- Fila 1: Topoplots ---
    subplot(2, 3, col);
    topo_v  = topo_data{col}(:);
    sig_idx = find(topo_sig{col});   % indices de electrodos sig (FDR)
    if ~isempty(chanlocs)
        try
            % topoplot (EEGLAB) soporta emarker2 para marcar electrodos
            if ~isempty(sig_idx)
                topoplot(topo_v, chanlocs, ...
                    'electrodes', 'off', ...
                    'emarker2',   {sig_idx, '.', 'k', 18, 2});
            else
                topoplot(topo_v, chanlocs, 'electrodes', 'off');
            end
        catch ME
            % fallback si topoplot falla
            try
                topoplotIndie(topo_v, chanlocs, 'electrodes','off');
            catch
            end
            text(0.5, 0.5, sprintf('sig=%d', numel(sig_idx)), ...
                'HorizontalAlignment','center','Units','normalized','Color','k');
        end
    else
        bar(topo_v); title('sin chanlocs');
    end
    set(gca, 'CLim', [-clim_topo clim_topo]);
    colormap(gca, jet); cb = colorbar; cb.Label.String = 'θ Δ power (dB)';
    % Subtitle con labels de electrodos sig
    if ~isempty(sig_idx) && ~isempty(chanlocs)
        sig_labels = strjoin({chanlocs(sig_idx).labels}, ', ');
        title({titles_topo{col}, sprintf('sig: %s', sig_labels)}, ...
            'FontSize', FIG_FS, 'FontWeight','bold');
    else
        title(titles_topo{col}, 'FontSize', FIG_FS+1, 'FontWeight','bold');
    end
    axis square; set(gca,'FontName',FIG_FONT,'FontSize',FIG_FS);

    %% --- Fila 2: TF maps ---
    subplot(2, 3, col + 3);
    imagesc(times_anal, FREX_TF, tf_maps{col});
    axis xy; colormap(gca, jet);
    set(gca, 'CLim', [-clim_tf clim_tf]);
    cb2 = colorbar; cb2.Label.String = 'Δ power (dB)';
    hold on;

    % Contorno del cluster (blanco si sig, gris si no)
    mask = tf_masks{col};
    if ~isempty(mask) && any(mask(:)) && tf_pvals{col} < ALPHA_CLUST
        contour(times_anal, FREX_TF, double(mask), 1, 'w-', 'LineWidth', 2.5);
    end

    % Líneas de referencia de bandas
    yline(BAND_THETA(2), 'w:', 'LineWidth', 1);   % 7 Hz
    yline(BAND_ALPHA(2), 'w:', 'LineWidth', 1);   % 13 Hz
    % Líneas temporales
    xline(0,              'k--', 'LineWidth', 1.5);
    xline(WIN_EARLY(1),   'w:', 'LineWidth', 1);
    xline(WIN_LATE(1),    'w:', 'LineWidth', 1);
    xline(WIN_LATE(2),    'w:', 'LineWidth', 1);

    xlabel('Time (ms)', 'FontSize', FIG_FS);
    ylabel('Frequency (Hz)', 'FontSize', FIG_FS);
    title(titles_tf{col}, 'FontSize', FIG_FS+1, 'FontWeight','bold');
    xlim([times_anal(1) times_anal(end)]);
    ylim([FREX_TF(1) FREX_TF(end)]);
    set(gca, 'FontName', FIG_FONT, 'FontSize', FIG_FS);
end

sgtitle(sprintf(['TF 2D CBPT  |  α_{voxel}=%.2f  |  N_{PERM}=%d  |  ' ...
                 'WIN\\_ANALYSIS=[%.0f,%.0f]ms'], ...
    ALPHA_THRESH, N_PERM, WIN_ANALYSIS(1), WIN_ANALYSIS(2)), ...
    'FontSize', FIG_FS+2, 'FontWeight','bold');

saveas(fig, fullfile(OUT_FIGS, 'S2c_TF_GroupFigure.png'));
print(fig, fullfile(OUT_FIGS, 'S2c_TF_GroupFigure'), '-dpng', '-r300');
fprintf('\nFigura guardada (PNG 300dpi).\n');

%% ── Guardar estadísticos ─────────────────────────────────────────────────
save(fullfile(OUT_STATS,'S2c_TF_GroupFigure.mat'), ...
    'p_cas','stat_cas','mask_cas','t_cas','map_cas', ...
    'p_ctr','stat_ctr','mask_ctr','t_ctr','map_ctr', ...
    'p_int','stat_int','mask_int','t_int','map_int', ...
    'topo_cas','topo_ctr','topo_int', ...
    't_topo_cas','t_topo_ctr','t_topo_int', ...
    'p_topo_cas','p_topo_ctr','p_topo_int', ...
    'sig_cas_fdr','sig_ctr_fdr','sig_int_fdr', ...
    'fdr_cas','fdr_ctr','fdr_int', ...
    'times_anal','FREX_TF','ALPHA_THRESH','N_PERM');

%% ── Reporte ───────────────────────────────────────────────────────────────
fid = fopen(fullfile(OUT_REVIEWER,'S2c_TF_GroupFigure_result.txt'),'w');
fprintf(fid,'=== S2c TF Group Figure — CBPT 2D fullband ===\n');
fprintf(fid,'%s | N_PERM=%d | alpha_voxel=%.2f\n\n', datestr(now,'yyyy-mm-dd'), N_PERM, ALPHA_THRESH);
fprintf(fid,'Cases  (N=%d): cluster stat=%.2f, p=%.4f\n', N_CASES, stat_cas, p_cas);
fprintf(fid,'Controls (N=%d): cluster stat=%.2f, p=%.4f\n', N_CONTROLS, stat_ctr, p_ctr);
fprintf(fid,'Interaction: cluster stat=%.2f, p=%.4f\n\n', stat_int, p_int);
fprintf(fid,'Topos — electrodos sig (FDR q<0.05):\n');
fprintf(fid,'  Cases: %d  (%s)\n', sum(sig_cas_fdr), strjoin(all_labels(sig_cas_fdr),', '));
fprintf(fid,'  Controls: %d\n', sum(sig_ctr_fdr));
fprintf(fid,'  Interaction: %d  (%s)\n', sum(sig_int_fdr), strjoin(all_labels(sig_int_fdr),', '));
fclose(fid);

fprintf('✓ S2c completo.\n\n');

%% ── FDR Benjamini-Hochberg (inline) ──────────────────────────────────────
function p_adj = fdr_bh(p_vals)
    p_vals = p_vals(:);
    n = numel(p_vals);
    [ps, si] = sort(p_vals);
    pa = ps .* n ./ (1:n)';
    for k = n-1:-1:1, pa(k) = min(pa(k), pa(k+1)); end
    p_adj = zeros(n,1);
    p_adj(si) = pa;
end
