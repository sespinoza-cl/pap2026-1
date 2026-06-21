%% S2b_CBPT_electrodes.m  —  FASE 1 Paso 2: Electrodos significativos -> ROI_CBPT
% Requiere: S2a_CBPT_timewindow.mat (para WIN_CBPT)
%           v1_S2_TF_data.mat (theta_topo_ch, theta_topo_nc, [N x 64])
% Output: outputs/stats/roi_cbpt.mat  (ROI_LABELS, ROI_IDX, WIN_CBPT, etc.)
%         outputs/reviewer/S2b_CBPT_electrodes_result.txt
%
% Algoritmo (permutation max-t sobre electrodos):
%   1. diff_topo = theta_topo_ch - theta_topo_nc  -> [N x 64]
%      (theta_topo computada sobre WIN_LATE=[900,1300]ms en v1_S2_TF_data.mat)
%   2. t-test univariado en cada uno de los 64 electrodos -> t_obs [1 x 64]
%   3. Permutacion max-t: N_PERM inversiones de signo, max|t| nulo
%   4. Electrodos sig.: t_obs > percentil 95 de null (FWER)
%   5. Adicionalmente: FDR q<0.05 (Benjamini-Hochberg) como resultado secundario
% Nota: theta_topo fue calculada con WIN_LATE. Si WIN_CBPT difiere
%   sustancialmente, recomputar topo con topo_band_ep.m (marcado en reporte).

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file')
        run('S0_config.m');
    else
        error('S0_config.m no encontrado. Correr desde Analysis_V1_Final/');
    end
end

fprintf('\n=== S2b_CBPT_electrodes.m — FASE 1 Paso 2 ===\n');

%% ── Cargar WIN_CBPT desde S2a ────────────────────────────────────────────
s2a_file = fullfile(OUT_STATS, 'S2a_CBPT_timewindow.mat');
assert(exist(s2a_file,'file')==2, ...
    'ERROR: %s no encontrado. Correr S2a_CBPT_timewindow.m primero.', s2a_file);
load(s2a_file, 'WIN_CBPT', 'p_cas_theta');
fprintf('WIN_CBPT cargado: [%.0f, %.0f] ms (p_cas_theta=%.4f)\n', ...
    WIN_CBPT(1), WIN_CBPT(2), p_cas_theta);

% Advertencia si WIN_CBPT difiere mucho de WIN_LATE
if isempty(WIN_CBPT)
    WIN_CBPT = WIN_LATE;
    fprintf('WARN: WIN_CBPT vacia, usando WIN_LATE=[%.0f,%.0f] ms\n', WIN_LATE(1), WIN_LATE(2));
end
delta_start = abs(WIN_CBPT(1) - WIN_LATE(1));
delta_end   = abs(WIN_CBPT(2) - WIN_LATE(2));
if delta_start > 200 || delta_end > 200
    fprintf(['WARN: WIN_CBPT difiere >200ms de WIN_LATE.\n' ...
             '      theta_topo fue computada con WIN_LATE=[%.0f,%.0f]ms.\n' ...
             '      Considerar recomputar theta_topo para WIN_CBPT exacta.\n'], ...
             WIN_LATE(1), WIN_LATE(2));
end

%% ── Cargar theta_topo ────────────────────────────────────────────────────
fprintf('\nCargando theta_topo de %s ...\n', FILE_TF);
load(FILE_TF, 'theta_topo_ch','theta_topo_nc', ...
              'theta_topo_ch_ctr','theta_topo_nc_ctr');

% theta_topo_ch y _nc: [N_cases x 64]
assert(size(theta_topo_ch,1) == N_CASES, ...
    'ERROR: dim 1 de theta_topo_ch (%d) != N_CASES (%d)', size(theta_topo_ch,1), N_CASES);
assert(size(theta_topo_ch,2) == EEG_N, ...
    'ERROR: dim 2 de theta_topo_ch (%d) != EEG_N (%d)', size(theta_topo_ch,2), EEG_N);
fprintf('Dims OK: theta_topo_ch [%d x %d]\n', size(theta_topo_ch,1), size(theta_topo_ch,2));

%% ── Diferencia por sujeto y electrodo ────────────────────────────────────
diff_topo    = theta_topo_ch    - theta_topo_nc;     % [N_cases x 64]
diff_ctr     = theta_topo_ch_ctr - theta_topo_nc_ctr; % [N_controls x 64]

% Estadistico observado: t-test en cada electrodo (Chew-NoChew > 0)
[~, p_obs, ~, stats_obs] = ttest(diff_topo);  % [1 x 64]
t_obs = stats_obs.tstat;                       % [1 x 64]

[~, p_ctr, ~, stats_ctr] = ttest(diff_ctr);
t_ctr = stats_ctr.tstat;

fprintf('\nEstadisticos observados:\n');
fprintf('  Casos:     t max=%.2f (elec %d)\n', max(t_obs), find(t_obs==max(t_obs),1));
fprintf('  Controles: t max=%.2f (elec %d)\n', max(t_ctr), find(t_ctr==max(t_ctr),1));

%% ── Permutation max-t (FWER) ────────────────────────────────────────────
fprintf('\nPermutation max-t (N_PERM=%d, seed=%d)...\n', N_PERM, RNG_SEED);
rng(RNG_SEED, 'twister');
signs_mat  = (rand(N_PERM, N_CASES) > 0.5)*2 - 1;   % [N_PERM x N_cases]
null_max_t = zeros(N_PERM, 1);

for pp = 1:N_PERM
    diff_perm  = diff_topo .* signs_mat(pp,:)';      % [N x 64]
    [~, ~, ~, s_p] = ttest(diff_perm);
    null_max_t(pp) = max(s_p.tstat);
end

t_thresh_fwer = prctile(null_max_t, 95);             % umbral para p_fwer < 0.05
sig_fwer = t_obs > t_thresh_fwer;                    % [1 x 64] logical

% p-valor por electrodo vs distribucion nula
p_fwer = zeros(1, EEG_N);
for e = 1:EEG_N
    p_fwer(e) = mean(null_max_t >= t_obs(e));
end

fprintf('Umbral max-t FWER: t > %.3f (95th percentil nulo)\n', t_thresh_fwer);
fprintf('Electrodos sig. FWER (p<0.05): %d / %d\n', sum(sig_fwer), EEG_N);

%% ── FDR (Benjamini-Hochberg) como resultado secundario ───────────────────
% FDR sobre p-valores del t-test observado
[p_fdr_adj] = fdr_bh(p_obs);
sig_fdr = p_fdr_adj < 0.05;

fprintf('Electrodos sig. FDR  (q<0.05): %d / %d\n', sum(sig_fdr), EEG_N);

%% ── Obtener labels de electrodos ─────────────────────────────────────────
% Necesita cargar chanlocs desde una epoca de referencia
epoch_ref = fullfile(DIR_EPOCHS, [CASES{1} '_Nc_ep.set']);
if exist(epoch_ref, 'file')
    EEG_ref = pop_loadset(epoch_ref);
    chanlocs = EEG_ref.chanlocs(1:EEG_N);
    all_labels = {chanlocs.labels};
    fprintf('Chanlocs cargados desde: %s\n', epoch_ref);
else
    % Fallback: etiquetas genericas
    all_labels = arrayfun(@(k) sprintf('Ch%d',k), 1:EEG_N, 'UniformOutput',false);
    chanlocs = [];
    fprintf('WARN: %s no encontrado. Usando labels genericos.\n', epoch_ref);
end

%% ── Definir ROI_CBPT ─────────────────────────────────────────────────────
% Decision jerarquica: FWER primero, FDR como fallback
if sum(sig_fwer) >= 1
    roi_mask    = sig_fwer;
    roi_method  = 'FWER-permutation (max-t)';
    fprintf('\nROI_CBPT definido por FWER (metodo primario)\n');
elseif sum(sig_fdr) >= 1
    roi_mask    = sig_fdr;
    roi_method  = 'FDR (BH, q<0.05) — FWER no significativo';
    fprintf('\nWARN: FWER n.s., usando FDR como fallback para ROI_CBPT\n');
else
    % Ninguno significativo: top-5 electrodos por t-stat
    [~, sort_idx] = sort(t_obs,'descend');
    roi_mask = false(1, EEG_N);
    roi_mask(sort_idx(1:5)) = true;
    roi_method = 'TOP-5 (ninguno sig. FWER ni FDR — resultado exploratorio)';
    fprintf('\nWARN: ninguna correccion sig. ROI_CBPT = top-5 electrodos.\n');
end

ROI_IDX    = find(roi_mask);
ROI_LABELS = all_labels(roi_mask);

fprintf('ROI_CBPT (%s):\n', roi_method);
fprintf('  Electrodes (%d): %s\n', numel(ROI_IDX), strjoin(ROI_LABELS,', '));
for k = 1:numel(ROI_IDX)
    fprintf('    %s: t=%.3f, p_fwer=%.4f, p_fdr=%.4f\n', ...
        ROI_LABELS{k}, t_obs(ROI_IDX(k)), p_fwer(ROI_IDX(k)), p_fdr_adj(ROI_IDX(k)));
end

%% ── Guardar roi_cbpt.mat (usado por todos los scripts de Fases 2-6) ──────
roi_file = fullfile(OUT_STATS, 'roi_cbpt.mat');
save(roi_file, ...
    'ROI_IDX','ROI_LABELS','roi_method', ...
    'WIN_CBPT','BAND_THETA', ...
    't_obs','p_obs','p_fwer','p_fdr_adj','sig_fwer','sig_fdr', ...
    't_ctr','p_ctr', ...
    'null_max_t','t_thresh_fwer', ...
    'N_CASES','N_CONTROLS','N_PERM','ALPHA_VOXEL','RNG_SEED');
fprintf('\nGuardado roi_cbpt.mat: %s\n', roi_file);

% Guardar tambien el archivo completo
out_file = fullfile(OUT_STATS, 'S2b_CBPT_electrodes.mat');
save(out_file, ...
    'diff_topo','diff_ctr','t_obs','t_ctr','p_obs','p_ctr', ...
    'p_fwer','p_fdr_adj','sig_fwer','sig_fdr', ...
    'ROI_IDX','ROI_LABELS','WIN_CBPT','roi_method', ...
    'null_max_t','t_thresh_fwer');
fprintf('Guardado: %s\n', out_file);

%% ── Figura ───────────────────────────────────────────────────────────────
fig = figure('Name','S2b CBPT Electrodos','Units','normalized',...
             'Position',[0.02 0.05 0.85 0.75]);

% Panel 1: Distribucion nula max-t vs observado
subplot(1,3,1);
histogram(null_max_t, 40, 'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');
hold on;
xline(t_thresh_fwer, 'r--', sprintf('p<0.05 (t>%.2f)',t_thresh_fwer), 'LineWidth',2);
for e = ROI_IDX
    xline(t_obs(e), 'b-', 'LineWidth', 1.5);
end
xlabel('Max t-stat (null)'); ylabel('Count');
title(sprintf('Null distribution max-t\nN_PERM=%d', N_PERM));
box off;

% Panel 2: t-stats por electrodo con umbral
subplot(1,3,2);
bar(1:EEG_N, t_obs, 'FaceColor',[0.7 0.7 0.9]);
hold on;
bar(ROI_IDX, t_obs(ROI_IDX), 'FaceColor', COL_SIG);
yline(t_thresh_fwer, 'r--', 'FWER', 'LineWidth', 2);
yline(0, 'k-');
xlabel('Electrode index'); ylabel('t-statistic');
title(sprintf('t-stat por electrodo\nROI=%s', strjoin(ROI_LABELS,', ')));
xlim([0 EEG_N+1]);
box off;

% Panel 3: Topografia si tenemos chanlocs
subplot(1,3,3);
if ~isempty(chanlocs)
    try
        topoplotIndie(t_obs, chanlocs);
        title(sprintf('Topografia t-stat\n(ROI=rojo, n=%d)', numel(ROI_IDX)));
        colorbar;
    catch ME
        fprintf('WARN topoplot: %s\n', ME.message);
        text(0.5,0.5,sprintf('chanlocs OK\n%d sig', numel(ROI_IDX)), ...
            'HorizontalAlignment','center');
    end
else
    text(0.5,0.5,'chanlocs no disponibles','HorizontalAlignment','center');
end
axis square;

sgtitle(sprintf('CBPT Electrodos — Fase 1 Paso 2 (N=%d casos)', N_CASES), 'FontWeight','bold');
saveas(fig, fullfile(OUT_FIGS,'S2b_CBPT_electrodes.png'));
fprintf('Figura guardada.\n');

%% ── Reporte para el revisor ───────────────────────────────────────────────
fid = fopen(fullfile(OUT_REVIEWER, 'S2b_CBPT_electrodes_result.txt'), 'w');
fprintf(fid,'=== S2b CBPT Electrodos — P2V1 Revision ===\n');
fprintf(fid,'Fecha: %s\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf(fid,'N_CASES=%d | N_PERM=%d | seed=%d\n\n', N_CASES, N_PERM, RNG_SEED);
fprintf(fid,'WIN_CBPT = [%.0f, %.0f] ms\n', WIN_CBPT(1), WIN_CBPT(2));
fprintf(fid,'theta_topo computada con WIN_LATE=[%.0f,%.0f]ms\n\n', WIN_LATE(1), WIN_LATE(2));
fprintf(fid,'--- CASOS ---\n');
fprintf(fid,'t_thresh_FWER = %.3f | Electrodos FWER: %d | FDR: %d\n', ...
    t_thresh_fwer, sum(sig_fwer), sum(sig_fdr));
fprintf(fid,'ROI_CBPT (%s):\n', roi_method);
fprintf(fid,'  N electrodos: %d\n', numel(ROI_IDX));
fprintf(fid,'  Labels: %s\n', strjoin(ROI_LABELS,', '));
for k=1:numel(ROI_IDX)
    fprintf(fid,'  %s: t=%.3f, p_fwer=%.4f, p_fdr=%.4f\n', ...
        ROI_LABELS{k}, t_obs(ROI_IDX(k)), p_fwer(ROI_IDX(k)), p_fdr_adj(ROI_IDX(k)));
end
fprintf(fid,'\n--- CONTROLES ---\n');
fprintf(fid,'t max=%.3f | ninguno esperado sig.\n', max(t_ctr));
fclose(fid);

fprintf('\n✓ S2b completo. ROI_CBPT: %d electrodos (%s)\n', ...
    numel(ROI_IDX), strjoin(ROI_LABELS,', '));
fprintf('Siguiente: S3_BandPower_ROI.m\n\n');

%% ── Funcion auxiliar: FDR Benjamini-Hochberg ─────────────────────────────
function p_adj = fdr_bh(p_vals)
    % Retorna p-valores ajustados por FDR (BH). Input/output: [1 x N] o [N x 1]
    p_vals = p_vals(:);
    n = numel(p_vals);
    [p_sorted, sort_idx] = sort(p_vals);
    p_adj_sorted = p_sorted .* n ./ (1:n)';
    % Asegurar monotonia
    for k = n-1:-1:1
        p_adj_sorted(k) = min(p_adj_sorted(k), p_adj_sorted(k+1));
    end
    p_adj = zeros(n, 1);
    p_adj(sort_idx) = p_adj_sorted;
    p_adj = p_adj';
end
