%% ============================================================================
%  S2_TF_plot.m - PIPELINE ELECTROFISIOLÓGICO (ERSP + ITPC)
%  Incluye: 6 TFs, 6 Topos (con electrodos sig por grupo), Ortogonalización EMG,
%           Barras Multi-banda, 2 LMEs y Correlaciones Cerebro-Conducta.
%
%  CAMBIO v2: Topoplots individuales (Cases / Controls) ahora marcan sus
%             propios electrodos significativos (one-sample t-test vs 0,
%             p < 0.05 sin corrección — estándar de visualización en EEG).
% ============================================================================
clear; clc; close all;

% Configuración de archivo de reporte
fid_stats = 1;

% Colores Institucionales y Estética
color_base     = [0.6 0.6 0.6];
color_controls = [0.55, 0.63, 0.80];
color_cases    = [0.35, 0.65, 0.55];
colors_groups  = {color_controls, color_cases};
sz_title = 16; sz_label = 14; sz_tick  = 13; sz_text = 12;
line_w = 1.5; f_name = 'Arial';
set(0,'DefaultAxesFontName',f_name,'DefaultAxesFontSize',sz_tick,'DefaultFigureColor','w');

% Parámetros TF, Permutación y Auditoría
time2 = [300 900];
freq2 = [4 7];
freq_emg = [70 90];
p_thresh      = 0.05;
n_perm        = 5000;
head_r        = 0.5;

% Rutas — fuente única (S0_paths garantiza coherencia con el resto del pipeline)
P         = S0_paths();
f_tf      = [P.dir_tf '\'];           % C:\...\Desktop\Exp2\EEG\TF\ (sólo lectura)
f_epoch   = [P.dir_epoch '\'];        % C:\...\Desktop\Exp2\EEG\Epochs\
f_beh     = P.file_beh;               % D:\Exp2\Version_Mayo\Pipelines\data_beh_tb_45.mat (N=45)
save_path = P.fig02;                  % D:\Exp2\Version_Mayo\Plots\Paper\Figure02_TF
eeglab_path = P.eeglab_path;
addpath(eeglab_path);
% S0_paths ya crea la carpeta; verificar por si se llama sin ella
if ~exist(save_path, 'dir'), mkdir(save_path); end

%% ======================================================================
%  1. CARGA DE MATRICES 6D Y EXTRACCIÓN (PWR vs ITPC)
% ======================================================================
cprintf(fid_stats, '--- 1. CARGA DE MATRICES 6D ---\n');
temp_file = dir(fullfile(f_epoch, '*_Ch.set'));
if isempty(temp_file), temp_file = dir(fullfile(f_epoch, '*_Nc.set')); end
EEG = pop_loadset('filename', temp_file(1).name, 'filepath', f_epoch);

s_b1 = load([f_tf 'Group_Cases_Nc_tfraw.mat']);
s_b2 = load([f_tf 'Group_Cases_Ch_tfraw.mat']);
c_b1 = load([f_tf 'Group_Controls_Nc_tfraw.mat']);
c_b2 = load([f_tf 'Group_Controls_Ch_tfraw.mat']);

% EXTRACCIÓN A: ERSP POWER (Inducido puro: Dim1=2, Dim6=1)
casos_b1_pwr = squeeze(s_b1.tfraw_pre_g(2,:,:,:,:,1));
casos_b2_pwr = squeeze(s_b2.tfraw_pre_g(2,:,:,:,:,1));
ctrl_b1_pwr  = squeeze(c_b1.tfraw_pre_g(2,:,:,:,:,1));
ctrl_b2_pwr  = squeeze(c_b2.tfraw_pre_g(2,:,:,:,:,1));

% EXTRACCIÓN B: ITPC (Total Phase: Dim1=1, Dim6=2)
casos_b1_itpc = squeeze(s_b1.tfraw_pre_g(1,:,:,:,:,2));
casos_b2_itpc = squeeze(s_b2.tfraw_pre_g(1,:,:,:,:,2));
ctrl_b1_itpc  = squeeze(c_b1.tfraw_pre_g(1,:,:,:,:,2));
ctrl_b2_itpc  = squeeze(c_b2.tfraw_pre_g(1,:,:,:,:,2));

% Metadata
channels = {'AFz','F1','Fz','F2','FC1','FCz','FC2'};
[~,electrodes] = ismember(channels,{EEG.chanlocs.labels});
times = s_b1.eeg_times; frex = s_b1.frex;
tidx2 = dsearchn(times', time2'); fidx2 = dsearchn(frex',  freq2');

nCases = size(casos_b1_pwr,4); nCtrl = size(ctrl_b1_pwr, 4); nElec = numel(EEG.chanlocs(1:64));

% Diferencias y ROIs
d_cases_pwr = casos_b2_pwr - casos_b1_pwr; d_ctrl_pwr = ctrl_b2_pwr - ctrl_b1_pwr;
tf_casos_roi_pwr = squeeze(mean(d_cases_pwr(electrodes,:,:,:),1));
tf_ctrl_roi_pwr  = squeeze(mean(d_ctrl_pwr( electrodes,:,:,:),1));

d_cases_itpc = casos_b2_itpc - casos_b1_itpc; d_ctrl_itpc = ctrl_b2_itpc - ctrl_b1_itpc;
tf_casos_roi_itpc = squeeze(mean(d_cases_itpc(electrodes,:,:,:),1));
tf_ctrl_roi_itpc  = squeeze(mean(d_ctrl_itpc( electrodes,:,:,:),1));

cprintf(fid_stats, '  Datos cargados: %d Casos, %d Controles.\n', nCases, nCtrl);

%% ======================================================================
%  2. TEST DE PERMUTACIONES CLUSTER-BASED (POWER & ITPC)
% ======================================================================
cprintf(fid_stats, '\n--- 2. PERMUTACIONES CLUSTER-CORRECTED ---\n');

% FUNCIÓN ANÓNIMA PARA PERMUTAR (Ahorra cientos de líneas)
run_perm = @(d1, d2, is_paired) perm_cluster(d1, d2, is_paired, n_perm, p_thresh);

% --- PERMUTACIONES POWER ---
chk_pwr = fullfile(save_path,'chk_all_pwr_clean.mat');
if isfile(chk_pwr), load(chk_pwr); else
    cprintf(fid_stats, '  Calculando Permutaciones POWER...\n');
    [sig_m_cas_pwr, th_cas_pwr] = run_perm(tf_casos_roi_pwr, [], true);
    [sig_m_ctr_pwr, th_ctr_pwr] = run_perm(tf_ctrl_roi_pwr, [], true);
    [sig_m_int_pwr, th_int_pwr] = run_perm(tf_casos_roi_pwr, tf_ctrl_roi_pwr, false);
    save(chk_pwr, 'sig_m_cas_pwr','sig_m_ctr_pwr','sig_m_int_pwr');
end

% --- PERMUTACIONES ITPC ---
chk_itpc = fullfile(save_path,'chk_all_itpc_clean.mat');
if isfile(chk_itpc), load(chk_itpc); else
    cprintf(fid_stats, '  Calculando Permutaciones ITPC...\n');
    [sig_m_cas_itpc, th_cas_itpc] = run_perm(tf_casos_roi_itpc, [], true);
    [sig_m_ctr_itpc, th_ctr_itpc] = run_perm(tf_ctrl_roi_itpc, [], true);
    [sig_m_int_itpc, th_int_itpc] = run_perm(tf_casos_roi_itpc, tf_ctrl_roi_itpc, false);
    save(chk_itpc, 'sig_m_cas_itpc','sig_m_ctr_itpc','sig_m_int_itpc');
end

%% ======================================================================
%  3. PLOTEO DE MAPAS TF (6 GRÁFICOS)
% ======================================================================
cprintf(fid_stats, '--- 3. Ploteando Mapas TF...\n');
lim_pwr = [-2 2]; lim_itpc = [-0.15 0.15];

% Función local para graficar TF
plot_tf = @(data, mask, lims, lbl, tit, fn) plot_tf_map(data, mask, times, frex, lims, lbl, tit, time2, freq2, fullfile(save_path, fn));

% MAPAS POWER
plot_tf(mean(tf_casos_roi_pwr,3), sig_m_cas_pwr, lim_pwr, 'Power (dB)', 'Cases (Chew-Nochew)', 'TF_Cases_PWR.png');
plot_tf(mean(tf_ctrl_roi_pwr,3), sig_m_ctr_pwr, lim_pwr, 'Power (dB)', 'Controls (Chew-Nochew)', 'TF_Controls_PWR.png');
plot_tf(mean(tf_casos_roi_pwr,3) - mean(tf_ctrl_roi_pwr,3), sig_m_int_pwr, lim_pwr, 'Interaction (dB)', 'Interaction (Cases-Ctrl)', 'TF_Interact_PWR.png');

% MAPAS ITPC
plot_tf(mean(tf_casos_roi_itpc,3), sig_m_cas_itpc, lim_itpc, 'ITPC (\Delta)', 'Cases (Chew-Nochew)', 'TF_Cases_ITPC.png');
plot_tf(mean(tf_ctrl_roi_itpc,3), sig_m_ctr_itpc, lim_itpc, 'ITPC (\Delta)', 'Controls (Chew-Nochew)', 'TF_Controls_ITPC.png');
plot_tf(mean(tf_casos_roi_itpc,3) - mean(tf_ctrl_roi_itpc,3), sig_m_int_itpc, lim_itpc, 'Interaction (\Delta)', 'Interaction (Cases-Ctrl)', 'TF_Interact_ITPC.png');

%% ======================================================================
%  4. PLOTEO DE TOPOPLOTS (6 GRÁFICOS)
%
%  CAMBIO: Cada topo individual (Cases / Controls) ahora muestra sus
%  propios electrodos significativos calculados con one-sample t-test vs 0
%  (p < 0.05, uncorrected — estándar de visualización en literatura EEG).
%  La interacción sigue usando el two-sample t-test entre grupos.
% ======================================================================
cprintf(fid_stats, '--- 4. Calculando Topoplots...\n');
lim_t_pwr = [-1 1]; lim_t_itpc = [-0.08 0.08];

% Función de extracción de ventana (theta early por defecto: time2, freq2)
% t_cas_pwr / t_ctr_pwr → shape (nElec, nSubj)
get_topo = @(d) squeeze(mean(mean(d(:,fidx2(1):fidx2(2),tidx2(1):tidx2(2),:),2),3));
t_cas_pwr  = get_topo(d_cases_pwr);  t_ctr_pwr  = get_topo(d_ctrl_pwr);
t_cas_itpc = get_topo(d_cases_itpc); t_ctr_itpc = get_topo(d_ctrl_itpc);

% ---- Electrodos sig: INTERACCIÓN (two-sample t-test, uncorrected) --------
sig_elec_pwr  = find(abs(tstat2(t_cas_pwr,  t_ctr_pwr))  > abs(tinv(0.05/2, nCases+nCtrl-2)));
sig_elec_itpc = find(abs(tstat2(t_cas_itpc, t_ctr_itpc)) > abs(tinv(0.05/2, nCases+nCtrl-2)));

% ---- Electrodos sig: CASOS vs cero (one-sample t-test por electrodo) -----
% t_cas_pwr es (nElec x nCases); ttest a lo largo de dim 2 (sujetos)
[~,~,~,st_cas_pwr]  = ttest(t_cas_pwr,  0, 'dim', 2);
[~,~,~,st_cas_itpc] = ttest(t_cas_itpc, 0, 'dim', 2);
sig_cas_pwr  = find(abs(st_cas_pwr.tstat)  > abs(tinv(0.05/2, nCases-1)));
sig_cas_itpc = find(abs(st_cas_itpc.tstat) > abs(tinv(0.05/2, nCases-1)));

% ---- Electrodos sig: CONTROLES vs cero (one-sample t-test por electrodo) -
[~,~,~,st_ctr_pwr]  = ttest(t_ctr_pwr,  0, 'dim', 2);
[~,~,~,st_ctr_itpc] = ttest(t_ctr_itpc, 0, 'dim', 2);
sig_ctr_pwr  = find(abs(st_ctr_pwr.tstat)  > abs(tinv(0.05/2, nCtrl-1)));
sig_ctr_itpc = find(abs(st_ctr_itpc.tstat) > abs(tinv(0.05/2, nCtrl-1)));

% ---- Reporte de electrodos significativos --------------------------------
cprintf(fid_stats, '  POWER — Electrodos sig (p<0.05 uncorr):\n');
cprintf(fid_stats, '    Cases:       %d / %d electrodes\n', numel(sig_cas_pwr),  nElec);
cprintf(fid_stats, '    Controls:    %d / %d electrodes\n', numel(sig_ctr_pwr),  nElec);
cprintf(fid_stats, '    Interaction: %d / %d electrodes\n', numel(sig_elec_pwr), nElec);
cprintf(fid_stats, '  ITPC — Electrodos sig (p<0.05 uncorr):\n');
cprintf(fid_stats, '    Cases:       %d / %d electrodes\n', numel(sig_cas_itpc),  nElec);
cprintf(fid_stats, '    Controls:    %d / %d electrodes\n', numel(sig_ctr_itpc),  nElec);
cprintf(fid_stats, '    Interaction: %d / %d electrodes\n', numel(sig_elec_itpc), nElec);

% ---- Función local para Topoplots ----------------------------------------
plot_topo = @(data, sigE, lims, lbl, tit, fn) ...
    plot_topo_map(data, sigE, EEG.chanlocs(1:64), head_r, lims, lbl, tit, fullfile(save_path, fn));

% TOPO POWER
plot_topo(mean(t_cas_pwr,2),  sig_cas_pwr,  lim_t_pwr, 'Power (dB)',    'Cases PWR',       'Topo_Cases_PWR.png');
plot_topo(mean(t_ctr_pwr,2),  sig_ctr_pwr,  lim_t_pwr, 'Power (dB)',    'Controls PWR',    'Topo_Controls_PWR.png');
plot_topo(mean(t_cas_pwr,2)-mean(t_ctr_pwr,2), sig_elec_pwr, lim_t_pwr, 'Int. Power (dB)', 'Interaction PWR', 'Topo_Interact_PWR.png');

% TOPO ITPC
plot_topo(mean(t_cas_itpc,2), sig_cas_itpc, lim_t_itpc, 'ITPC (\Delta)', 'Cases ITPC',      'Topo_Cases_ITPC.png');
plot_topo(mean(t_ctr_itpc,2), sig_ctr_itpc, lim_t_itpc, 'ITPC (\Delta)', 'Controls ITPC',   'Topo_Controls_ITPC.png');
plot_topo(mean(t_cas_itpc,2)-mean(t_ctr_itpc,2), sig_elec_itpc, lim_t_itpc, 'Int. ITPC', 'Interaction ITPC', 'Topo_Interact_ITPC.png');

%% ======================================================================
%  5. BLINDAJE EMG (ORTOGONALIZACIÓN ESPACIAL) - SOLO POWER
% ======================================================================
cprintf(fid_stats, '\n--- 5. BLINDAJE EMG (POWER) ---\n');
fidx_emg = dsearchn(frex', freq_emg');
topo_emg_cases = squeeze(mean(d_cases_pwr(:, fidx_emg(1):fidx_emg(2), tidx2(1):tidx2(2), :), [2 3 4]));
topo_emg_ctrl  = squeeze(mean(d_ctrl_pwr(:,  fidx_emg(1):fidx_emg(2), tidx2(1):tidx2(2), :), [2 3 4]));

interact_emg   = double(topo_emg_cases - topo_emg_ctrl);
interact_theta = double(mean(t_cas_pwr,2) - mean(t_ctr_pwr,2));

p_reg = polyfit(interact_emg, interact_theta, 1);
theta_pred = polyval(p_reg, interact_emg);
interact_theta_residual = interact_theta - theta_pred;

[r_emg, p_val_emg] = corr(interact_theta_residual, interact_emg);
cprintf(fid_stats, 'θ Residual vs EMG: r = %.3f (p = %.3f)\n', r_emg, p_val_emg);

figEMG = figure('Color','w','Position',[200 200 850 450]);

% Usar tiledlayout en lugar de subplot para evitar deformaciones
t = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Gráfico 1 ---
ax1 = nexttile;
topoplotIndie(interact_theta, EEG.chanlocs(1:64), 'plotrad',0.6, 'shading','flat', 'electrodes','on');
colormap(ax1, jet); clim(ax1, lim_t_pwr);
title(ax1, 'Original (Theta)', 'FontSize', sz_title);

% --- Gráfico 2 ---
ax2 = nexttile;
topoplotIndie(interact_theta_residual, EEG.chanlocs(1:64), 'plotrad',0.6, 'shading','flat', 'electrodes','on');
colormap(ax2, jet); clim(ax2, lim_t_pwr);
title(ax2, 'Neural Effect (Residual)', 'FontSize', sz_title, 'FontWeight', 'bold');

% --- Colorbar Global ---
cb = colorbar(ax2);
cb.Layout.Tile = 'east';
ylabel(cb, 'Residual Power (dB)', 'FontSize', sz_label);

exportgraphics(figEMG, fullfile(save_path, 'Topo_Interact_Residual_EMG.png'), 'Resolution', 300); close(figEMG);

%% ======================================================================
%  6. BARRAS UNILATERALES, ORIGEN CERO Y EJES DINÁMICOS (TF-GUIADAS)
% ======================================================================
cprintf(fid_stats, '\n--- 6. COMPARACIONES FOCALES Y PLOTEOS UNILATERALES ---\n');

% Definición estricta de ventanas: {Nombre, Freq_Hz, Time_ms, Nombre_Archivo}
% Bandas desde S0_paths (θ4-7 / α8-12 / β13-30); ventanas Early/Late/Active
bTh = P.bands_hz{1};  bAl = P.bands_hz{2};  bBe = P.bands_hz{3};
sTh = sprintf('Theta (%d-%d Hz)', bTh);
sAl = sprintf('Alpha (%d-%d Hz)', bAl);
sBe = sprintf('Beta (%d-%d Hz)',  bBe);
tf_rois = {
    sTh, bTh, [0 300],   'Theta_WinEarly';
    sTh, bTh, [300 900], 'Theta_WinLate';
    sTh, bTh, [200 700], 'Theta_WinActive';

    sAl, bAl, [0 300],   'Alpha_WinEarly';
    sAl, bAl, [300 900], 'Alpha_WinLate';
    sAl, bAl, [200 700], 'Alpha_WinActive';

    sBe, bBe, [0 300],   'Beta_WinEarly';
    sBe, bBe, [300 900], 'Beta_WinLate';
    sBe, bBe, [200 700], 'Beta_WinActive'
    };

% Paleta de colores exacta
color_cases_ch = [0.35, 0.65, 0.55]; % Verde oscuro (Cases Chew)
color_cases_nc = [0.70, 0.85, 0.80]; % Verde claro (Cases No-Chew)
color_ctrls_ch = [0.55, 0.63, 0.80]; % Azul oscuro (Controls Chew)
color_ctrls_nc = [0.80, 0.85, 0.92]; % Azul claro (Controls No-Chew)
colors_ch = {color_ctrls_ch, color_cases_ch};
colors_nc = {color_ctrls_nc, color_cases_nc};

% Función anónima para colapsar espacio, frecuencia y tiempo
get_pwr = @(data, t, f) squeeze(mean(data(electrodes, f(1):f(2), t(1):t(2), :), [1 2 3]));

for r = 1:size(tf_rois, 1)
    roi_label = tf_rois{r,1};
    f_hz      = tf_rois{r,2};
    t_ms      = tf_rois{r,3};
    plot_name = tf_rois{r,4};

    % Etiqueta dinámica para consola
    if t_ms(1) == 0 && t_ms(2) == 300
        win_label = 'Early (0-300ms)';
    elseif t_ms(1) == 300 && t_ms(2) == 900
        win_label = 'Late (300-900ms)';
    else
        win_label = sprintf('Active (%d-%dms)', t_ms(1), t_ms(2));
    end

    roi_full_name = sprintf('%s - %s', roi_label, win_label);
    fidx = dsearchn(frex', f_hz');
    tidx = dsearchn(times', t_ms');

    c1 = get_pwr(ctrl_b1_pwr, tidx, fidx); c2 = get_pwr(ctrl_b2_pwr, tidx, fidx);
    s1 = get_pwr(casos_b1_pwr, tidx, fidx); s2 = get_pwr(casos_b2_pwr, tidx, fidx);

    plot_data_bars = {{c1, c2}, {s1, s2}};
    groups = {'Controls', 'Cases'};
    cprintf(fid_stats, '\n >> ROI: %s\n', roi_full_name);

    for g = 1:2
        d = plot_data_bars{g};

        % 1. ESTADÍSTICA (Wilcoxon Signed-Rank)
        [p_val, ~, stats_w] = signrank(d{1}, d{2}, 'tail', 'both');
        if isfield(stats_w, 'zval')
            stat_str = sprintf('Z = %.3f', stats_w.zval);
        else
            stat_str = sprintf('W = %.1f', stats_w.signedrank);
        end
        cprintf(fid_stats, '  [%s] %s: p = %.4f (%s)\n', plot_name, groups{g}, p_val, stat_str);

        % 2. GENERACIÓN DEL GRÁFICO (Unilateral & Dynamic Limits)
        fig_bar = figure('Color','w','Position',[100 100 350 400], 'Visible', 'off');
        hold on;

        % Media y SEM
        m1 = nanmean(d{1}); m2 = nanmean(d{2});
        se1 = nanstd(d{1})/sqrt(length(d{1})); se2 = nanstd(d{2})/sqrt(length(d{2}));
        means = [m1, m2]; sems = [se1, se2];

        % Dibujar Barras (parten de 0 por defecto)
        b1 = bar(1, m1, 0.6, 'FaceColor', colors_nc{g}, 'EdgeColor', 'none');
        b2 = bar(2, m2, 0.6, 'FaceColor', colors_ch{g}, 'EdgeColor', 'none');

        % Dibujar Barras de Error UNILATERALES proyeccion "hacia afuera"
        if m1 >= 0
            errorbar(1, m1, 0, se1, 'Color', colors_nc{g}*0.7, 'LineStyle', 'none', 'LineWidth', 2.5, 'CapSize', 8);
        else
            errorbar(1, m1, se1, 0, 'Color', colors_nc{g}*0.7, 'LineStyle', 'none', 'LineWidth', 2.5, 'CapSize', 8);
        end

        if m2 >= 0
            errorbar(2, m2, 0, se2, 'Color', colors_ch{g}*0.7, 'LineStyle', 'none', 'LineWidth', 2.5, 'CapSize', 8);
        else
            errorbar(2, m2, se2, 0, 'Color', colors_ch{g}*0.7, 'LineStyle', 'none', 'LineWidth', 2.5, 'CapSize', 8);
        end

        max_abs_val = 1.5;
        current_lims = [-max_abs_val, max_abs_val];

        ylim(gca, current_lims);
        xlim(gca, [0.3 2.7]);
        yline(gca, 0, 'k-', 'LineWidth', 1.0);

        set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 18, 'LineWidth', 1.2, ...
            'XTick', [1 2], 'XTickLabel', {'No-Chew', 'Chew'}, 'YGrid', 'on', 'XGrid', 'off', 'GridAlpha', 0.15);
        ylabel('Power (dB)', 'FontWeight', 'bold', 'FontSize', 22);

        if p_val < 0.001, sig_txt = '***';
        elseif p_val < 0.01, sig_txt = '**';
        elseif p_val < 0.05, sig_txt = '*';
        else, sig_txt = 'n.s.';
        end

        if mean(means) >= 0
            y_pos = current_lims(2) * 0.88;
            v_align = 'top';
        else
            y_pos = current_lims(1) * 0.88;
            v_align = 'bottom';
        end
        text(1.5, y_pos, sig_txt, 'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', v_align);

        % 3. GUARDAR GRÁFICO
        fn_out = sprintf('Bar_Band%s_Group%s.png', plot_name, groups{g});
        exportgraphics(fig_bar, fullfile(save_path, fn_out), 'Resolution', 300);
        close(fig_bar);
    end
end

cprintf(fid_stats, '\n  -> %d Gráficos de barras unilaterales exportados.\n', size(tf_rois,1)*2);

%% ======================================================================
%  7. MODELOS DE EFECTOS MIXTOS (LME) CON CONTROL EMG ORTOGONAL
% ======================================================================
fprintf('\n--- 7. LME - ANÁLISIS ROBUSTO CON COVARIABLE EMG ---\n');

% Preparar vectores de diseño
group_labels = [repmat({'Control'}, nCtrl*2, 1); repmat({'Cases'}, nCases*2, 1)];
block_labels = [repmat({'B1'}, nCtrl, 1); repmat({'B2'}, nCtrl, 1); repmat({'B1'}, nCases, 1); repmat({'B2'}, nCases, 1)];
subj_ids = [ (1:nCtrl)'; (1:nCtrl)'; ( (1:nCases) + 100 )'; ( (1:nCases) + 100 )' ];

freq_emg_high = [70 90];
fidx_emg = dsearchn(frex', freq_emg_high');

for r = 1:size(tf_rois, 1)
    roi_name = tf_rois{r,1};
    f_hz = tf_rois{r,2};
    t_ms = tf_rois{r,3};

    fidx = dsearchn(frex', f_hz');
    tidx = dsearchn(times', t_ms');

    all_power = [get_pwr(ctrl_b1_pwr, tidx, fidx); get_pwr(ctrl_b2_pwr, tidx, fidx); get_pwr(casos_b1_pwr, tidx, fidx); get_pwr(casos_b2_pwr, tidx, fidx)];
    all_itpc  = [get_pwr(ctrl_b1_itpc, tidx, fidx); get_pwr(ctrl_b2_itpc, tidx, fidx); get_pwr(casos_b1_itpc, tidx, fidx); get_pwr(casos_b2_itpc, tidx, fidx)];
    all_emg   = [get_pwr(ctrl_b1_pwr, tidx, fidx_emg); get_pwr(ctrl_b2_pwr, tidx, fidx_emg); get_pwr(casos_b1_pwr, tidx, fidx_emg); get_pwr(casos_b2_pwr, tidx, fidx_emg)];

    tbl = table(categorical(subj_ids), categorical(group_labels), categorical(block_labels), all_power, all_itpc, all_emg, ...
        'VariableNames', {'ID', 'Group', 'Block', 'Power', 'ITPC', 'EMG_Power'});

    lme_pwr  = fitlme(tbl, 'Power ~ Group * Block + EMG_Power + (1|ID)');
    lme_itpc = fitlme(tbl, 'ITPC ~ Group * Block + (1|ID)');

    fprintf('\n======================================================\n');
    fprintf('>> LME RESULTADOS: %s\n', upper(roi_name));
    fprintf('======================================================\n');
    fprintf(' [POWER] ERSP (Covariable EMG 70-90Hz incluida):\n');
    disp(dataset2table(anova(lme_pwr)));
    fprintf(' [ITPC] Phase Reset:\n');
    disp(dataset2table(anova(lme_itpc)));
end
fprintf('\n>>> PIPELINE ESTADÍSTICO FINALIZADO.\n');

%% ======================================================================
%  8. CORRELACIONES CEREBRO-CONDUCTA (BLOQUE 2 - CHEW)
% ======================================================================
fprintf('\n--- CORRELACIONES CEREBRO-CONDUCTA (BLOQUE 2: CHEW) ---\n');

load(f_beh, 'tb_data_45');
cas_ies_ch = double(tb_data_45.casos.chew.ies(:));
ctr_ies_ch = double(tb_data_45.controles.chew.ies(:));

for r = 1:size(tf_rois, 1)
    roi_name  = tf_rois{r,1};
    f_hz      = tf_rois{r,2};
    t_ms      = tf_rois{r,3};
    plot_name = tf_rois{r,4};

    fidx = dsearchn(frex', f_hz');
    tidx = dsearchn(times', t_ms');

    c2_pwr = get_pwr(ctrl_b2_pwr, tidx, fidx);
    s2_pwr = get_pwr(casos_b2_pwr, tidx, fidx);

    [r_cas, p_cas] = corr(s2_pwr, cas_ies_ch, 'Type', 'Pearson');
    [r_ctr, p_ctr] = corr(c2_pwr, ctr_ies_ch, 'Type', 'Pearson');

    cprintf(fid_stats, '  [%s - Bloque Chew] Casos:    r = % .3f, p = %.3f\n', plot_name, r_cas, p_cas);
    cprintf(fid_stats, '  [%s - Bloque Chew] Controles: r = % .3f, p = %.3f\n', plot_name, r_ctr, p_ctr);

    fig_corr = figure('Color','w','Position',[150 150 550 500], 'Visible', 'off');
    hold on;

    scatter(c2_pwr, ctr_ies_ch, 60, color_controls, 'filled', 'MarkerEdgeColor', 'w');
    scatter(s2_pwr, cas_ies_ch, 60, color_cases, 'filled', 'MarkerEdgeColor', 'w');

    p_fit_ctr = polyfit(c2_pwr, ctr_ies_ch, 1);
    x_ctr_line = linspace(min(c2_pwr), max(c2_pwr), 100);
    plot(x_ctr_line, polyval(p_fit_ctr, x_ctr_line), '-', 'Color', color_controls, 'LineWidth', 2.5);

    p_fit_cas = polyfit(s2_pwr, cas_ies_ch, 1);
    x_cas_line = linspace(min(s2_pwr), max(s2_pwr), 100);
    plot(x_cas_line, polyval(p_fit_cas, x_cas_line), '-', 'Color', color_cases, 'LineWidth', 2.5);

    set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, 'LineWidth', 1.2);
    xlabel('EEG Power (dB) [Chew]', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('IES (ms) [Chew]', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('Brain-Behavior (Chew Block): %s\n%d-%d ms', roi_name, t_ms(1), t_ms(2)), 'FontSize', 15);

    text_cas = sprintf('Cases: r=%.2f (p=%.3f)', r_cas, p_cas);
    text_ctr = sprintf('Controls: r=%.2f (p=%.3f)', r_ctr, p_ctr);
    legend({'Controls', 'Cases', text_ctr, text_cas}, 'Location', 'best', 'Box', 'off', 'FontSize', 11);

    fn_out_corr = sprintf('Corr_BrainBeh_B2_%s.png', plot_name);
    exportgraphics(fig_corr, fullfile(save_path, fn_out_corr), 'Resolution', 300);
    close(fig_corr);
end
fprintf('\n>>> CORRELACIONES B2 FINALIZADAS.\n');

% --- CORRELACIÓN ALFA vs CONDUCTA (Spearman, Delta) ----------------------
f_hz = P.bands_hz{2}; t_ms = [300 900];   % alpha desde config (8-12)
fidx = dsearchn(frex', f_hz');
tidx = dsearchn(times', t_ms');

pwr_Ch   = abs(get_pwr(casos_b2_pwr, tidx, fidx));
pwr_Nc   = abs(get_pwr(casos_b1_pwr, tidx, fidx));
Delta_dB = 10 * log10((pwr_Ch + eps) ./ (pwr_Nc + eps));
Delta_IES = tb_data_45.casos.chew.ies_m - tb_data_45.casos.nochew.ies_m;
mask = isfinite(Delta_dB) & isfinite(Delta_IES) & isreal(Delta_dB);

[rho, p] = corr(Delta_dB(mask), Delta_IES(mask), 'Type', 'Spearman', 'Rows', 'pairwise');
fprintf('\n>>> RESULTADO: ALFA (8-13 Hz) VS CONDUCTA <<<\n');
fprintf('Spearman rho: %.3f, p-value: %.4f (N=%d)\n', rho, p, sum(mask));

fig_alpha = figure('Color','w','Position',[200 200 500 450], 'Visible', 'off');
scatter(Delta_dB(mask), Delta_IES(mask), 80, color_cases, 'filled', 'MarkerFaceAlpha', 0.7);
lsline;
xlabel('\Delta Alpha Power (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('\Delta IES (ms)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('\\alpha vs Conducta: \\rho=%.3f, p=%.4f', rho, p), 'FontSize', 14);
set(gca, 'Box', 'off', 'TickDir', 'out');
grid on;
exportgraphics(fig_alpha, fullfile(save_path, 'Corr_Alpha_vs_DeltaIES.png'), 'Resolution', 300);
close(fig_alpha);

fprintf('\n✓ S2_TF_plot.m completado. Outputs en:\n  %s\n', save_path);

%% ======================================================================
%  FUNCIONES LOCALES
% ======================================================================

function tval = tstat2(d1, d2)
    [~,~,~,st] = ttest2(d1, d2, 'dim', 2); tval = st.tstat;
end

function plot_tf_map(data, mask, times, frex, lims, lbl, tit, time2, freq2, save_name)
    fig = figure('Color','w','Position',[100 100 560 460], 'Visible', 'off');
    ax  = axes('Color','w');
    contourf(ax, times, frex, data, 100, 'linecolor','none'); hold(ax,'on');
    if any(mask(:))
        contour(ax, times, frex, mask, 1, 'linecolor','k','LineWidth',2.5);
    end
    colormap(ax, jet); clim(ax, lims);
    cb = colorbar(ax); ylabel(cb, lbl, 'FontSize', 14);
    xlim(ax, [-200 1300]); ylim(ax, [1 30]);
    xlabel(ax, 'Time (ms)',      'FontSize', 14, 'FontWeight', 'bold');
    ylabel(ax, 'Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    title(ax,  tit,              'FontSize', 16, 'FontWeight', 'bold');
    exportgraphics(fig, save_name, 'Resolution', 300); close(fig);
end

function plot_topo_map(data, sigE, chanlocs, head_r, lims, lbl, tit, save_name)
    fig = figure('Color','w','Position',[100 100 480 460], 'Visible', 'off');
    % Llamada base — emarker2 como primer intento (puede fallar silenciosamente)
    if ~isempty(sigE)
        topoplotIndie(data, chanlocs, 'plotrad', 0.6, 'shading', 'flat', ...
                      'electrodes', 'on', 'emarker2', {sigE, 'o', 'k', 14, 3});
    else
        topoplotIndie(data, chanlocs, 'plotrad', 0.6, 'shading', 'flat', ...
                      'electrodes', 'on');
    end
    colormap(jet); clim(lims);
    % ---- Marcadores explícitos (garantiza visibilidad en todas las versiones) ----
    % EEGLAB coord: x = r*sin(theta_deg), y = r*cos(theta_deg), escalado por squeezefac
    if ~isempty(sigE)
        hold on;
        plotrad    = 0.6;
        squeezefac = head_r / plotrad;   % ≈ 0.833
        for ei = 1:numel(sigE)
            idx    = sigE(ei);
            th_rad = chanlocs(idx).theta * pi / 180;
            r_sc   = chanlocs(idx).radius * squeezefac;
            xp     = r_sc * sin(th_rad);
            yp     = r_sc * cos(th_rad);
            plot(xp, yp, 'ko', 'MarkerSize', 12, 'LineWidth', 3, ...
                 'MarkerFaceColor', 'k');
        end
        hold off;
    end
    cb = colorbar; ylabel(cb, lbl, 'FontSize', 14);
    title(tit, 'FontSize', 16, 'FontWeight', 'bold');
    exportgraphics(fig, save_name, 'Resolution', 300); close(fig);
end

function cprintf(fid, varargin)
    fprintf(1, varargin{:});
    if fid > 1, fprintf(fid, varargin{:}); end
end

function [sig_mask, thresh95] = perm_cluster(d1, d2, is_paired, n_perm, p_thresh)
    all_d = []; Ntot = 0; nC = 0;

    if is_paired
        df    = size(d1,3) - 1;
        t_thr = abs(tinv(p_thresh/2, df));
        [~,~,~,st] = ttest(d1, 0, 'dim', 3);
        t_map = st.tstat;
    else
        df    = size(d1,3) + size(d2,3) - 2;
        t_thr = abs(tinv(p_thresh/2, df));
        [~,~,~,st] = ttest2(d1, d2, 'dim', 3);
        t_map = st.tstat;
        all_d = cat(3, d1, d2);
        Ntot  = size(all_d,3);
        nC    = size(d1,3);
    end

    max_cls = zeros(n_perm,1);

    parfor pm = 1:n_perm
        if is_paired
            signs = (randi(2,1,size(d1,3))*2 - 3);
            [~,~,~,sp] = ttest(d1 .* reshape(signs,1,1,[]), 0, 'dim', 3);
        else
            idx = randperm(Ntot);
            [~,~,~,sp] = ttest2(all_d(:,:,idx(1:nC)), all_d(:,:,idx(nC+1:end)), 'dim', 3);
        end
        bm = abs(sp.tstat) > t_thr;
        cl = bwconncomp(bm, 8);
        if cl.NumObjects > 0
            max_cls(pm) = max(cellfun(@(n) sum(abs(sp.tstat(n))), cl.PixelIdxList));
        end
    end

    thresh95  = prctile(max_cls, 95);
    sig_mask  = zeros(size(t_map));
    bo = abs(t_map) > t_thr;
    co = bwconncomp(bo, 8);

    if co.NumObjects > 0
        su = cellfun(@(n) sum(abs(t_map(n))), co.PixelIdxList);
        for ci = 1:co.NumObjects
            if su(ci) > thresh95
                sig_mask(co.PixelIdxList{ci}) = 1;
            end
        end
    end
end