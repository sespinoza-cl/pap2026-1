%% S2c_replot.m — Regenera la figura S2c sin rehacer el CBPT
% Carga S2c_TF_GroupFigure.mat y redibuja con electrodos sig marcados.
% Correr desde Analysis_V1_Final/

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

%% Cargar resultados del CBPT
stats_file = fullfile(OUT_STATS, 'S2c_TF_GroupFigure.mat');
assert(exist(stats_file,'file')==2, 'Falta %s — correr S2c primero', stats_file);
load(stats_file);

%% Cargar theta_topo y chanlocs
load(FILE_TF, 'theta_topo_ch','theta_topo_nc', ...
              'theta_topo_ch_ctr','theta_topo_nc_ctr');

ep_ref = fullfile(DIR_EPOCHS, [CASES{1} '_Nc_ep.set']);
if exist(ep_ref,'file')
    EEG_ref  = pop_loadset(ep_ref);
    chanlocs = EEG_ref.chanlocs(1:EEG_N);
else
    error('No se encontró %s para cargar chanlocs', ep_ref);
end

%% Recompute topografías (por si no estaban en el .mat)
if ~exist('topo_cas','var')
    topo_cas = median(theta_topo_ch    - theta_topo_nc,    1);
    topo_ctr = median(theta_topo_ch_ctr - theta_topo_nc_ctr, 1);
    topo_int = topo_cas - topo_ctr;
end

%% Colorescalas
clim_tf   = 1.5;
clim_topo = max(abs([topo_cas topo_ctr topo_int])) * 1.1;
if clim_topo == 0, clim_topo = 1; end

times_anal = times_anal(:)';
FREX_TF    = FREX_TF(:)';

%% Datos para el loop de figura
titles_topo = {sprintf('Cases  Ch > Nc  (N=%d)', N_CASES), ...
               sprintf('Controls  Ch > Nc  (N=%d)', N_CONTROLS), ...
               'Interaction  (Cases \minus Controls)'};
titles_tf   = {sprintf('Cases  p_{cluster}=%.4f', p_cas), ...
               sprintf('Controls  p_{cluster}=%.4f', p_ctr), ...
               sprintf('Interaction  p_{cluster}=%.4f', p_int)};

topo_data = {topo_cas(:), topo_ctr(:), topo_int(:)};
topo_sig  = {sig_cas_fdr(:), sig_ctr_fdr(:), sig_int_fdr(:)};
tf_maps   = {map_cas, map_ctr, map_int};
tf_masks  = {mask_cas, mask_ctr, mask_int};
tf_pvals  = {p_cas, p_ctr, p_int};

%% Figura
fig = figure('Name','S2c TF+Topo (replot)','Units','normalized',...
             'Position',[0.01 0.02 0.97 0.90],'Color','w');

for col = 1:3
    % ── Fila 1: Topo ──────────────────────────────────────────────────────
    subplot(2,3,col);
    topo_v  = topo_data{col};
    sig_idx = find(topo_sig{col});

    try
        if ~isempty(sig_idx)
            topoplot(topo_v, chanlocs, 'electrodes','off', ...
                     'emarker2', {sig_idx, '.', 'k', 18, 2});
        else
            topoplot(topo_v, chanlocs, 'electrodes','off');
        end
    catch
        topoplotIndie(topo_v, chanlocs, 'electrodes','off');
    end

    set(gca,'CLim',[-clim_topo clim_topo]);
    colormap(gca, jet);
    cb = colorbar; cb.Label.String = '\Delta\theta power (dB)';

    if ~isempty(sig_idx)
        sig_labels = strjoin({chanlocs(sig_idx).labels}, ', ');
        title({titles_topo{col}, sprintf('\\bullet sig: %s', sig_labels)}, ...
            'FontSize', FIG_FS, 'FontWeight','bold');
    else
        title({titles_topo{col}, 'no sig. electrodes'}, ...
            'FontSize', FIG_FS, 'FontWeight','bold');
    end
    axis square;
    set(gca,'FontName',FIG_FONT,'FontSize',FIG_FS);

    % ── Fila 2: TF map ───────────────────────────────────────────────────
    subplot(2,3,col+3);
    imagesc(times_anal, FREX_TF, tf_maps{col});
    axis xy;
    colormap(gca, jet);
    set(gca,'CLim',[-clim_tf clim_tf]);
    cb2 = colorbar; cb2.Label.String = '\Delta power (dB)';
    hold on;

    % Contorno del cluster significativo
    mask = tf_masks{col};
    if ~isempty(mask) && any(mask(:)) && tf_pvals{col} < ALPHA_CLUST
        contour(times_anal, FREX_TF, double(mask), 1, 'w-', 'LineWidth', 2.5);
    end

    % Bandas de frecuencia (líneas horizontales)
    yline(BAND_THETA(2), 'w:', 'LineWidth', 1.2);   % 7 Hz
    yline(BAND_ALPHA(2), 'w:', 'LineWidth', 1.2);   % 13 Hz

    % Ventanas temporales (líneas verticales)
    xline(0,           'k--', 'LineWidth', 1.5);
    xline(WIN_EARLY(1),'w:',  'LineWidth', 1);
    xline(WIN_LATE(1), 'w:',  'LineWidth', 1);
    xline(WIN_LATE(2), 'w:',  'LineWidth', 1);

    xlabel('Time (ms)','FontSize',FIG_FS);
    ylabel('Frequency (Hz)','FontSize',FIG_FS);
    title(titles_tf{col},'FontSize',FIG_FS+1,'FontWeight','bold');
    xlim([times_anal(1) times_anal(end)]);
    ylim([FREX_TF(1) FREX_TF(end)]);
    set(gca,'FontName',FIG_FONT,'FontSize',FIG_FS);
end

sgtitle(sprintf('TF 2D CBPT | \\alpha_{voxel}=%.2f | N_{PERM}=%d | [%.0f, %.0f] ms', ...
    ALPHA_THRESH, N_PERM, WIN_ANALYSIS(1), WIN_ANALYSIS(2)), ...
    'FontSize', FIG_FS+2, 'FontWeight','bold');

saveas(fig, fullfile(OUT_FIGS,'S2c_TF_GroupFigure.png'));
print(fig, fullfile(OUT_FIGS,'S2c_TF_GroupFigure'), '-dpng', '-r300');
fprintf('✓ Figura guardada: %s\n', fullfile(OUT_FIGS,'S2c_TF_GroupFigure.png'));
