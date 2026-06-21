%% S2c_export_panels.m — Paneles multibanda para el paper
% Genera 12 figuras limpias (sin títulos ni stats):
%   9 topos : theta (4-7Hz) | alpha (8-13Hz) | theta+alpha (4-13Hz)
%             × {Cases Ch-Nc | Controls Ch-Nc | Interaction}
%   3 TF maps: Cases | Controls | Interaction
%
% Método banda: eegfilt (EEGLAB) + Hilbert → potencia media en WIN_LATE
% Los topos de theta se recomputan igual que alpha para consistencia.
% Destino: outputs/figures_ok/  |  300 dpi
% Correr desde Analysis_V1_Final/

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

OUT_OK = fullfile(ROOT_FINAL, 'outputs', 'figures_ok');
if ~exist(OUT_OK,'dir'), mkdir(OUT_OK); end

% Suprimir warning cosmético de polyshape (vértices duplicados en el borde del círculo)
warning('off', 'MATLAB:polyshape:repairedBySimplify');

%% ── Cargar TF stats (para los 3 TF maps) ────────────────────────────────
stats_file = fullfile(OUT_STATS, 'S2c_TF_GroupFigure.mat');
assert(exist(stats_file,'file')==2, 'Falta %s — correr S2c primero.', stats_file);
load(stats_file);   % map_*, mask_*, p_*, times_anal, FREX_TF

times_anal = times_anal(:)';
FREX_TF    = FREX_TF(:)';

%% ── Chanlocs ─────────────────────────────────────────────────────────────
ep_ref   = fullfile(DIR_EPOCHS, [CASES{1} '_Nc_ep.set']);
EEG_ref  = pop_loadset(ep_ref);
chanlocs = EEG_ref.chanlocs(1:EEG_N);
[xe, ye] = elec_xy(chanlocs);

%% ── Definición de bandas ─────────────────────────────────────────────────
bands = struct( ...
    'name',  {'theta',      'alpha',       'theta_alpha'}, ...
    'range', {[4 7],        [8 13],        [4 13]       }, ...
    'label', {'\theta 4-7Hz', '\alpha 8-13Hz', '\theta+\alpha 4-13Hz'} );
N_BANDS = numel(bands);

%% ── Computar potencia de banda por electrodo desde épocas ───────────────
% Matriz de salida: bp_cas_ch{b}(s, e) = potencia banda b, caso s, electrodo e, cond Ch
fprintf('\n[S2c_export_panels] Computando potencia de banda...\n');

bp_cas_ch = cell(N_BANDS,1);
bp_cas_nc = cell(N_BANDS,1);
bp_ctr_ch = cell(N_BANDS,1);
bp_ctr_nc = cell(N_BANDS,1);
for b = 1:N_BANDS
    bp_cas_ch{b} = zeros(N_CASES, EEG_N);
    bp_cas_nc{b} = zeros(N_CASES, EEG_N);
    bp_ctr_ch{b} = zeros(N_CONTROLS, EEG_N);
    bp_ctr_nc{b} = zeros(N_CONTROLS, EEG_N);
end

% Casos — cargar cada sujeto una sola vez (ambas condiciones y todas las bandas)
for s = 1:N_CASES
    subj = CASES{s};
    fprintf('  Casos %d/%d: %s\r', s, N_CASES, subj);
    [pw_ch, pw_nc] = epoch_multibandpow(DIR_EPOCHS, subj, bands, WIN_LATE, WIN_BASE, EEG_N);
    for b = 1:N_BANDS
        bp_cas_ch{b}(s,:) = pw_ch(b,:);
        bp_cas_nc{b}(s,:) = pw_nc(b,:);
    end
end
fprintf('\n');

% Controles
for s = 1:N_CONTROLS
    subj = CONTROLS{s};
    fprintf('  Controles %d/%d: %s\r', s, N_CONTROLS, subj);
    [pw_ch, pw_nc] = epoch_multibandpow(DIR_EPOCHS, subj, bands, WIN_LATE, WIN_BASE, EEG_N);
    for b = 1:N_BANDS
        bp_ctr_ch{b}(s,:) = pw_ch(b,:);
        bp_ctr_nc{b}(s,:) = pw_nc(b,:);
    end
end
fprintf('\n');

%% ── Contrastes y FDR por banda ───────────────────────────────────────────
topo  = struct();  % topo.(cas/ctr/int){b} = [1×64]
sigfdr = struct(); % sigfdr.(cas/ctr/int){b} = logical [64×1]

for b = 1:N_BANDS
    % bp_*{b} ya son dB normalizados a baseline (10·log10(WIN_LATE/WIN_BASE))
    % El contraste es la diferencia de esos dB entre condiciones
    diff_cas = bp_cas_ch{b} - bp_cas_nc{b};  % [31×64]
    diff_ctr = bp_ctr_ch{b} - bp_ctr_nc{b};  % [15×64]

    topo.cas{b} = median(diff_cas, 1);
    topo.ctr{b} = median(diff_ctr, 1);
    topo.int{b} = topo.cas{b} - topo.ctr{b};

    % FDR Benjamini-Hochberg por electrodo
    % NOTA: variables con prefijo 'pe_' para no sobreescribir p_cas/p_ctr/p_int
    % escalares del cluster TF cargados desde S2c_TF_GroupFigure.mat
    [~, pe_cas] = ttest(diff_cas);
    [~, pe_ctr] = ttest(diff_ctr);
    [~, pe_int] = ttest2(diff_cas, diff_ctr);

    sigfdr.cas{b} = fdr_bh(pe_cas(:)) < 0.05;
    sigfdr.ctr{b} = fdr_bh(pe_ctr(:)) < 0.05;
    sigfdr.int{b} = fdr_bh(pe_int(:)) < 0.05;

    fprintf('  %s — sig FDR: Cases=%d, Controls=%d, Interaction=%d\n', ...
        bands(b).name, sum(sigfdr.cas{b}), sum(sigfdr.ctr{b}), sum(sigfdr.int{b}));
end

%% ── Colorescalas ─────────────────────────────────────────────────────────
clim_tf = 1.5;

% Una escala por banda (compartida entre los 3 contrastes de esa banda)
clim_band = zeros(N_BANDS, 1);
for b = 1:N_BANDS
    v = abs([topo.cas{b}(:); topo.ctr{b}(:); topo.int{b}(:)]);
    clim_band(b) = max(v) * 1.1;
    if clim_band(b) == 0, clim_band(b) = 1; end
end

%% ── Tabla de paneles ──────────────────────────────────────────────────────
% {tipo, datos, mascara, pval, sig_logical, band_label, clim, fname}
panels = {};

for b = 1:N_BANDS
    bl  = bands(b).label;
    bn  = bands(b).name;
    cl  = clim_band(b);

    panels(end+1,:) = {'topo', topo.cas{b}(:), [], [], sigfdr.cas{b}, bl, cl, ...
        sprintf('topo_%s_cases_chew_vs_nochew', bn)};
    panels(end+1,:) = {'topo', topo.ctr{b}(:), [], [], sigfdr.ctr{b}, bl, cl, ...
        sprintf('topo_%s_controls_chew_vs_nochew', bn)};
    panels(end+1,:) = {'topo', topo.int{b}(:), [], [], sigfdr.int{b}, bl, cl, ...
        sprintf('topo_%s_interaction_cases_minus_controls', bn)};
end

panels(end+1,:) = {'tf', map_cas, mask_cas, p_cas, [], 'Cases Ch-Nc',    clim_tf, 'tf_cases_chew_vs_nochew'};
panels(end+1,:) = {'tf', map_ctr, mask_ctr, p_ctr, [], 'Controls Ch-Nc', clim_tf, 'tf_controls_chew_vs_nochew'};
panels(end+1,:) = {'tf', map_int, mask_int, p_int, [], 'Interaction',     clim_tf, 'tf_interaction_cases_minus_controls'};

%% ── Generar paneles ───────────────────────────────────────────────────────
n_panels = size(panels, 1);
fprintf('\nGenerando %d paneles...\n', n_panels);

for k = 1:n_panels
    tipo      = panels{k,1};
    datos     = panels{k,2};
    mascara   = panels{k,3};
    pval      = panels{k,4};
    sig_log   = panels{k,5};
    blabel    = panels{k,6};
    clim_val  = panels{k,7};
    fname     = panels{k,8};

    if strcmp(tipo, 'topo')

        fig = figure('Visible','off','Color','w', ...
                     'Units','centimeters','Position',[0 0 12 8]);

        sig_idx = find(sig_log);
        try
            if ~isempty(sig_idx)
                topoplot(datos, chanlocs, 'electrodes','off', 'whitebk','on', ...
                         'emarker2', {sig_idx, '.', 'k', 22, 2});
            else
                topoplot(datos, chanlocs, 'electrodes','off', 'whitebk','on');
            end
        catch
            topoplotIndie(datos, chanlocs, 'electrodes','off');
            if ~isempty(sig_idx)
                hold on;
                plot(xe(sig_idx), ye(sig_idx), 'k.', 'MarkerSize', 22);
            end
        end

        ax = gca;
        set(ax, 'CLim', [-clim_val clim_val]);
        colormap(ax, jet);

        % Fondo blanco: máscara alpha en la superficie + polyshape sobre esquinas
        topo_background_fix(ax, 0.511);
        try
            xl = xlim(ax); yl = ylim(ax);
            ang = linspace(0, 2*pi, 361);
            ps_o = polyshape([xl(1) xl(2) xl(2) xl(1)], [yl(1) yl(1) yl(2) yl(2)]);
            ps_i = polyshape(cos(ang)*0.511, sin(ang)*0.511);
            hold on;
            ph = plot(subtract(ps_o, ps_i));
            ph.FaceColor = 'w'; ph.EdgeColor = 'none'; ph.FaceAlpha = 1;
        catch; end
        set(ax,'Color','w'); set(fig,'Color','w');

        % Etiqueta de banda (esquina inferior izquierda del topoplot)
        hold on;
        text(-0.45, -0.47, blabel, 'FontSize', FIG_FS+1, 'FontWeight','bold', ...
             'Color','k', 'FontName', FIG_FONT, 'HorizontalAlignment','left');

        % Colorbar delgado
        drawnow;
        pos0 = ax.Position;
        new_w = pos0(3) * 0.76;
        ax.Position = [pos0(1), pos0(2), new_w, pos0(4)];
        cb = colorbar(ax, 'eastoutside');
        cb.Position = [pos0(1)+new_w+0.04, pos0(2)+pos0(4)*0.10, 0.04, pos0(4)*0.80];
        ax.Position = [pos0(1), pos0(2), new_w, pos0(4)];
        set(ax, 'FontName',FIG_FONT, 'FontSize',FIG_FS);

    else  % TF map

        fig = figure('Visible','off','Color','w', ...
                     'Units','centimeters','Position',[0 0 12 7]);

        imagesc(times_anal, FREX_TF, datos);
        axis xy;
        colormap(gca, jet);
        set(gca, 'CLim', [-clim_val clim_val]);
        colorbar;
        hold on;

        if ~isempty(mascara) && any(mascara(:)) && pval < ALPHA_CLUST
            contour(times_anal, FREX_TF, double(mascara), 1, 'w-', 'LineWidth', 2.5);
        end

        yline(BAND_THETA(2), 'w:', 'LineWidth', 1.2);
        yline(BAND_ALPHA(2), 'w:', 'LineWidth', 1.2);
        xline(0,             'k--', 'LineWidth', 1.5);
        xline(WIN_EARLY(1),  'w:',  'LineWidth', 1.0);
        xline(WIN_LATE(1),   'w:',  'LineWidth', 1.0);
        xline(WIN_LATE(2),   'w:',  'LineWidth', 1.0);

        xlabel('Time (ms)', 'FontSize', FIG_FS);
        ylabel('Frequency (Hz)', 'FontSize', FIG_FS);
        xlim([times_anal(1) times_anal(end)]);
        ylim([FREX_TF(1) FREX_TF(end)]);
        set(gca, 'FontName',FIG_FONT, 'FontSize',FIG_FS);
    end

    out_path = fullfile(OUT_OK, fname);
    print(fig, out_path, '-dpng', '-r300');
    close(fig);
    fprintf('  OK (%d/%d): %s\n', k, n_panels, fname);
end

fprintf('\n[S2c_export_panels] %d paneles en figures_ok/\n', n_panels);

%% ========== Funciones locales ==========

function [pw_ch, pw_nc] = epoch_multibandpow(dir_ep, subj, bands, win_ms, win_base_ms, n_ch)
% Carga épocas Ch y Nc y devuelve dB de banda normalizados a baseline.
% pw_ch(b, e) = 10·log10(pow_WIN_LATE / pow_WIN_BASE) para banda b, electrodo e.
% Usa eegfilt + Hilbert envelope, promedio sobre trials.

    n_bands = numel(bands);
    pw_ch   = NaN(n_bands, n_ch);
    pw_nc   = NaN(n_bands, n_ch);

    fname_ch = fullfile(dir_ep, [subj '_Ch_ep.set']);
    fname_nc = fullfile(dir_ep, [subj '_Nc_ep.set']);

    if ~exist(fname_ch,'file') || ~exist(fname_nc,'file')
        warning('Archivos no encontrados para %s', subj);
        return;
    end

    EEG_ch = pop_loadset(fname_ch);
    EEG_nc = pop_loadset(fname_nc);
    fs     = EEG_ch.srate;

    tidx_ch   = EEG_ch.times >= win_ms(1)      & EEG_ch.times <= win_ms(2);
    tidx_nc   = EEG_nc.times >= win_ms(1)      & EEG_nc.times <= win_ms(2);
    bidx_ch   = EEG_ch.times >= win_base_ms(1) & EEG_ch.times <= win_base_ms(2);
    bidx_nc   = EEG_nc.times >= win_base_ms(1) & EEG_nc.times <= win_base_ms(2);

    for b = 1:n_bands
        flo = bands(b).range(1);
        fhi = bands(b).range(2);

        pw_ch(b,:) = bandpow_db(EEG_ch, flo, fhi, fs, n_ch, tidx_ch, bidx_ch);
        pw_nc(b,:) = bandpow_db(EEG_nc, flo, fhi, fs, n_ch, tidx_nc, bidx_nc);
    end
end

function pow_db = bandpow_db(EEG, flo, fhi, fs, n_ch, tidx_win, tidx_base)
% Potencia de banda normalizada a baseline por electrodo.
% Resultado: 10·log10(mean_pow_WIN / mean_pow_BASE), media sobre trials.
% Devuelve [1 × n_ch].

    n_trials = EEG.trials;
    n_times  = EEG.pnts;

    % Filtrar todos los trials concatenados (más estable que filtrar por separado)
    data_concat = reshape(EEG.data(1:n_ch, :, :), n_ch, n_trials * n_times);

    try
        data_filt = eegfilt(data_concat, fs, flo, fhi);
    catch
        [b_c, a_c] = butter(4, [flo fhi] / (fs/2), 'bandpass');
        data_filt  = filtfilt(b_c, a_c, data_concat')';
    end

    % Envelope de potencia instantánea [n_ch × n_times × n_trials]
    env   = abs(hilbert(data_filt'))';
    env3d = reshape(env, n_ch, n_times, n_trials);
    pow2  = env3d .^ 2;

    % Media sobre la ventana de análisis y sobre trials
    p_win  = squeeze(mean(mean(pow2(:, tidx_win,  :), 2), 3));  % [n_ch × 1]
    p_base = squeeze(mean(mean(pow2(:, tidx_base, :), 2), 3));  % [n_ch × 1]

    % dB normalizados: 10·log10(WIN / BASE)
    pow_db = (10 * log10((p_win + eps) ./ (p_base + eps)))';    % [1 × n_ch]
end

function topo_background_fix(ax, head_r)
% Aplica alpha=0 fuera del círculo de la cabeza en la Surface de topoplot.
    for h = ax.Children'
        try
            if strcmp(class(h), 'matlab.graphics.primitive.Surface') && numel(h.XData) > 16
                dist = sqrt(h.XData.^2 + h.YData.^2);
                h.AlphaData = double(dist <= head_r);
                h.FaceAlpha = 'flat';
                break;
            end
        catch; end
    end
end

function p_adj = fdr_bh(p_vals)
% FDR Benjamini-Hochberg. Input/output: vector columna.
    p_vals = p_vals(:);
    n = numel(p_vals);
    [ps, si] = sort(p_vals);
    pa = ps .* n ./ (1:n)';
    for i = n-1:-1:1, pa(i) = min(pa(i), pa(i+1)); end
    p_adj = zeros(n,1);
    p_adj(si) = pa;
end

function [xp, yp] = elec_xy(clocs)
% Coordenadas de electrodos en espacio topoplot EEGLAB.
    headrad = 0.511;
    n = numel(clocs);
    xp = zeros(1,n); yp = zeros(1,n);
    for i = 1:n
        th = clocs(i).theta * pi/180;
        rd = clocs(i).radius;
        xp(i) = sin(th) * rd / headrad * 0.5;
        yp(i) = cos(th) * rd / headrad * 0.5;
    end
end
