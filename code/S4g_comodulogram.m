%% S4g_comodulogram.m
% Comodulogram PAC: grilla EMG-fase × EEG-amplitud
% Condición: Cases chewing (N=31)
% MI de Tort (18 bins) para cada par de frecuencias.
%
% DECISIONES DE ANÁLISIS:
%   ROI    : 18 electrodos frontales, theta-interacción FDR (S2c) — igual que S4b
%   Datos  : epochs concatenados de _Ch_ep.set (-1000:1500 ms)
%             → solo ventana tarea activa; EMG continuo se reconstruye
%   Grilla : Fase 0.4–2.4 Hz (paso 0.2, BW ±0.15, Butter 2)
%            Amp  4–30 Hz   (paso 1,   BW ±2,    Butter 4)
% Colormap: parula
%
% Salidas:
%   stats/  S4g_comodulogram.mat
%   figures/S4g_comodulogram.png  (MI raw)
%   figures/S4g_comodulogram_zscore.png

clc;   % (sin 'clear': permite pipelining desde run_S4gh.m)

%% ── Config ───────────────────────────────────────────────────────────────
THIS_DIR  = fileparts(mfilename('fullpath'));
ROOT_V1F  = fileparts(THIS_DIR);
ROOT_P2V1 = fileparts(ROOT_V1F);

cfg_path = fullfile(ROOT_V1F, 'S0_config.m');
run(cfg_path);

% Fallbacks si S0_config no resolvió bien las rutas
if ~exist('DATA_PAC','var') || ~exist(DATA_PAC,'dir')
    DATA_PAC = fullfile(ROOT_P2V1, 'Data_PAC');
end
if ~exist('OUT_STATS','var') || ~exist(OUT_STATS,'dir')
    OUT_STATS = fullfile(ROOT_V1F, 'outputs', 'stats');
    if ~exist(OUT_STATS,'dir'), mkdir(OUT_STATS); end
end
if ~exist('OUT_FIGS','var') || ~exist(OUT_FIGS,'dir')
    OUT_FIGS = fullfile(ROOT_V1F, 'outputs', 'figures');
    if ~exist(OUT_FIGS,'dir'), mkdir(OUT_FIGS); end
end
if ~exist('EMG_CHANS','var'), EMG_CHANS = [65 66]; end
if ~exist('MI_BINS','var'),   MI_BINS   = 18;       end

% Suffix: continuo del bloque chewing (con EMG limpio)
CONT_SUFFIX = '_Ch_clean_emg.set';

% Diagnóstico de rutas
fprintf('DATA_PAC : %s  (existe: %d)\n', DATA_PAC, exist(DATA_PAC,'dir'));
fname_test = fullfile(DATA_PAC, [CASES{1} CONT_SUFFIX]);
fprintf('Test file: %s\n  (existe: %d)\n', fname_test, exist(fname_test,'file'));

%% ── ROI EEG: 18 electrodos frontales theta-interacción (S0_config → S2c FDR) ──
% ROI_CBPT viene de S0_config.m — igual que S4b_PAC_ROI.m
% NO cargar roi_cbpt.mat (ese tiene 3 electrodos FWER de S2b, no se usa aquí)
ROI_LABELS = ROI_CBPT;
fprintf('ROI CBPT (%d electrodos): %s\n', numel(ROI_LABELS), strjoin(ROI_LABELS, ' '));

%% ── Grilla de frecuencias ────────────────────────────────────────────────
F_PHASE = 0.4 : 0.2 : 2.4;    % 11 pasos — rango masticatorio
F_AMP   = 4   : 1  : 30;      % 27 pasos — theta → beta
BW_P    = 0.15;                % ±Hz filtro fase
BW_A    = 2.0;                 % ±Hz filtro amplitud
N_P     = numel(F_PHASE);
N_A     = numel(F_AMP);
N_BINS  = MI_BINS;
TRIM    = 512;   % 2 s a 256 Hz — descarta bordes post-filtrado

%% ── Verificar EEGLAB ─────────────────────────────────────────────────────
if isempty(which('pop_loadset'))
    error('pop_loadset no encontrado. Ejecuta ''eeglab'' primero.');
end
fprintf('EEGLAB: %s\n\n', fileparts(which('pop_loadset')));

%% ── Pre-allocate ─────────────────────────────────────────────────────────
MI_all = nan(N_P, N_A, N_CASES);

fprintf('Comodulogram: %d sujetos × %d fase × %d amp\n', N_CASES, N_P, N_A);
fprintf('Grilla fase: %.1f–%.1f Hz | BW ±%.2f Hz | Butter 2\n', ...
    F_PHASE(1), F_PHASE(end), BW_P);
fprintf('Grilla amp:  %d–%d Hz   | BW ±%.1f Hz | Butter 4\n\n', ...
    F_AMP(1), F_AMP(end), BW_A);

%% ── Loop principal ───────────────────────────────────────────────────────
for s = 1:N_CASES
    fname = fullfile(DATA_PAC, [CASES{s} CONT_SUFFIX]);
    if ~exist(fname, 'file')
        fprintf('  [%s] archivo no encontrado — omitido\n', CASES{s});
        continue
    end

    EEG = pop_loadset(fname);
    fs  = EEG.srate;
    nyq = fs / 2;

    % ── EMG: criterio único bilateral + fallback (R6, igual que S4b/S4h) ──
    ch_emg = EMG_CHANS(EMG_CHANS <= EEG.nbchan);
    if isempty(ch_emg), continue; end
    emg_raw = emg_bilateral(EEG.data, ch_emg, fs);   % [1 × T]

    % ── ROI EEG: 18 electrodos theta-interacción (ROI_CBPT, S0_config) ───
    n_eeg = min(64, EEG.nbchan);
    [~, roi_idx] = ismember(ROI_LABELS, {EEG.chanlocs(1:n_eeg).labels});
    roi_idx(roi_idx == 0) = [];
    if isempty(roi_idx)
        fprintf('  [%s] ROI no encontrado en chanlocs — omitido\n', CASES{s});
        continue
    end
    eeg_roi = mean(double(EEG.data(roi_idx, :)), 1);   % [1 × T]

    T_sig  = numel(emg_raw);
    trim_s = min(TRIM, floor(T_sig * 0.05));

    % ── MI para cada par (fp, fa) ─────────────────────────────────────────
    MI_s = nan(N_P, N_A);

    for ip = 1:N_P
        fp   = F_PHASE(ip);
        lo_p = max(0.05, fp - BW_P);
        hi_p = min(nyq  - 0.5, fp + BW_P);

        [b_p, a_p] = butter(2, [lo_p hi_p] / nyq, 'bandpass');
        emg_filt   = filtfilt(b_p, a_p, emg_raw);
        phase_sig  = angle(hilbert(emg_filt));

        for ia = 1:N_A
            fa   = F_AMP(ia);
            lo_a = max(1.0,       fa - BW_A);
            hi_a = min(nyq - 0.5, fa + BW_A);

            [b_a, a_a] = butter(4, [lo_a hi_a] / nyq, 'bandpass');
            eeg_filt   = filtfilt(b_a, a_a, eeg_roi);
            amp_sig    = abs(hilbert(eeg_filt));

            t1 = trim_s + 1;
            t2 = T_sig  - trim_s;
            if t2 - t1 < fs * 2   % mínimo 2 s útiles
                continue
            end
            MI_s(ip, ia) = tort_mi(phase_sig(t1:t2)', amp_sig(t1:t2)', N_BINS);
        end
    end

    MI_all(:, :, s) = MI_s;
    n_ok = sum(~isnan(MI_s(:)));
    fprintf('  [%d/%d] %s (fs=%d, T=%.0fs) — %d/%d pares OK\n', ...
        s, N_CASES, CASES{s}, fs, T_sig/fs, n_ok, N_P*N_A);
end

%% ── Estadísticas grupo ───────────────────────────────────────────────────
MI_mean = mean(MI_all, 3, 'omitnan');
MI_sem  = std(MI_all,  0, 3, 'omitnan') / sqrt(N_CASES);
MI_z    = (MI_mean - mean(MI_mean(:))) / std(MI_mean(:));

n_valid = sum(~isnan(MI_all(1,1,:)));
fprintf('\nSujetos con datos válidos: %d/%d\n', n_valid, N_CASES);

%% ── Guardar ──────────────────────────────────────────────────────────────
save(fullfile(OUT_STATS, 'S4g_comodulogram.mat'), ...
    'MI_all','MI_mean','MI_sem','MI_z','F_PHASE','F_AMP', ...
    'BW_P','BW_A','N_BINS','ROI_LABELS','N_CASES');
fprintf('Guardado: S4g_comodulogram.mat\n');

%% ── Figura 1: MI raw ─────────────────────────────────────────────────────
valid_mi = MI_mean(~isnan(MI_mean(:)));
clim_raw = [];
if numel(valid_mi) > 1 && range(valid_mi) > 0
    clim_raw = [min(valid_mi), prctile(valid_mi, 98)];
end

fh = figure('Units','inches','Position',[1 1 5.5 4.5],'Color','w');
ax = axes(fh);
draw_como(ax, F_AMP, F_PHASE, MI_mean, 'MI (Tort)', clim_raw);
title(ax, sprintf('PAC Comodulogram — Cases (N=%d, Chewing)', n_valid), ...
    'FontSize', 12, 'FontWeight', 'bold');
set(fh, 'PaperPositionMode', 'auto');
print(fh, fullfile(OUT_FIGS, 'S4g_comodulogram'), '-dpng', '-r300');
close(fh);
fprintf('Guardada: S4g_comodulogram.png\n');

%% ── Figura 2: z-score ────────────────────────────────────────────────────
fh2 = figure('Units','inches','Position',[1 1 5.5 4.5],'Color','w');
ax2 = axes(fh2);
draw_como(ax2, F_AMP, F_PHASE, MI_z, 'MI (z-score)', []);
title(ax2, sprintf('PAC Comodulogram — z-score (N=%d)', n_valid), ...
    'FontSize', 12, 'FontWeight', 'bold');
set(fh2, 'PaperPositionMode', 'auto');
print(fh2, fullfile(OUT_FIGS, 'S4g_comodulogram_zscore'), '-dpng', '-r300');
close(fh2);
fprintf('Guardada: S4g_comodulogram_zscore.png\n');
fprintf('\n=== FIN comodulogram ===\n');

%% ════════════════════════════════════════════════════════════════════════
%  FUNCIONES LOCALES  (deben estar al final del script)
%% ════════════════════════════════════════════════════════════════════════

function draw_como(ax, F_AMP, F_PHASE, MI_map, cbar_lbl, clim_in)
% Dibuja un comodulogram en el axes dado con etiquetas de banda internas.
    imagesc(ax, F_AMP, F_PHASE, MI_map);
    axis(ax, 'xy');
    colormap(ax, parula);
    cb = colorbar(ax);
    cb.Label.String   = cbar_lbl;
    cb.Label.FontSize = 10;
    hold(ax, 'on');

    band_edges = [4 7; 8 13; 13 30];
    band_names = {'\theta', '\alpha', '\beta'};
    lc = [0.9 0.9 0.9];   % color líneas de referencia

    for b = 1:3
        for e = 1:2
            xv = band_edges(b, e);
            plot(ax, [xv xv], [F_PHASE(1) F_PHASE(end)], ...
                 '--', 'Color', [lc 0.55], 'LineWidth', 0.9);
        end
        xmid = mean(band_edges(b, :));
        text(ax, xmid, F_PHASE(end) - 0.10, band_names{b}, ...
             'HorizontalAlignment', 'center', 'FontSize', 10, ...
             'Color', 'w', 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end

    % Rango fase masticatoria original (0.5–2.0 Hz)
    for yr = [0.5, 2.0]
        if yr >= F_PHASE(1) && yr <= F_PHASE(end)
            plot(ax, [F_AMP(1) F_AMP(end)], [yr yr], ...
                 '--', 'Color', [lc 0.55], 'LineWidth', 0.9);
        end
    end

    % CLim
    if ~isempty(clim_in) && numel(clim_in) == 2 && diff(clim_in) > 0
        set(ax, 'CLim', clim_in);
    else
        valid = MI_map(~isnan(MI_map(:)));
        if numel(valid) > 1 && range(valid) > 0
            set(ax, 'CLim', [min(valid), prctile(valid, 98)]);
        end
    end

    xlabel(ax, 'Amplitude frequency (Hz)', 'FontSize', 11);
    ylabel(ax, 'Phase frequency (Hz)',      'FontSize', 11);
    ax.XTick      = [4 7 8 13 15 20 25 30];
    ax.YTick      = F_PHASE;
    ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), F_PHASE, 'UniformOutput', false);
    ax.TickDir    = 'out';
    ax.FontSize   = 9;
    ax.Box        = 'off';
    ax.XLim       = [F_AMP(1) - 0.3,  F_AMP(end) + 0.3];
    ax.YLim       = [F_PHASE(1) - 0.05, F_PHASE(end) + 0.05];
end
