%% S4h_comodulogram_stats.m
% Comodulogram estadístico: z-map por sujeto + CBPT 2D
%
% LÓGICA:
%   1. Por sujeto (parfor): MI observado en grilla N_P × N_A
%   2. Por sujeto (parfor): N_SURR surrogates por circular-shift de amplitud
%      → z_map(ip,ia) = (MI_obs - mean_null) / std_null
%   3. Grupo: t-test one-sample (z > 0) celda a celda → t_map [N_P × N_A]
%   4. CBPT 2D: sign-flip permutation sobre z_maps, estadístico max-cluster-sum
%   5. Figura: MI_mean + contorno cluster significativo
%
% DECISIONES DE ANÁLISIS:
%   ROI EEG  : ROI_CBPT (18 electrodos frontales, S0_config) — igual que S4b
%   Fase EMG : sweep 0.4–2.4 Hz paso 0.2, BW ±0.15 Hz, Butter 2
%              (grilla exploratoria — la f_chew individual se marca como ref.)
%              Para el PAC principal (S4b) sí se usa f_chew_ind ± 0.5 Hz.
%   Amp EEG  : sweep 4–30 Hz paso 1, BW ±2 Hz, Butter 4
%   Null     : N_SURR circular-shifts de amplitud EEG; el mismo shift
%              se aplica a todas las f_amp de un surrogate dado (preserva
%              covarianza entre bandas y evita N_SURR × N_A × N_P filtrados)
%   CBPT     : sign-flip sobre z_maps centrados, max-cluster-sum, N_PERM=5000
%
% Salidas:
%   stats/  S4h_comodulogram_stats.mat
%   figures/S4h_MI_cluster.png    (MI_mean + contorno cluster sig.)
%   figures/S4h_zmap_cluster.png  (z-map grupal + contorno cluster sig.)

clc;   % (sin 'clear': permite pipelining desde run_S4gh.m)

%% ── Config ───────────────────────────────────────────────────────────────────
THIS_DIR = fileparts(mfilename('fullpath'));
ROOT_V1F = fileparts(THIS_DIR);
run(fullfile(ROOT_V1F, 'S0_config.m'));

CONT_SUFFIX = '_Ch_clean_emg.set';

%% ── Parámetros grilla ────────────────────────────────────────────────────────
F_PHASE = 0.4 : 0.2 : 2.4;    % [Hz] 11 puntos — rango masticatorio
F_AMP   = 4   : 1  : 30;      % [Hz] 27 puntos — theta → beta
BW_P    = 0.15;                % [Hz] semi-ancho filtro fase
BW_A    = 2.0;                 % [Hz] semi-ancho filtro amplitud
N_P     = numel(F_PHASE);      % 11
N_A     = numel(F_AMP);        % 27
N_BINS  = MI_BINS;             % 18 (desde S0_config)
TRIM    = 512;                 % muestras descartadas en bordes post-filtrado

%% ── Parámetros estadísticos ──────────────────────────────────────────────────
% N_SURR heredado de S0_config (=500) — unificado con S4b; antes hardcoded 200 (I5)
N_PERM      = 5000;   % iteraciones CBPT sign-flip
ALPHA_VOXEL = 0.05;   % umbral voxel para definir clusters (p uncorrected)
MIN_SHIFT_S = 12.5;   % [s] 5 ciclos a 0.4 Hz (f_fase más lenta del sweep) — R7
rng(RNG_SEED);        % reproducibilidad

%% ── EEGLAB ───────────────────────────────────────────────────────────────────
if isempty(which('pop_loadset'))
    error('Ejecuta "eeglab" primero y vuelve a correr este script.');
end

%% ── f_chew individual (solo para referencia en figura) ───────────────────────
tmp      = load(FILE_CHEW, 'T_freq');
T_freq   = tmp.T_freq;
f_chew_v = nan(N_CASES, 1);
for s_ = 1:N_CASES
    idx_ = strcmp(T_freq.Sujeto, CASES{s_});
    if any(idx_)
        f_chew_v(s_) = mean([T_freq.Freq_Left(idx_), T_freq.Freq_Right(idx_)], 'omitnan');
    end
end
f_chew_mean = mean(f_chew_v, 'omitnan');
fprintf('f_chew individual: M=%.3f ± %.3f Hz (N=%d)\n', ...
    f_chew_mean, std(f_chew_v,'omitnan'), sum(~isnan(f_chew_v)));

%% ── Pool paralelo ────────────────────────────────────────────────────────────
if isempty(gcp('nocreate'))
    try, parpool('local', min(8, feature('numcores'))); catch; end
end
% Replicar el path del cliente (ya con eeglab nogui) a los workers.
% NO usar genpath(eeglab) — añade octave-compat que shadowa builtins y atasca parfor.
cur_path = path();
try, spmd, addpath(cur_path); end, catch; end

%% ── Variables broadcast (evitar overhead de closures en parfor) ──────────────
cases_l    = CASES;
roi_l      = ROI_CBPT;          % 18 electrodos — fuente única S0_config
pac_dir_l  = DATA_PAC;
suf_l      = CONT_SUFFIX;
emg_l      = EMG_CHANS;
fp_l       = F_PHASE;
fa_l       = F_AMP;
np_l       = N_P;
na_l       = N_A;
bwp_l      = BW_P;
bwa_l      = BW_A;
nb_l       = N_BINS;
ns_l       = N_SURR;
trim_l     = TRIM;
min_sh_s_l = MIN_SHIFT_S;

%% ── Pre-allocate ─────────────────────────────────────────────────────────────
z_maps_out    = nan(N_P, N_A, N_CASES);
MI_obs_out    = nan(N_P, N_A, N_CASES);

fprintf('\nIniciando parfor: %d sujetos × %d fases × %d amps × %d surr\n', ...
    N_CASES, N_P, N_A, N_SURR);
fprintf('Estimado: ~%d llamadas vectorizadas tort_mi_batch por sujeto\n\n', N_P*(1+N_SURR));

%% ── Loop principal (parfor sobre sujetos) ────────────────────────────────────
parfor s = 1:N_CASES
    fname = fullfile(pac_dir_l, [cases_l{s} suf_l]);
    if ~exist(fname, 'file')
        fprintf('  [SKIP] %s no encontrado\n', cases_l{s});
        continue
    end

    EEG = pop_loadset(fname); %#ok<PFBNS>
    fs  = EEG.srate;
    T   = EEG.pnts;
    nyq = fs / 2;

    % ── EMG: criterio único bilateral + fallback (R6, igual que S4b) ───────
    ch_emg = emg_l(emg_l <= EEG.nbchan);
    if isempty(ch_emg), continue; end
    emg_raw = emg_bilateral(EEG.data, ch_emg, fs);   % [1 × T]

    % ── ROI EEG: promedio de 18 electrodos (ROI_CBPT) ─────────────────────
    n_eeg = min(64, EEG.nbchan);
    [~, roi_idx] = ismember(roi_l, {EEG.chanlocs(1:n_eeg).labels});
    roi_idx(roi_idx == 0) = [];
    if isempty(roi_idx)
        fprintf('  [SKIP] %s — ROI no encontrado en chanlocs\n', cases_l{s});
        continue
    end
    eeg_roi = mean(double(EEG.data(roi_idx, :)), 1);  % [1 × T]

    % ── Pre-computar señales de amplitud: UNA VEZ por sujeto ─────────────
    % [N_A × T] — surrogates solo shiftean esta matriz, sin re-filtrar
    amp_all = zeros(na_l, T);
    for ia = 1:na_l
        fa   = fa_l(ia);
        lo_a = max(1.0,       fa - bwa_l);
        hi_a = min(nyq - 0.5, fa + bwa_l);
        [b_a, a_a] = butter(4, [lo_a hi_a] / nyq, 'bandpass');
        amp_all(ia, :) = abs(hilbert(filtfilt(b_a, a_a, eeg_roi)));
    end

    % ── Ventana útil (descartar bordes por artefacto de filtro) ──────────
    trim_s = min(trim_l, floor(T * 0.05));
    t1 = trim_s + 1;
    t2 = T - trim_s;
    if (t2 - t1) < fs * 5   % mínimo 5 s útiles
        fprintf('  [SKIP] %s — señal demasiado corta tras trim\n', cases_l{s});
        continue
    end
    T_trim   = t2 - t1 + 1;
    amp_trim = amp_all(:, t1:t2);      % [N_A × T_trim]

    % ── Pre-generar shifts (un shift por surrogate; mismo para todas las fa) ─
    min_sh = round(min_sh_s_l * fs);
    max_sh = T_trim - min_sh;
    if max_sh <= min_sh
        fprintf('  [SKIP] %s — señal demasiado corta para surrogates\n', cases_l{s});
        continue
    end
    shifts = randi([min_sh, max_sh], ns_l, 1);  % [N_SURR × 1]

    % ── Loop sobre frecuencias de fase ────────────────────────────────────
    loc_MI_obs = nan(np_l, na_l);
    loc_z_map  = nan(np_l, na_l);

    for ip = 1:np_l
        fp   = fp_l(ip);
        lo_p = max(0.05, fp - bwp_l);
        hi_p = min(nyq  - 0.5, fp + bwp_l);

        [b_p, a_p] = butter(2, [lo_p hi_p] / nyq, 'bandpass');
        phase_full = angle(hilbert(filtfilt(b_p, a_p, emg_raw)));
        ph_trim    = phase_full(t1:t2);    % [1 × T_trim]

        % MI observado — vectorizado sobre N_A
        loc_MI_obs(ip, :) = tort_mi_batch(ph_trim, amp_trim, nb_l)';

        % Null: N_SURR circular-shifts de amplitud (fase fija)
        mi_surr = nan(ns_l, na_l);
        for k = 1:ns_l
            amp_surr  = circshift(amp_trim, shifts(k), 2);  % [N_A × T_trim]
            mi_surr(k, :) = tort_mi_batch(ph_trim, amp_surr, nb_l)';
        end

        % z-score celda a celda para este fp
        null_m = mean(mi_surr, 1);           % [1 × N_A]
        null_s = std(mi_surr, 0, 1);         % [1 × N_A]
        loc_z_map(ip, :) = (loc_MI_obs(ip,:) - null_m) ./ (null_s + eps);
    end

    z_maps_out(:, :, s)  = loc_z_map;
    MI_obs_out(:, :, s)  = loc_MI_obs;
    fprintf('  [%d/%d] %s OK\n', s, N_CASES, cases_l{s});
end

%% ── Conteo de sujetos válidos ────────────────────────────────────────────────
valid_mask = ~isnan(squeeze(z_maps_out(1, 1, :)));
n_valid    = sum(valid_mask);
fprintf('\nSujetos con datos válidos: %d / %d\n', n_valid, N_CASES);

z_valid  = z_maps_out(:, :, valid_mask);   % [N_P × N_A × n_valid]
MI_valid = MI_obs_out(:, :, valid_mask);

%% ── Estadísticas grupales celda a celda ──────────────────────────────────────
% t-test one-sample (z > 0, right-tail) para cada (fp, fa)
t_map = nan(N_P, N_A);
p_map = nan(N_P, N_A);

for ip = 1:N_P
    for ia = 1:N_A
        z_v = squeeze(z_valid(ip, ia, :));
        if sum(~isnan(z_v)) < 5, continue; end
        [~, p_map(ip,ia), ~, stats] = ttest(z_v, 0, 'tail', 'right');
        t_map(ip, ia) = stats.tstat;
    end
end

%% ── CBPT 2D: sign-flip permutation ──────────────────────────────────────────
fprintf('\nCBPT 2D (sign-flip, N_PERM=%d)...\n', N_PERM);

t_crit = tinv(1 - ALPHA_VOXEL, n_valid - 1);
fprintf('t_crit (alpha_voxel=%.2f, df=%d) = %.3f\n', ALPHA_VOXEL, n_valid-1, t_crit);

% Clusters observados
mask_obs = t_map > t_crit;
mask_obs(isnan(mask_obs)) = false;
[clust_label, n_clust] = bwlabel(mask_obs, 4);   % 4-connectivity en (fp × fa)

clust_stats_obs = nan(n_clust, 1);
for c = 1:n_clust
    clust_stats_obs(c) = sum(t_map(clust_label == c));
end
max_clust_obs = max([0; clust_stats_obs]);
fprintf('Clusters observados: %d  |  max_stat = %.1f\n', n_clust, max_clust_obs);

% Null sign-flip
% z_valid ya es (MI_obs - mean_surr)/std_surr → bajo H0 la media es ≈ 0
% El sign-flip es válido porque bajo H0, z_i y -z_i son equiprobables.
max_clust_null = nan(N_PERM, 1);

for perm = 1:N_PERM
    signs  = 2 * (rand(1, 1, n_valid) > 0.5) - 1;   % [1×1×n_valid] ±1
    z_flip = z_valid .* signs;                        % [N_P × N_A × n_valid]

    sd_flip = std(z_flip, 0, 3);
    sd_flip(sd_flip < eps) = eps;
    t_perm  = mean(z_flip, 3) ./ (sd_flip / sqrt(n_valid));

    mask_perm = t_perm > t_crit;
    mask_perm(isnan(mask_perm)) = false;
    [cp, np_c] = bwlabel(mask_perm, 4);
    if np_c == 0
        max_clust_null(perm) = 0;
    else
        cs = arrayfun(@(c) sum(t_perm(cp == c)), 1:np_c);
        max_clust_null(perm) = max(cs);
    end
end

% p-valor CBPT por cluster y máscara de significancia
p_clusters = nan(n_clust, 1);
sig_mask   = false(N_P, N_A);
fprintf('\n=== Resultados CBPT 2D ===\n');
for c = 1:n_clust
    p_clusters(c) = mean(max_clust_null >= clust_stats_obs(c));
    sig_str = '';
    if p_clusters(c) < 0.05,  sig_str = ' *';   end
    if p_clusters(c) < 0.01,  sig_str = ' **';  end
    if p_clusters(c) < 0.001, sig_str = ' ***'; end
    if p_clusters(c) < 0.05
        sig_mask(clust_label == c) = true;
        % Rango de frecuencias del cluster
        [ip_c, ia_c] = find(clust_label == c);
        fprintf('  Cluster %d: stat=%.1f  p_cbpt=%.4f%s  fp=[%.1f-%.1f]Hz  fa=[%d-%d]Hz\n', ...
            c, clust_stats_obs(c), p_clusters(c), sig_str, ...
            F_PHASE(min(ip_c)), F_PHASE(max(ip_c)), F_AMP(min(ia_c)), F_AMP(max(ia_c)));
    else
        fprintf('  Cluster %d: stat=%.1f  p_cbpt=%.4f  (n.s.)\n', ...
            c, clust_stats_obs(c), p_clusters(c));
    end
end
if ~any(sig_mask(:))
    fprintf('  [!] Ningún cluster alcanza p_cbpt < 0.05\n');
end

%% ── Guardar ──────────────────────────────────────────────────────────────────
MI_mean = mean(MI_valid, 3, 'omitnan');
z_mean  = mean(z_valid,  3, 'omitnan');

prov = struct('script',mfilename('fullpath'),'date',datestr(now,'yyyy-mm-dd HH:MM:SS'),...
    'matlab',version,'rng_seed',RNG_SEED,'subjects',{CASES},'roi',{ROI_CBPT},...
    'f_phase',F_PHASE,'f_amp',F_AMP,'bw_phase',BW_P,'bw_amp',BW_A,...
    'n_surr',N_SURR,'n_perm',N_PERM,'min_shift_s',MIN_SHIFT_S,'alpha_voxel',ALPHA_VOXEL,...
    'emg_criterion','bilateral-average (R6)','null','circular-shift amplitude',...
    'n_valid',n_valid,'n_clusters',n_clust,'framing','DESCRIPTIVO (C2/A4)');
save(fullfile(OUT_STATS, 'S4h_comodulogram_stats.mat'), ...
    'z_maps_out','MI_obs_out','t_map','p_map','sig_mask', ...
    'clust_label','clust_stats_obs','p_clusters','max_clust_null', ...
    'F_PHASE','F_AMP','N_SURR','N_PERM','ALPHA_VOXEL','t_crit', ...
    'f_chew_v','f_chew_mean','n_valid','MI_mean','z_mean','prov', '-v7.3');
fid=fopen(fullfile(OUT_STATS,'S4h_comodulogram_stats.provenance.json'),'w');
fprintf(fid,'%s',jsonencode(prov)); fclose(fid);
fprintf('\nGuardado: S4h_comodulogram_stats.mat (+ provenance.json)\n');

%% ── Figuras ──────────────────────────────────────────────────────────────────
plot_specs = {
    MI_mean, 'MI (Tort, mean)', 'S4h_MI_cluster', 1;
    z_mean,  'z(MI vs null, mean)', 'S4h_zmap_cluster', 0
};

for fig_i = 1:size(plot_specs, 1)
    data_plot = plot_specs{fig_i, 1};
    clbl      = plot_specs{fig_i, 2};
    fsuf      = plot_specs{fig_i, 3};
    use_pct98 = plot_specs{fig_i, 4};

    fh = figure('Units','inches','Position',[1 1 6.5 5.5],'Color','w');
    ax = axes(fh);
    imagesc(ax, F_AMP, F_PHASE, data_plot);
    axis(ax, 'xy');
    colormap(ax, parula);
    cb = colorbar(ax);
    cb.Label.String   = clbl;
    cb.Label.FontSize = 10;

    if use_pct98
        v = data_plot(~isnan(data_plot(:)));
        if numel(v) > 1 && range(v) > 0
            set(ax, 'CLim', [min(v), prctile(v, 98)]);
        end
    end

    hold(ax, 'on');

    % Contorno cluster significativo
    if any(sig_mask(:))
        contour(ax, F_AMP, F_PHASE, double(sig_mask), [0.5 0.5], ...
            'w-', 'LineWidth', 2.5);
    end

    % Línea f_chew media (referencia)
    if f_chew_mean >= F_PHASE(1) && f_chew_mean <= F_PHASE(end)
        plot(ax, [F_AMP(1) F_AMP(end)], [f_chew_mean f_chew_mean], ...
             'r--', 'LineWidth', 1.8);
        text(ax, F_AMP(end) - 0.5, f_chew_mean + 0.06, ...
             sprintf('\\itf_{chew}\\rm = %.2f Hz', f_chew_mean), ...
             'Color','r','FontSize',8,'HorizontalAlignment','right', ...
             'VerticalAlignment','bottom');
    end

    % Bordes de banda de amplitud (referencia)
    band_edges = [4 7; 8 13; 13 30];
    band_names = {'\theta', '\alpha', '\beta'};
    lc = [0.85 0.85 0.85];
    for b = 1:3
        for e = 1:2
            plot(ax, [band_edges(b,e) band_edges(b,e)], [F_PHASE(1) F_PHASE(end)], ...
                 '--', 'Color', [lc 0.55], 'LineWidth', 0.9);
        end
        text(ax, mean(band_edges(b,:)), F_PHASE(end) - 0.09, band_names{b}, ...
             'Color','w','FontSize',10,'FontWeight','bold', ...
             'HorizontalAlignment','center','VerticalAlignment','top');
    end

    % Rango masticatorio (referencia fase)
    for yr = [0.5, 2.0]
        if yr >= F_PHASE(1) && yr <= F_PHASE(end)
            plot(ax, [F_AMP(1) F_AMP(end)], [yr yr], ...
                 '--', 'Color', [lc 0.4], 'LineWidth', 0.9);
        end
    end

    xlabel(ax, 'Amplitude frequency (Hz)', 'FontSize', 11);
    ylabel(ax, 'Phase frequency (Hz)',      'FontSize', 11);
    title(ax, sprintf('PAC Comodulogram — Cases (N=%d)  [%s]', n_valid, clbl), ...
          'FontSize', 10, 'FontWeight', 'bold');
    ax.XTick      = [4 7 8 13 15 20 25 30];
    ax.YTick      = F_PHASE;
    ax.YTickLabel = arrayfun(@(x) sprintf('%.1f', x), F_PHASE, 'UniformOutput', false);
    ax.TickDir    = 'out';
    ax.FontSize   = 9;
    ax.Box        = 'off';
    ax.XLim       = [F_AMP(1)-0.3, F_AMP(end)+0.3];
    ax.YLim       = [F_PHASE(1)-0.05, F_PHASE(end)+0.05];

    print(fh, fullfile(OUT_FIGS, fsuf), '-dpng', '-r300');
    close(fh);
    fprintf('Guardada: %s.png\n', fsuf);
end

fprintf('\n=== FIN S4h_comodulogram_stats ===\n');
