%% ============================================================================
%  S3b_FOOOF_fromPython.m  —  FOOOF via Python (specparam) para P2V1
%  ----------------------------------------------------------------------------
%  Mismo approach que S6_FOOOF_V3.m (paper2_v3):
%   1. MATLAB: calcula Welch PSD (ROI frontal) por sujeto × condición
%   2. MATLAB: guarda PSDs en .mat temporal
%   3. Python:  corre fooof_fit.py (specparam) via system()
%   4. MATLAB: carga resultados → FOOOF_Workspace_V1.mat
%   5. MATLAB: genera figuras al estilo Figure03_FOOOF
%
%  Figuras:
%    Fig A  — PSD broadband (media ± SEM) + curva aperiódica punteada, 4 grupos
%    Fig B  — Exponente aperiódico: Casos Ch vs Nc (boxplot + scatter paired)
%    Fig C  — Δ Residual periódico (Ch−Nc): Casos vs Controles (espectro)
%    Panel  — Panel combinado A+B (una figura, dos paneles)
%
%  Correr desde Analysis_V1_Final/ con EEGLAB y fooof_fit.py en code/
% ============================================================================

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

fprintf('\n=== S3b_FOOOF_fromPython.m ===\n');

%% ── 0. CONFIGURACIÓN FOOOF ───────────────────────────────────────────────
% !! Ajustar esta ruta al Python que tiene specparam/fooof instalado !!
PYTHON_EXE  = 'python';   % o 'C:\miniconda3\python.exe', etc.
FOOOF_PY    = fullfile(ROOT_FINAL, 'code', 'fooof_fit.py');
TMP_DIR     = fullfile(OUT_STATS, 'fooof_tmp');
if ~exist(TMP_DIR,'dir'), mkdir(TMP_DIR); end

% Parámetros PSD
NFFT        = 2 * FS;        % = 1000 @ 500Hz → resolución 0.5 Hz
WIN_TASK_MS = [0, 1500];     % período de tarea (evita baseline)
FIT_RANGE   = [3, 35];       % Hz para el ajuste FOOOF

% ROI para FOOOF: usar electrodos de la INTERACCION theta (S2c sig_int_fdr)
% OJO: roi_cbpt.mat tiene el efecto principal de Casos (AF3/Fpz/AFz), NO la interacción.
% La interacción theta tiene 18 electrodos frontales: Fp1, AF7, AF3, F1, F3, F5, F7,
% FC1, Fpz, Fp2, AF8, AF4, AFz, Fz, F2, F4, FC2, FCz (de S2c sig_int_fdr).
%
% Jerarquía: S2c_stats con chanlocs → fallback hardcoded interacción

% Intento 1: cargar electrodos de interacción desde S2c + chanlocs de referencia
f_s2c = fullfile(OUT_STATS, 'S2c_TF_GroupFigure.mat');
ROI_FOOOF = {};
if exist(f_s2c,'file')
    try
        S2c = load(f_s2c, 'sig_int_fdr');
        % Cargar chanlocs desde primera época disponible
        ep_ref = fullfile(DIR_EPOCHS, [CASES{1} EP_SUFFIX_NC]);
        if exist(ep_ref,'file')
            EEG_ref = pop_loadset('filename',[CASES{1} EP_SUFFIX_NC], 'filepath',DIR_EPOCHS);
            sig_idx = find(S2c.sig_int_fdr(:)' == 1);
            ROI_FOOOF = {EEG_ref.chanlocs(sig_idx).labels};
            fprintf('ROI desde S2c sig_int_fdr (%d electrodos): %s\n', ...
                numel(ROI_FOOOF), strjoin(ROI_FOOOF, ', '));
        end
    catch ME
        fprintf('[WARN] No se pudo cargar ROI desde S2c: %s\n', ME.message);
    end
end

% Fallback: hardcode de los 18 electrodos de interacción theta (verificados 2026-06-18)
if isempty(ROI_FOOOF)
    ROI_FOOOF = {'Fp1','AF7','AF3','F1','F3','F5','F7','FC1', ...
                 'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','FC2','FCz'};
    fprintf('[FALLBACK] ROI interacción theta hardcoded (%d ch): %s\n', ...
        numel(ROI_FOOOF), strjoin(ROI_FOOOF, ', '));
end

%% ── Verificar Python + specparam ─────────────────────────────────────────
[st1, ~] = system(sprintf('"%s" -c "import specparam" 2>&1', PYTHON_EXE));
[st2, ~] = system(sprintf('"%s" -c "import fooof" 2>&1', PYTHON_EXE));
if st1 ~= 0 && st2 ~= 0
    error(['Python o specparam/fooof no encontrado.\n' ...
           'Python: %s\n' ...
           'Instalar: pip install specparam  (o pip install fooof)'], PYTHON_EXE);
end
assert(exist(FOOOF_PY,'file')==2, 'fooof_fit.py no encontrado en: %s', FOOOF_PY);
fprintf('Python OK. fooof_fit.py OK.\n');

%% ── 1. CALCULAR PSD POR GRUPO × CONDICIÓN ────────────────────────────────
% Estructura de salida: GR.Cases / GR.Controls
%   .PSD_Ch [nFreq × nSubj]  — PSD lineal (µV²/Hz)
%   .PSD_Nc [nFreq × nSubj]
%   .f      [nFreq × 1]
%   .n      scalar

fprintf('\n--- Calculando PSDs (Welch) ---\n');
addpath(fileparts(which('eeglab.m')));  % asegurar EEGLAB en path

GR = struct();
grp_names = {'Cases','Controls'};
suj_lists = {CASES, CONTROLS};

for gc = 1:2
    grp  = grp_names{gc};
    subs = suj_lists{gc};
    nS   = numel(subs);
    fprintf('\n=== %s (N=%d) ===\n', grp, nS);

    % Obtener referencia para freqs, srate, ROI_IDX
    ref_EEG = [];
    for k = 1:nS
        fp_ref = fullfile(DIR_EPOCHS, [subs{k} EP_SUFFIX_NC]);
        if exist(fp_ref,'file')
            ref_EEG = pop_loadset('filename',[subs{k} EP_SUFFIX_NC], ...
                                  'filepath',DIR_EPOCHS);
            break;
        end
    end
    if isempty(ref_EEG)
        warning('%s: no se encontró ningún archivo de referencia. Saltando.', grp);
        continue;
    end

    srate  = ref_EEG.srate;
    times  = ref_EEG.times;
    t_task = find(times >= WIN_TASK_MS(1) & times <= WIN_TASK_MS(2));
    if isempty(t_task)
        error('WIN_TASK_MS=[%d,%d] ms no tiene índices en las épocas', WIN_TASK_MS);
    end

    % ROI por etiqueta (robusto a reordenación de canales)
    all_labels = {ref_EEG.chanlocs.labels};
    ROI_IDX    = find(ismember(all_labels, ROI_FOOOF));
    if isempty(ROI_IDX)
        warning('%s: ningún canal del ROI encontrado. Usando primeros 4 frontales.', grp);
        ROI_IDX = 1:min(4, EEG_N);
    end
    ROI_used = all_labels(ROI_IDX);
    fprintf('  ROI usado: %s\n', strjoin(ROI_used, ', '));

    % Template de frecuencias (ventana no puede exceder la longitud de un trial)
    nfft_use = min(NFFT, numel(t_task));
    [~, freqs] = pwelch(zeros(nfft_use * 4, 1), hann(nfft_use), floor(nfft_use/2), ...
                        nfft_use, srate);
    nFreq = numel(freqs);

    psds_ch = nan(nFreq, nS);
    psds_nc = nan(nFreq, nS);

    for si = 1:nS
        suj = subs{si};
        for ic = 1:2
            if ic==1, suf = EP_SUFFIX;    cond = 'Ch';
            else,     suf = EP_SUFFIX_NC; cond = 'Nc'; end
            fp = fullfile(DIR_EPOCHS, [suj suf]);
            if ~exist(fp,'file')
                fprintf('  [SKIP] %s %s — archivo no encontrado\n', suj, cond);
                continue;
            end
            try
                EP = pop_loadset('filename',[suj suf], 'filepath',DIR_EPOCHS);
                if EP.trials < 3
                    fprintf('  [SKIP] %s %s — trials=%d\n', suj, cond, EP.trials);
                    continue;
                end
                % ROI_IDX puede diferir por sujeto; recalcular si cambia
                lbs_s = {EP.chanlocs.labels};
                roi_s = find(ismember(lbs_s, ROI_FOOOF));
                if isempty(roi_s), roi_s = ROI_IDX; end

                % Concatenar trials → vector continuo
                seg = [];
                for e = 1:EP.trials
                    roi_avg = mean(double(EP.data(roi_s, t_task, e)), 1);
                    seg = [seg, roi_avg(:)']; %#ok<AGROW>
                end
                [pxx, ~] = pwelch(seg(:), hann(nfft_use), floor(nfft_use/2), ...
                                  nfft_use, srate);
                if ic==1, psds_ch(:,si) = pxx;
                else,      psds_nc(:,si) = pxx; end
            catch ME
                fprintf('  [WARN] %s %s: %s\n', suj, cond, ME.message);
            end
        end
        fprintf('  %-10s Ch=%s Nc=%s\n', suj, ...
            ternario(~isnan(psds_ch(1,si)),'OK','---'), ...
            ternario(~isnan(psds_nc(1,si)),'OK','---'));
    end

    GR.(grp).f       = freqs;
    GR.(grp).PSD_Ch  = psds_ch;
    GR.(grp).PSD_Nc  = psds_nc;
    GR.(grp).n       = nS;
    GR.(grp).roi     = ROI_used;
end

%% ── 2. LLAMAR PYTHON PARA EL AJUSTE FOOOF ────────────────────────────────
fprintf('\n--- Corriendo Python specparam ---\n');

for gc = 1:2
    grp = grp_names{gc};
    if ~isfield(GR, grp), continue; end
    freqs = GR.(grp).f;

    for ic = 1:2
        if ic==1, cond = 'Ch'; psds = GR.(grp).PSD_Ch;
        else,     cond = 'Nc'; psds = GR.(grp).PSD_Nc; end

        % Excluir sujetos sin datos
        valid = all(~isnan(psds), 1);
        if ~any(valid), warning('%s/%s: sin PSDs válidas', grp, cond); continue; end
        psds_v = psds(:, valid);

        mat_in  = fullfile(TMP_DIR, sprintf('%s_%s_psd_in.mat',  grp, cond));
        mat_out = fullfile(TMP_DIR, sprintf('%s_%s_psd_out.mat', grp, cond));

        save(mat_in, 'psds_v', 'freqs', '-v7.3');
        % Renombrar psds_v → psds para que fooof_fit.py lo encuentre
        % (fooof_fit.py busca variable 'psds')
        tmp = load(mat_in); tmp.psds = tmp.psds_v; tmp = rmfield(tmp,'psds_v');
        save(mat_in, '-struct', 'tmp', '-v7.3');

        cmd = sprintf('"%s" "%s" --input "%s" --output "%s" --fit_low %.1f --fit_high %.1f', ...
            PYTHON_EXE, FOOOF_PY, mat_in, mat_out, FIT_RANGE(1), FIT_RANGE(2));
        fprintf('\n--- %s / %s ---\n', grp, cond);
        [st, out_txt] = system(cmd);
        fprintf('%s\n', out_txt);
        if st ~= 0
            warning('fooof_fit.py error %d para %s/%s', st, grp, cond);
            continue;
        end
        if ~exist(mat_out,'file')
            warning('No se generó %s', mat_out); continue;
        end

        R = load(mat_out);
        nAll = GR.(grp).n;
        % Rellenar NaN donde los sujetos no tenían datos
        exp_full = nan(nAll,1);   off_full = nan(nAll,1);
        r2_full  = nan(nAll,1);
        ap_full  = nan(numel(freqs), nAll);
        res_full = nan(numel(freqs), nAll);
        psd_log_full = nan(numel(freqs), nAll);

        vi = find(valid);
        exp_full(vi)    = R.exponents(:);
        off_full(vi)    = R.offsets(:);
        r2_full(vi)     = R.r_squared(:);
        ap_full(:,vi)   = R.ap_fits;
        res_full(:,vi)  = R.residuals;
        psd_log_full(:,vi) = R.psd_log;

        if ic==1
            GR.(grp).exp_Ch = exp_full; GR.(grp).off_Ch = off_full;
            GR.(grp).r2_Ch  = r2_full;  GR.(grp).AP_Ch  = ap_full;
            GR.(grp).Res_Ch = res_full;  GR.(grp).PSD_log_Ch = psd_log_full;
        else
            GR.(grp).exp_Nc = exp_full; GR.(grp).off_Nc = off_full;
            GR.(grp).r2_Nc  = r2_full;  GR.(grp).AP_Nc  = ap_full;
            GR.(grp).Res_Nc = res_full;  GR.(grp).PSD_log_Nc = psd_log_full;
        end
    end
end

%% ── 3. GUARDAR WORKSPACE ──────────────────────────────────────────────────
fw = fullfile(OUT_STATS, 'FOOOF_Workspace_V1.mat');
save(fw, 'GR', 'FIT_RANGE', 'WIN_TASK_MS', 'ROI_FOOOF', '-v7.3');
fprintf('\n>>> FOOOF_Workspace_V1.mat guardado: %s\n', fw);

%% ── 4. ESTADÍSTICOS ──────────────────────────────────────────────────────
fprintf('\n--- Estadísticos exponent ---\n');

exp_cas_ch = GR.Cases.exp_Ch;
exp_cas_nc = GR.Cases.exp_Nc;
exp_ctr_ch = GR.Controls.exp_Ch;
exp_ctr_nc = GR.Controls.exp_Nc;

valid_cas = ~isnan(exp_cas_ch) & ~isnan(exp_cas_nc);
valid_ctr = ~isnan(exp_ctr_ch) & ~isnan(exp_ctr_nc);

[p_cas, ~]  = signrank(exp_cas_ch(valid_cas), exp_cas_nc(valid_cas));
[p_ctr, ~]  = signrank(exp_ctr_ch(valid_ctr), exp_ctr_nc(valid_ctr));
d_exp_cas   = exp_cas_nc(valid_cas) - exp_cas_ch(valid_cas);  % Nc-Ch: pos=aplanamiento

fprintf('  Cases  Nc: %.3f ± %.3f\n', nanmean(exp_cas_nc), nanstd(exp_cas_nc));
fprintf('  Cases  Ch: %.3f ± %.3f   Wilcoxon p=%.4f %s\n', ...
    nanmean(exp_cas_ch), nanstd(exp_cas_ch), p_cas, p2s(p_cas));
fprintf('  Controls Nc: %.3f ± %.3f\n', nanmean(exp_ctr_nc), nanstd(exp_ctr_nc));
fprintf('  Controls Ch: %.3f ± %.3f   Wilcoxon p=%.4f %s\n', ...
    nanmean(exp_ctr_ch), nanstd(exp_ctr_ch), p_ctr, p2s(p_ctr));

%% ── 5. FIGURA A: PSD BROADBAND + FIT APERIÓDICO ──────────────────────────
fprintf('\n--- Generando figuras ---\n');
freqs_c = GR.Cases.f;
freqs_k = GR.Controls.f;
assert(numel(freqs_c)==numel(freqs_k) && max(abs(freqs_c-freqs_k))<0.01, ...
    'Vectors de frecuencia inconsistentes');
f = freqs_c;
f_mask = f >= FIT_RANGE(1) & f <= FIT_RANGE(2);

% Promediar en dB: 10*log10(PSD)
dB = @(x) 10*log10(max(x, eps));
psd_db = @(M) dB(M);   % [nFreq × nSubj] → misma escala

mn = @(M) nanmean(M, 2)';
se = @(M) nanstd(M,0,2)' ./ sqrt(sum(~isnan(M),2)');

cas_ch_psd = psd_db(GR.Cases.PSD_Ch);   cas_nc_psd = psd_db(GR.Cases.PSD_Nc);
ctr_ch_psd = psd_db(GR.Controls.PSD_Ch);ctr_nc_psd = psd_db(GR.Controls.PSD_Nc);

% No se normaliza — se grafica la PSD en dB tal como está (sin alineación artificial)

% Fits aperiódicos promedio (ya en log10 µV²/Hz)
ap_cas_ch = GR.Cases.AP_Ch;   ap_cas_nc = GR.Cases.AP_Nc;
ap_ctr_ch = GR.Controls.AP_Ch;ap_ctr_nc = GR.Controls.AP_Nc;

figA = figure('Name','FOOOF PSD','Units','inches','Position',[1 1 7 5],'Color','w');
axA  = axes(figA); hold(axA,'on');

% Bandas de frecuencia (fondo gris)
band_ranges = {BAND_THETA, BAND_ALPHA, BAND_BETA};
band_names  = {'\theta', '\alpha', '\beta'};
y_lim_tmp = [-18, 2];  % límites provisionales, se reajustan
for bi = 1:3
    xr = band_ranges{bi};
    patch(axA, [xr(1) xr(2) xr(2) xr(1)], ...
          [y_lim_tmp(1) y_lim_tmp(1) y_lim_tmp(2) y_lim_tmp(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
    text(axA, mean(xr), y_lim_tmp(2)-0.4, band_names{bi}, ...
        'HorizontalAlignment','center','FontSize',FIG_FS-1,'Color',[0.5 0.5 0.5]);
end

% Media ± SEM con relleno
clr_cas = COL_CASE;  clr_ctr = COL_CTRL;

fill_sem(axA, f, mn(cas_nc_psd), se(cas_nc_psd), clr_ctr, 0.20);
fill_sem(axA, f, mn(ctr_nc_psd), se(ctr_nc_psd), clr_ctr, 0.10);
fill_sem(axA, f, mn(cas_ch_psd), se(cas_ch_psd), clr_cas, 0.20);
fill_sem(axA, f, mn(ctr_ch_psd), se(ctr_ch_psd), clr_cas, 0.10);

hL(1) = plot(axA, f, mn(ctr_nc_psd), '-',  'Color', clr_ctr,    'LineWidth',FIG_LW);
hL(2) = plot(axA, f, mn(ctr_ch_psd), '--', 'Color', clr_ctr,    'LineWidth',FIG_LW);
hL(3) = plot(axA, f, mn(cas_nc_psd), '-',  'Color', clr_cas,    'LineWidth',FIG_LW);
hL(4) = plot(axA, f, mn(cas_ch_psd), '--', 'Color', clr_cas,    'LineWidth',FIG_LW);

% Curvas aperiódicas (punteadas)
ap_to_db = @(A) 10 * A;   % ap_fits es log10 en la misma escala dB
% (fooof_fit.py devuelve residuals en log10 power, que equivale a dB/10)
plot(axA, f(f_mask), 10*nanmean(ap_cas_nc(f_mask,:),2), ':', ...
    'Color', clr_cas*0.7+0.3, 'LineWidth',1,'HandleVisibility','off');
plot(axA, f(f_mask), 10*nanmean(ap_cas_ch(f_mask,:),2), ':', ...
    'Color', clr_cas,   'LineWidth',1,'HandleVisibility','off');
plot(axA, f(f_mask), 10*nanmean(ap_ctr_nc(f_mask,:),2), ':', ...
    'Color', clr_ctr*0.7+0.3, 'LineWidth',1,'HandleVisibility','off');
plot(axA, f(f_mask), 10*nanmean(ap_ctr_ch(f_mask,:),2), ':', ...
    'Color', clr_ctr,   'LineWidth',1,'HandleVisibility','off');

hLdot = plot(axA, NaN, NaN, 'k:', 'LineWidth',1);  % leyenda punteada

legend([hL hLdot], {'Controls – No Chew','Controls – Chew', ...
    'Cases – No Chew','Cases – Chew','Dotted: aperiodic fit'}, ...
    'Location','southwest','FontSize',FIG_FS-1,'Box','off');

xlim(axA, [FIT_RANGE(1) FIT_RANGE(2)]);
xlabel(axA, 'Frequency (Hz)', 'FontSize',FIG_FS);
ylabel(axA, 'Power (dB)',     'FontSize',FIG_FS);
title(axA,  sprintf('Broadband PSD (pre-FOOOF) — ROI: %s', strjoin(ROI_FOOOF,', ')), ...
    'FontSize',FIG_FS+1,'FontWeight','bold');
set(axA,'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS);

exportgraphics(figA, fullfile(OUT_FIGS,'S3b_FOOOF_PSD.png'), 'Resolution',FIG_DPI);
close(figA);
fprintf('  Fig A guardada.\n');

%% ── 6. FIGURA B: EXPONENT BOXPLOT (CASOS Ch vs Nc) ──────────────────────
figB = figure('Name','FOOOF Exponent','Units','inches','Position',[1 1 5 6],'Color','w');
axB  = axes(figB);  hold(axB,'on');

nCas  = N_CASES;
xNc   = 1; xCh = 2;
jit   = 0.12;   bw = 0.40;
clrNc_box = COL_NC;  clrCh_box = clr_cas;

% Cajas IQR
draw_box(axB, exp_cas_nc, xNc, bw, clrNc_box, 0.55);
draw_box(axB, exp_cas_ch, xCh, bw, clrCh_box, 0.55);

% Líneas pareadas
rng(RNG_SEED);
for i = 1:nCas
    if valid_cas(i)
        line(axB, [xNc xCh], [exp_cas_nc(i) exp_cas_ch(i)], ...
             'Color',[0.5 0.5 0.5 0.30],'LineWidth',0.8);
    end
end

% Scatter individual
scatter(axB, xNc + randn(nCas,1)*jit, exp_cas_nc, 28, clrNc_box, ...
        'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
scatter(axB, xCh + randn(nCas,1)*jit, exp_cas_ch, 28, clrCh_box, ...
        'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');

% Media ± SEM
sem_nc_exp = nanstd(exp_cas_nc)/sqrt(sum(valid_cas));
sem_ch_exp = nanstd(exp_cas_ch)/sqrt(sum(valid_cas));
errorbar(axB, xNc, nanmean(exp_cas_nc), sem_nc_exp, 'k^', ...
         'MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',6);
errorbar(axB, xCh, nanmean(exp_cas_ch), sem_ch_exp, 'k^', ...
         'MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',6);

% Bracket significancia
ymax_exp = max([exp_cas_nc(valid_cas); exp_cas_ch(valid_cas)]) + 0.05;
yb = ymax_exp + 0.05;
line(axB, [xNc xNc xCh xCh], [yb-0.03 yb yb yb-0.03],'Color','k','LineWidth',1.2);
text(axB, mean([xNc xCh]), yb+0.03, p2s(p_cas), ...
     'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');

set(axB,'XTick',[xNc xCh],'XTickLabel',{'No Chew','Chew'},'XLim',[0.4 2.6], ...
        'YGrid','on','GridAlpha',0.25,'GridColor',[0.6 0.6 0.6], ...
        'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS+1);
ylabel(axB, 'Aperiodic Exponent', 'FontSize',FIG_FS+2);
title(axB, sprintf('FOOOF Exponent — Cases (N=%d)', N_CASES), ...
    'FontSize',FIG_FS+2,'FontWeight','bold');

exportgraphics(figB, fullfile(OUT_FIGS,'S3b_FOOOF_Exponent.png'), 'Resolution',FIG_DPI);
close(figB);
fprintf('  Fig B guardada.\n');

%% ── 7. FIGURA C: Δ RESIDUAL PERIÓDICO (Ch−Nc) ────────────────────────────
figC = figure('Name','FOOOF Delta Residual','Units','inches','Position',[1 1 7 5],'Color','w');
axC  = axes(figC);  hold(axC,'on');

% Δ residual: mean(Ch) - mean(Nc) por sujeto → promedio grupo
delta_cas = GR.Cases.Res_Ch  - GR.Cases.Res_Nc;   % [nFreq × nCas]
delta_ctr = GR.Controls.Res_Ch - GR.Controls.Res_Nc; % [nFreq × nCtr]

% Bandas de referencia
for bi = 1:3
    xr = band_ranges{bi};
    patch(axC, [xr(1) xr(2) xr(2) xr(1)], ...
          [-2.5 -2.5 2.5 2.5], [0.9 0.9 0.9],'EdgeColor','none', ...
          'FaceAlpha',0.5,'HandleVisibility','off');
    text(axC, mean(xr), 2.2, band_names{bi}, ...
        'HorizontalAlignment','center','FontSize',FIG_FS-1,'Color',[0.5 0.5 0.5]);
end

yline(axC, 0, 'k--', 'LineWidth',1.0,'HandleVisibility','off');

fill_sem(axC, f, mn(delta_ctr), se(delta_ctr), clr_ctr, 0.25);
fill_sem(axC, f, mn(delta_cas), se(delta_cas), clr_cas, 0.25);

plot(axC, f, mn(delta_ctr), '-', 'Color',clr_ctr, 'LineWidth',FIG_LW, ...
    'DisplayName','Controls');
plot(axC, f, mn(delta_cas), '-', 'Color',clr_cas, 'LineWidth',FIG_LW, ...
    'DisplayName','Cases');

% Marcar significancia por banda (Wilcoxon Δ vs 0)
bands_test = {BAND_THETA, BAND_ALPHA, BAND_BETA};
y_sig = 2.1;
for bi = 1:3
    fmask_b = f >= bands_test{bi}(1) & f <= bands_test{bi}(2);
    if ~any(fmask_b), continue; end
    delta_cas_b = nanmean(delta_cas(fmask_b,:), 1)';
    delta_ctr_b = nanmean(delta_ctr(fmask_b,:), 1)';
    [p_cas_b] = signrank(delta_cas_b(valid_cas));
    [p_ctr_b] = signrank(delta_ctr_b(valid_ctr));
    txt_cas = sprintf('Cases:%s\nCtrl:%s', p2s(p_cas_b), p2s(p_ctr_b));
    text(axC, mean(bands_test{bi}), y_sig, txt_cas, ...
        'HorizontalAlignment','center','FontSize',FIG_FS-2, ...
        'Color','k','VerticalAlignment','bottom');
end

xlim(axC, [FIT_RANGE(1) FIT_RANGE(2)]);
ylim(axC, [-2.5 2.5]);
xlabel(axC,'Frequency (Hz)','FontSize',FIG_FS);
ylabel(axC,'\Delta Periodic Power Ch−Nc (dB)','FontSize',FIG_FS);
title(axC,'Chewing-induced oscillatory change (FOOOF residual)', ...
    'FontSize',FIG_FS+1,'FontWeight','bold');
legend(axC,'Location','southeast','FontSize',FIG_FS-1,'Box','off');
set(axC,'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS);

exportgraphics(figC, fullfile(OUT_FIGS,'S3b_FOOOF_DeltaPeriodic.png'), 'Resolution',FIG_DPI);
close(figC);
fprintf('  Fig C guardada.\n');

%% ── 8. PANEL COMBINADO (A + B) ────────────────────────────────────────────
figP = figure('Name','FOOOF Panel','Units','inches','Position',[1 1 12 5],'Color','w');

% Panel izquierdo: PSD (re-dibujar)
axP1 = subplot(1,2,1,'Parent',figP);  hold(axP1,'on');
for bi=1:3
    xr=band_ranges{bi};
    patch(axP1,[xr(1) xr(2) xr(2) xr(1)],[-18 -18 2 2],[0.9 0.9 0.9],...
        'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
    text(axP1,mean(xr),-0.5,band_names{bi},'HorizontalAlignment','center',...
        'FontSize',FIG_FS-1,'Color',[0.5 0.5 0.5]);
end
fill_sem(axP1,f,mn(cas_nc_psd),se(cas_nc_psd),clr_cas,0.15);
fill_sem(axP1,f,mn(cas_ch_psd),se(cas_ch_psd),clr_cas,0.15);
fill_sem(axP1,f,mn(ctr_nc_psd),se(ctr_nc_psd),clr_ctr,0.15);
fill_sem(axP1,f,mn(ctr_ch_psd),se(ctr_ch_psd),clr_ctr,0.15);
hPL(1)=plot(axP1,f,mn(ctr_nc_psd),'-', 'Color',clr_ctr,'LineWidth',FIG_LW);
hPL(2)=plot(axP1,f,mn(ctr_ch_psd),'--','Color',clr_ctr,'LineWidth',FIG_LW);
hPL(3)=plot(axP1,f,mn(cas_nc_psd),'-', 'Color',clr_cas,'LineWidth',FIG_LW);
hPL(4)=plot(axP1,f,mn(cas_ch_psd),'--','Color',clr_cas,'LineWidth',FIG_LW);
plot(axP1,f(f_mask),10*nanmean(ap_cas_nc(f_mask,:),2),':', 'Color',clr_cas*0.6+0.4,'LineWidth',1,'HandleVisibility','off');
plot(axP1,f(f_mask),10*nanmean(ap_cas_ch(f_mask,:),2),':', 'Color',clr_cas,'LineWidth',1,'HandleVisibility','off');
plot(axP1,f(f_mask),10*nanmean(ap_ctr_nc(f_mask,:),2),':', 'Color',clr_ctr*0.6+0.4,'LineWidth',1,'HandleVisibility','off');
plot(axP1,f(f_mask),10*nanmean(ap_ctr_ch(f_mask,:),2),':', 'Color',clr_ctr,'LineWidth',1,'HandleVisibility','off');
hPdot=plot(axP1,NaN,NaN,'k:','LineWidth',1);
legend(axP1,[hPL hPdot],{'Controls–Nc','Controls–Ch','Cases–Nc','Cases–Ch','Dotted: aperiodic fit'},...
    'Location','southwest','FontSize',FIG_FS-2,'Box','off');
xlim(axP1,[FIT_RANGE(1) FIT_RANGE(2)]);
xlabel(axP1,'Frequency (Hz)','FontSize',FIG_FS);
ylabel(axP1,'Power (dB)','FontSize',FIG_FS);
title(axP1,'Broadband PSD (pre-FOOOF)','FontSize',FIG_FS+1,'FontWeight','bold');
set(axP1,'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS);

% Panel derecho: Exponent (re-dibujar)
axP2 = subplot(1,2,2,'Parent',figP); hold(axP2,'on');
draw_box(axP2,exp_cas_nc,1,0.40,clrNc_box,0.55);
draw_box(axP2,exp_cas_ch,2,0.40,clrCh_box,0.55);
rng(RNG_SEED);
for i=1:nCas
    if valid_cas(i)
        line(axP2,[1 2],[exp_cas_nc(i) exp_cas_ch(i)],'Color',[0.5 0.5 0.5 0.30],'LineWidth',0.8);
    end
end
scatter(axP2,1+randn(nCas,1)*jit,exp_cas_nc,28,clrNc_box,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
scatter(axP2,2+randn(nCas,1)*jit,exp_cas_ch,28,clrCh_box,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
errorbar(axP2,1,nanmean(exp_cas_nc),sem_nc_exp,'k^','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',6);
errorbar(axP2,2,nanmean(exp_cas_ch),sem_ch_exp,'k^','MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',6);
ymax_p2=max([exp_cas_nc(valid_cas);exp_cas_ch(valid_cas)])+0.05;
yb2=ymax_p2+0.05;
line(axP2,[1 1 2 2],[yb2-0.03 yb2 yb2 yb2-0.03],'Color','k','LineWidth',1.2);
text(axP2,1.5,yb2+0.03,p2s(p_cas),'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
set(axP2,'XTick',[1 2],'XTickLabel',{'No Chew','Chew'},'XLim',[0.4 2.6],...
    'YGrid','on','GridAlpha',0.25,'GridColor',[0.6 0.6 0.6],...
    'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS+1);
ylabel(axP2,'Aperiodic Exponent','FontSize',FIG_FS+2);
title(axP2,'Aperiodic Exponent','FontSize',FIG_FS+1,'FontWeight','bold');

exportgraphics(figP, fullfile(OUT_FIGS,'S3b_FOOOF_Panel.png'), 'Resolution',FIG_DPI);
close(figP);
fprintf('  Panel guardado.\n');

%% ── 9. REPORTE ───────────────────────────────────────────────────────────
fid = fopen(fullfile(OUT_REVIEWER,'S3b_FOOOF_result.txt'),'w');
fprintf(fid,'=== FOOOF (S3b — Python specparam) — P2V1 ===\n');
fprintf(fid,'N_casos=%d | N_controles=%d | ROI=%s\n', N_CASES, N_CONTROLS, strjoin(ROI_FOOOF,', '));
fprintf(fid,'FIT_RANGE=[%.0f,%.0f] Hz | WIN_TASK=[%.0f,%.0f] ms\n\n', ...
    FIT_RANGE, WIN_TASK_MS);
fprintf(fid,'Exponente aperiódico:\n');
fprintf(fid,'  Cases  Nc: %.3f ± %.3f\n', nanmean(exp_cas_nc), nanstd(exp_cas_nc));
fprintf(fid,'  Cases  Ch: %.3f ± %.3f   Wilcoxon p=%.4f %s\n', ...
    nanmean(exp_cas_ch), nanstd(exp_cas_ch), p_cas, p2s(p_cas));
fprintf(fid,'  Controls Nc: %.3f ± %.3f\n', nanmean(exp_ctr_nc), nanstd(exp_ctr_nc));
fprintf(fid,'  Controls Ch: %.3f ± %.3f   Wilcoxon p=%.4f %s\n', ...
    nanmean(exp_ctr_ch), nanstd(exp_ctr_ch), p_ctr, p2s(p_ctr));
fprintf(fid,'\nFiguras:\n');
fprintf(fid,'  S3b_FOOOF_PSD.png       — PSD broadband\n');
fprintf(fid,'  S3b_FOOOF_Exponent.png  — Boxplot exponent Cases\n');
fprintf(fid,'  S3b_FOOOF_DeltaPeriodic.png — Delta residual\n');
fprintf(fid,'  S3b_FOOOF_Panel.png     — Panel combinado\n');
fclose(fid);

fprintf('\n✓ S3b_FOOOF_fromPython.m completado.\n');
fprintf('Archivos en: %s\n\n', OUT_FIGS);

%% ════════════════════════════════════════════════════════════════════
%  FUNCIONES LOCALES
%% ════════════════════════════════════════════════════════════════════

function draw_box(ax, data, xc, bw, clr, alpha)
    data = data(~isnan(data));
    if isempty(data), return; end
    q   = quantile(data, [0.25 0.50 0.75]);
    iqr_v = q(3) - q(1);
    wlo = min(data(data >= q(1) - 1.5*iqr_v));
    whi = max(data(data <= q(3) + 1.5*iqr_v));
    hw  = bw/2;
    fill(ax, [xc-hw xc+hw xc+hw xc-hw xc-hw], ...
             [q(1)  q(1)  q(3)  q(3)  q(1)], ...
        clr, 'FaceAlpha',alpha, 'EdgeColor',clr*0.65,'LineWidth',1.2);
    line(ax, [xc-hw xc+hw], [q(2) q(2)], 'Color','w','LineWidth',2.5);
    line(ax, [xc xc], [q(1) wlo], 'Color',clr*0.65,'LineWidth',1.2);
    line(ax, [xc xc], [q(3) whi], 'Color',clr*0.65,'LineWidth',1.2);
    % Outliers
    out_v = data(data < wlo | data > whi);
    if ~isempty(out_v)
        scatter(ax, repmat(xc,size(out_v)), out_v, 20, clr*0.65, ...
            'filled','MarkerEdgeColor','none','HandleVisibility','off');
    end
end

function fill_sem(ax, x, mn_v, se_v, clr, alp)
    xi = [x(:)', fliplr(x(:)')];
    yi = [mn_v-se_v, fliplr(mn_v+se_v)];
    fill(ax, xi, yi, clr, 'FaceAlpha',alp,'EdgeColor','none','HandleVisibility','off');
end

function s = p2s(p)
    if     isnan(p),  s = 'n/a';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'NS';
    end
end

function out = ternario(cond, a, b)
    if cond, out = a; else, out = b; end
end
