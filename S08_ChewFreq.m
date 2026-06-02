%% ============================================================================
%  S5_ChewFreq.m  —  Frecuencia masticatoria EMG vs conducta y EEG
%  ----------------------------------------------------------------------------
%  Carga la frecuencia masticatoria individual (peak_hz, Hz) desde el
%  workspace PAC (pre-calculada en S5b_PAC_Continuous.m) y la relaciona con:
%    (a) Distribución de frecuencias en el grupo Casos (descriptivo)
%    (b) ΔIES (Ch−Nc) — correlación Spearman
%    (c) ΔzMI Beta-Early (PAC) — correlación Spearman
%
%  Fuente peak_hz: PAC_Continuous_Workspace.mat → W.peak_hz (1×30, Hz)
%  Alternativa:    chew_metrics.mat si se prefiere cargar por separado.
%
%  Figuras:
%    Fig 5A — Distribución peak_hz (histograma + violín)
%    Fig 5B — peak_hz vs ΔIES (scatter Spearman)
%    Fig 5C — peak_hz vs ΔzMI Beta-Early (scatter Spearman)
%  ============================================================================
clear; clc; close all;

%% -- 1. RUTAS ----------------------------------------------------------------
P       = S0_paths();
dir_out = P.fig05;
f_rep   = fullfile(dir_out, 'Reporte_ChewFreq.txt');

%% -- 2. ESTÉTICA -------------------------------------------------------------
clrCas = P.clr_cas;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',12,'DefaultFigureColor','w');
dpi = 300;

%% -- 3. CARGA DE DATOS -------------------------------------------------------
% 3A. Frecuencia masticatoria individual (desde workspace PAC)
fprintf('>>> Cargando peak_hz desde PAC workspace...\n');
W_pac    = load(fullfile(P.dir_pac, 'PAC_Continuous_Workspace.mat'), ...
                'peak_hz', 'B', 'casos');
peak_hz  = W_pac.peak_hz(:);   % N=30 Hz

% 3B. Conducta (IES casos)
load(P.file_beh, 'tb_data_45');
ies_ch  = double(tb_data_45.casos.chew.ies(:));
ies_nc  = double(tb_data_45.casos.nochew.ies(:));
d_ies   = ies_ch - ies_nc;    % ΔIES (positivo = empeora, negativo = mejora)

% 3C. PAC Beta-Early (Ch condition)
% Workspace: B.Ch.Beta es (30 sujetos × 3 ventanas) → col 1 = Early
zmi_beta_early_ch = W_pac.B.Ch.Beta(:,1);   % (30×1), Early window

nS = numel(peak_hz);
fprintf('    N = %d Casos | f_chew: %.2f ± %.2f Hz\n', nS, mean(peak_hz), std(peak_hz));

%% -- 4. ESTADÍSTICOS ---------------------------------------------------------
[rho_ies,  p_ies]  = corr(peak_hz, d_ies,            'type','Spearman');
[rho_pac,  p_pac]  = corr(peak_hz, zmi_beta_early_ch, 'type','Spearman');

%% -- 5. REPORTE ESCRITO ------------------------------------------------------
fid = fopen(f_rep, 'w');
fprintf(fid, '================================================================\n');
fprintf(fid, '  Frecuencia Masticatoria EMG — Resultados paper  (S5_ChewFreq.m)\n');
fprintf(fid, '  Generado: %s\n', datestr(now));
fprintf(fid, '  N = %d Casos\n', nS);
fprintf(fid, '================================================================\n\n');

fprintf(fid, '── Descriptiva peak_hz ──────────────────────────────────────\n');
fprintf(fid, '  Media  : %.3f Hz\n', mean(peak_hz));
fprintf(fid, '  Mediana: %.3f Hz\n', median(peak_hz));
fprintf(fid, '  SD     : %.3f Hz\n', std(peak_hz));
fprintf(fid, '  Rango  : [%.3f – %.3f] Hz\n\n', min(peak_hz), max(peak_hz));

fprintf(fid, '── Correlaciones Spearman ───────────────────────────────────\n');
fprintf(fid, '  peak_hz vs ΔIES (Ch−Nc):         rho=%.3f  p=%.4f %s\n', ...
        rho_ies, p_ies, p2s(p_ies));
fprintf(fid, '  peak_hz vs ΔzMI_Beta_Early (Ch): rho=%.3f  p=%.4f %s\n\n', ...
        rho_pac, p_pac, p2s(p_pac));

fclose(fid);
fprintf('  → Reporte: %s\n', f_rep);

%% -- 6. FIGURA 5A: DISTRIBUCIÓN PEAK_HZ ------------------------------------
fprintf('>>> Generando Fig 5A: Distribución frecuencia masticatoria...\n');

figA = figure('Units','inches','Position',[1 1 5 5]);
axA  = axes(figA);  hold(axA,'on');

% Histograma
edges = linspace(floor(min(peak_hz)*10)/10, ceil(max(peak_hz)*10)/10, 12);
hh   = histogram(axA, peak_hz, edges, 'FaceColor',clrCas,'FaceAlpha',0.7,'EdgeColor','none');
ymax = max(hh.Values) * 1.15;   % máximo de conteos directamente del objeto

% Línea de mediana
xmed = double(median(peak_hz(:)));
yhi  = double(ymax) * 0.92;
plot(axA, [xmed xmed], [0 yhi], '--','Color',clrCas*0.7,'LineWidth',1.8);
text(axA, xmed+0.01, ymax*0.88, sprintf('Mediana = %.2f Hz', xmed), ...
     'FontSize',10,'Color',clrCas*0.7,'HorizontalAlignment','left');

% Estadísticos en gráfico
ylim(axA, [0 ymax]);
txt_d = sprintf('M = %.2f \\pm %.2f Hz\n[%.2f -- %.2f]', ...
        mean(peak_hz), std(peak_hz), min(peak_hz), max(peak_hz));
text(axA, 0.98, 0.98, txt_d, 'Units','normalized', ...
     'VerticalAlignment','top','HorizontalAlignment','right','FontSize',10);

set(axA,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.2,'GridColor',[.6 .6 .6]);
xlabel(axA,'Frecuencia masticatoria pico (Hz)');
ylabel(axA,'Número de sujetos');
title(axA,'Distribución f_{masticatoria} — Casos','FontWeight','bold','FontSize',13);

exportgraphics(figA, fullfile(dir_out,'Fig5A_ChewFreq_Distribution.png'),'Resolution',dpi);
close(figA);

%% -- 7. FIGURA 5B: PEAK_HZ vs ΔIES ----------------------------------------
fprintf('>>> Generando Fig 5B: peak_hz vs ΔIES...\n');

figB = figure('Units','inches','Position',[1 1 5 5]);
axB  = axes(figB);  hold(axB,'on');

scatter(axB, peak_hz, d_ies, 55, clrCas, ...
        'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');

% Regresión lineal (visual)
lm   = polyfit(peak_hz, d_ies, 1);
xrng = linspace(min(peak_hz)-0.05, max(peak_hz)+0.05, 100);
plot(axB, xrng, polyval(lm,xrng), '--','Color',[.4 .4 .4],'LineWidth',1.3);

% Línea de referencia y=0
yline(axB, 0, '-','Color',[.7 .7 .7],'LineWidth',1,'Alpha',0.5);

str = sprintf('\\rho = %.3f\np = %.4f  %s', rho_ies, p_ies, p2s(p_ies));
text(axB, 0.05, 0.95, str, 'Units','normalized', ...
     'VerticalAlignment','top','FontSize',11);

set(axB,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
xlabel(axB,'Frecuencia masticatoria pico (Hz)','FontSize',12);
ylabel(axB,'\DeltaIES (Ch − NoChew, ms)','FontSize',12);
title(axB,'f_{chew} vs \DeltaIES','FontWeight','bold','FontSize',13);

exportgraphics(figB, fullfile(dir_out,'Fig5B_ChewFreq_vs_DeltaIES.png'),'Resolution',dpi);
close(figB);

%% -- 8. FIGURA 5C: PEAK_HZ vs ΔZMI BETA EARLY -----------------------------
fprintf('>>> Generando Fig 5C: peak_hz vs zMI_Beta_Early...\n');

figC = figure('Units','inches','Position',[1 1 5 5]);
axC  = axes(figC);  hold(axC,'on');

scatter(axC, peak_hz, zmi_beta_early_ch, 55, clrCas, ...
        'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');

lm   = polyfit(peak_hz, zmi_beta_early_ch, 1);
xrng = linspace(min(peak_hz)-0.05, max(peak_hz)+0.05, 100);
plot(axC, xrng, polyval(lm,xrng), '--','Color',[.4 .4 .4],'LineWidth',1.3);

str = sprintf('\\rho = %.3f\np = %.4f  %s', rho_pac, p_pac, p2s(p_pac));
text(axC, 0.05, 0.95, str, 'Units','normalized', ...
     'VerticalAlignment','top','FontSize',11);

set(axC,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
xlabel(axC,'Frecuencia masticatoria pico (Hz)','FontSize',12);
ylabel(axC,'zMI \beta Early (Ch)','FontSize',12);
title(axC,'f_{chew} vs PAC \beta Early','FontWeight','bold','FontSize',13);

exportgraphics(figC, fullfile(dir_out,'Fig5C_ChewFreq_vs_PAC_Beta.png'),'Resolution',dpi);
close(figC);

fprintf('\n✓ S5_ChewFreq.m completado. Outputs en:\n  %s\n', dir_out);

%% ── FUNCIÓN AUXILIAR ──────────────────────────────────────────────────────

function s = p2s(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'NS';
    end
end
