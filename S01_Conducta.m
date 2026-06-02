%% Conducta


clear; close all; clc;

%% ── Rutas (desde S0_paths) ───────────────────────────────────────────────
P        = S0_paths();
FILE_TB  = P.file_beh;
DIR_OUT  = P.fig01;

%% ── Carga ────────────────────────────────────────────────────────────────
load(FILE_TB, 'tb_data_45');

cas_nc = tb_data_45.casos.nochew;
cas_ch = tb_data_45.casos.chew;
ctr_nc = tb_data_45.controles.nochew;
ctr_ch = tb_data_45.controles.chew;

% Vectores numéricos
CAS.rt_nc  = double(cas_nc.mean(:));
CAS.rt_ch  = double(cas_ch.mean(:));
CAS.acc_nc = double(cas_nc.acc(:));
CAS.acc_ch = double(cas_ch.acc(:));
CAS.ies_nc = double(cas_nc.ies(:));
CAS.ies_ch = double(cas_ch.ies(:));

CTR.rt_nc  = double(ctr_nc.mean(:));
CTR.rt_ch  = double(ctr_ch.mean(:));
CTR.acc_nc = double(ctr_nc.acc(:));
CTR.acc_ch = double(ctr_ch.acc(:));
CTR.ies_nc = double(ctr_nc.ies(:));
CTR.ies_ch = double(ctr_ch.ies(:));

nCas = numel(CAS.rt_nc);   % 30
nCtr = numel(CTR.rt_nc);   % 15
fprintf('Casos: n=%d | Controles: n=%d\n', nCas, nCtr);

%% ── Tests estadísticos ───────────────────────────────────────────────────
% Baseline (Nc): Mann-Whitney U
p_base_rt  = ranksum(CAS.rt_nc,  CTR.rt_nc);
p_base_ies = ranksum(CAS.ies_nc, CTR.ies_nc);

% Cambio Nc→Ch: Wilcoxon signed-rank
p_cas_rt  = signrank(CAS.rt_nc,  CAS.rt_ch);
p_cas_ies = signrank(CAS.ies_nc, CAS.ies_ch);
p_ctr_rt  = signrank(CTR.rt_nc,  CTR.rt_ch);
p_ctr_ies = signrank(CTR.ies_nc, CTR.ies_ch);

% Cohen's d (diferencias pareadas)
d_cas_ies = mean(CAS.ies_ch - CAS.ies_nc) / std(CAS.ies_ch - CAS.ies_nc);
d_ctr_ies = mean(CTR.ies_ch - CTR.ies_nc) / std(CTR.ies_ch - CTR.ies_nc);

%% ── Descomposición linealizada del ΔIES% ────────────────────────────────
%   IES = RT / acc
%   ΔIES% ≈ ΔRT%  +  (−Δacc%)      [1er orden, válido para cambios pequeños]
%   ΔRT%  = (RT_ch  − RT_nc)  / RT_nc  × 100
%   −Δacc%= −(acc_ch − acc_nc) / acc_nc × 100   (neg = mejora → reduce IES)

CAS.delta_ies_pct = (CAS.ies_ch - CAS.ies_nc) ./ CAS.ies_nc * 100;
CAS.rt_contrib    = (CAS.rt_ch  - CAS.rt_nc)  ./ CAS.rt_nc  * 100;
CAS.acc_contrib   = -(CAS.acc_ch - CAS.acc_nc) ./ CAS.acc_nc * 100;

CTR.delta_ies_pct = (CTR.ies_ch - CTR.ies_nc) ./ CTR.ies_nc * 100;
CTR.rt_contrib    = (CTR.rt_ch  - CTR.rt_nc)  ./ CTR.rt_nc  * 100;
CTR.acc_contrib   = -(CTR.acc_ch - CTR.acc_nc) ./ CTR.acc_nc * 100;

mu_cas_ies = mean(CAS.delta_ies_pct);
mu_ctr_ies = mean(CTR.delta_ies_pct);
mu_cas_rt  = mean(CAS.rt_contrib);
mu_ctr_rt  = mean(CTR.rt_contrib);
mu_cas_acc = mean(CAS.acc_contrib);
mu_ctr_acc = mean(CTR.acc_contrib);
se_cas_ies = std(CAS.delta_ies_pct) / sqrt(nCas);
se_ctr_ies = std(CTR.delta_ies_pct) / sqrt(nCtr);

%% ── Consola ──────────────────────────────────────────────────────────────
fprintf('\n════════════════════════════════════════════════════════\n');
fprintf(' RESUMEN — n=45 (30 Casos, 15 Controles)\n');
fprintf('════════════════════════════════════════════════════════\n');
fprintf('\nBASELINE (Nc):\n');
fprintf('  RT_mean: Casos=%.1f  Ctrl=%.1f  p=%.3f %s\n', ...
    mean(CAS.rt_nc), mean(CTR.rt_nc), p_base_rt, p2stars(p_base_rt));
fprintf('  IES:     Casos=%.1f  Ctrl=%.1f  p=%.3f %s\n', ...
    mean(CAS.ies_nc), mean(CTR.ies_nc), p_base_ies, p2stars(p_base_ies));
fprintf('\nCAMBIO Nc→Ch (Wilcoxon):\n');
fprintf('  Casos IES: Δ=%.1f ms (%.1f%%)  p=%.4f %s  d=%.2f\n', ...
    mean(CAS.ies_ch-CAS.ies_nc), mu_cas_ies, p_cas_ies, p2stars(p_cas_ies), d_cas_ies);
fprintf('  Ctrl  IES: Δ=%.1f ms (%.1f%%)  p=%.4f %s  d=%.2f\n', ...
    mean(CTR.ies_ch-CTR.ies_nc), mu_ctr_ies, p_ctr_ies, p2stars(p_ctr_ies), d_ctr_ies);
fprintf('\nDESCOMPOSICIÓN ΔIES (promedio grupal):\n');
fprintf('  Casos: RT=%.1f%%  ACC=%.1f%%  Total=%.1f%%\n', ...
    mu_cas_rt, mu_cas_acc, mu_cas_ies);
fprintf('  Ctrl:  RT=%.1f%%  ACC=%.1f%%  Total=%.1f%%\n', ...
    mu_ctr_rt, mu_ctr_acc, mu_ctr_ies);
fprintf('════════════════════════════════════════════════════════\n\n');

%% ════════════════════════════════════════════════
%  CONFIGURACIÓN GENERAL DE ESTÉTICA (Estilo Nature)
% ═════════════════════════════════════════════════
% Dimensiones idénticas para ambas figuras
figWidth  = 460;
figHeight = 520;

% Paleta de Colores: Verde (Casos) y Azul (Controles)
clrCas  = [0.15 0.55 0.30]; % Verde oscuro (RT Casos)
cCasACC = [0.60 0.85 0.70]; % Verde claro (ACC Casos)

clrCtr  = [0.18 0.42 0.78]; % Azul oscuro (RT Controles)
cCtrACC = [0.65 0.80 0.95]; % Azul claro (ACC Controles)

% ════════════════════════════════════════════════
%  FIGURA A — Baseline RT (Nc): Casos vs Controles
% ═════════════════════════════════════════════════
figA = figure('Name','Baseline RT','Color','w','Position',[60 80 figWidth figHeight]);
axA  = axes('Parent',figA);
hold(axA,'on');

xCas = 1;  xCtr = 2;
rng(42);
jCas = (rand(nCas,1)-0.5)*0.22;
jCtr = (rand(nCtr,1)-0.5)*0.22;

% Puntos individuales
scatter(axA, xCas+jCas, CAS.rt_nc, 30, clrCas, ...
    'filled','MarkerFaceAlpha',0.30,'MarkerEdgeColor','none');
scatter(axA, xCtr+jCtr, CTR.rt_nc, 30, clrCtr, ...
    'filled','MarkerFaceAlpha',0.30,'MarkerEdgeColor','none');

% Cajas (IQR) y Bigotes
boxdata = {CAS.rt_nc, xCas, clrCas; CTR.rt_nc, xCtr, clrCtr};
for k = 1:2
    d  = boxdata{k,1};
    xp = boxdata{k,2};
    c  = boxdata{k,3};
    
    q1 = prctile(d,25); q3 = prctile(d,75); med = median(d);
    iq = q3-q1;
    wlo = max(d(d >= q1-1.5*iq));
    whi = min(d(d <= q3+1.5*iq));
    
    % Caja principal
    fill(axA, xp+[-0.18 0.18 0.18 -0.18 -0.18], ...
              [q1 q1 q3 q3 q1], c, ...
        'FaceAlpha',0.20,'EdgeColor',c,'LineWidth',1.2);
        
    % Línea de la mediana
    plot(axA, xp+[-0.18 0.18],[med med],'-','Color',c,'LineWidth',2.5);
    
    % Líneas verticales de los bigotes
    plot(axA, [xp xp],[whi q3],'-','Color',c,'LineWidth',1.2);
    plot(axA, [xp xp],[q1 wlo],'-','Color',c,'LineWidth',1.2);
    
    % Topes horizontales de los bigotes (caps)
    cw = 0.05; % Ancho del tope
    plot(axA, [xp-cw xp+cw],[whi whi],'-','Color',c,'LineWidth',1.2);
    plot(axA, [xp-cw xp+cw],[wlo wlo],'-','Color',c,'LineWidth',1.2);
end

% Barra de significancia
ymax = max([CAS.rt_nc; CTR.rt_nc]);
ysig = ymax * 1.05;
plot(axA,[xCas xCtr],[ysig ysig],'k-','LineWidth',1.0);
plot(axA,[xCas xCas],[ysig-10 ysig],'k-','LineWidth',1.0);
plot(axA,[xCtr xCtr],[ysig-10 ysig],'k-','LineWidth',1.0);
text(axA, 1.5, ysig+8, 'n.s.','HorizontalAlignment','center','FontSize',10,'FontName','Arial');

% Etiquetas inferiores
text(axA, xCas, 340, sprintf('%.0f ± %.0f ms', ...
    mean(CAS.rt_nc), std(CAS.rt_nc)/sqrt(nCas)), ...
    'HorizontalAlignment','center','FontSize',9,'Color',clrCas,'FontName','Arial');
text(axA, xCtr, 340, sprintf('%.0f ± %.0f ms', ...
    mean(CTR.rt_nc), std(CTR.rt_nc)/sqrt(nCtr)), ...
    'HorizontalAlignment','center','FontSize',9,'Color',clrCtr,'FontName','Arial');

% Formato Estilo Nature
axA.XTick      = [xCas xCtr];
axA.XTickLabel = {'Cases','Controls'};
axA.XLim       = [0.5 2.5];
axA.YLim       = [300 ysig+25];
axA.FontName   = 'Arial';
axA.FontSize   = 10;
axA.Box        = 'off';
axA.TickDir    = 'out';
axA.LineWidth  = 1.0;
ylabel(axA,'Baseline RT (ms)','FontSize',11,'FontWeight','bold','FontName','Arial');
hold(axA,'off');

fA = fullfile(DIR_OUT,'Fig_Baseline_RT.png');
exportgraphics(figA, fA, 'Resolution',300);
fprintf('[OK] Fig A → %s\n', fA);

%% ════════════════════════════════════════════
%  FIGURA B — IES Shift Decomposition (Nc → Ch)
% ═════════════════════════════════════════════
figB = figure('Name','IES Decomposition','Color','w','Position',[560 80 figWidth figHeight]);
axB  = axes('Parent',figB);
hold(axB,'on');

xCas = 1;    xCtr = 2.4;   bw = 0.5;

% ── Función auxiliar para dibujar barras con bordes ─────────────────────
draw_bar = @(x, y0, y1, width, color) fill(axB, ...
    [x-width/2 x+width/2 x+width/2 x-width/2 x-width/2], ...
    [y0 y0 y1 y1 y0], color, 'FaceAlpha', 0.9, 'EdgeColor', 'k', 'LineWidth', 1);

% ── Barras apiladas ──────────────────────────────────────────────────────
% Casos: RT y ACC
draw_bar(xCas, 0, mu_cas_rt, bw, clrCas);
draw_bar(xCas, mu_cas_rt, mu_cas_rt + mu_cas_acc, bw, cCasACC);

% Controles: RT y ACC
draw_bar(xCtr, 0, mu_ctr_rt, bw, clrCtr);
draw_bar(xCtr, mu_ctr_rt, mu_ctr_rt + mu_ctr_acc, bw, cCtrACC);

% ── Línea de referencia (Cero) ───────────────────────────────────────────
yline(axB, 0, 'k-', 'LineWidth', 1.0);

% ── Barras de Error (Centradas) y Texto de Total ─────────────────────────
% Casos Total
y_cas_tot = mu_cas_rt + mu_cas_acc; 
plot(axB, [xCas xCas], [y_cas_tot, y_cas_tot - se_cas_ies], 'k-', 'LineWidth', 1.5); % Línea vertical
plot(axB, [xCas-0.08 xCas+0.08], [y_cas_tot - se_cas_ies, y_cas_tot - se_cas_ies], 'k-', 'LineWidth', 1.5); % Tope inferior

% Texto Total Casos (Ubicado debajo del bigote de error)
text(axB, xCas, (y_cas_tot - se_cas_ies) - 0.8, sprintf('Total:\n%.1f%%', y_cas_tot), ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',10,'FontName','Arial','Color','k');

% Controles Total
y_ctr_tot = mu_ctr_rt + mu_ctr_acc; 
plot(axB, [xCtr xCtr], [y_ctr_tot, y_ctr_tot - se_ctr_ies], 'k-', 'LineWidth', 1.5); % Línea vertical
plot(axB, [xCtr-0.08 xCtr+0.08], [y_ctr_tot - se_ctr_ies, y_ctr_tot - se_ctr_ies], 'k-', 'LineWidth', 1.5); % Tope inferior

% Texto Total Controles (Ubicado debajo del bigote de error)
text(axB, xCtr, (y_ctr_tot - se_ctr_ies) - 0.8, sprintf('Total:\n%.1f%%', y_ctr_tot), ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',10,'FontName','Arial','Color','k');

% ── Etiquetas numéricas dentro de las barras ─────────────────────────────
if abs(mu_cas_rt) > 1.0
    text(axB, xCas, mu_cas_rt/2, sprintf('%.1f%%', mu_cas_rt), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',10,'Color','w','FontName','Arial');
end
if abs(mu_cas_acc) > 1.0
    text(axB, xCas, mu_cas_rt + mu_cas_acc/2, sprintf('%.1f%%', mu_cas_acc), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',10,'Color','k','FontName','Arial');
end
if abs(mu_ctr_rt) > 1.0
    text(axB, xCtr, mu_ctr_rt/2, sprintf('%.1f%%', mu_ctr_rt), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',10,'Color','w','FontName','Arial');
end
if abs(mu_ctr_acc) > 1.0
    text(axB, xCtr, mu_ctr_rt + mu_ctr_acc/2, sprintf('%.1f%%', mu_ctr_acc), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',10,'Color','k','FontName','Arial');
end

% ── Anotación de Significancia ───────────────────────────────────────────
text(axB, xCas, 2.5, p2stars(p_cas_ies), ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',12,'FontName','Arial','Color','k');
text(axB, xCtr, 2.5, p2stars(p_ctr_ies), ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',12,'FontName','Arial','Color','k');

% ── Leyenda (Diamante eliminado) ─────────────────────────────────────────
h1 = fill(axB, NaN, NaN, clrCas,  'EdgeColor', 'k');
h2 = fill(axB, NaN, NaN, cCasACC, 'EdgeColor', 'k');
h3 = fill(axB, NaN, NaN, clrCtr,  'EdgeColor', 'k');
h4 = fill(axB, NaN, NaN, cCtrACC, 'EdgeColor', 'k');

legend(axB, [h1, h2, h3, h4], ...
    {'RT (Cases)', 'ACC (Cases)', ...
     'RT (Controls)', 'ACC (Controls)'}, ...
    'Location','southwest', 'FontSize', 9, 'FontName', 'Arial', 'Box', 'off');

% ── Formato Estilo Nature ────────────────────────────────────────────────
axB.XTick      = [xCas xCtr];
axB.XTickLabel = {'Cases', 'Controls'};
axB.XLim       = [0.3 3.3];
axB.YLim       = [-26 4]; 
axB.TickDir    = 'out';
axB.FontName   = 'Arial';
axB.FontSize   = 10;
axB.Box        = 'off';
axB.LineWidth  = 1.0;
ylabel(axB,'\Delta IES (%) (Chew - NoChew)','FontSize',11,'FontWeight','bold','FontName','Arial');
hold(axB,'off');

fB = fullfile(DIR_OUT,'Fig_IES_Shift_Decomposition.png');
exportgraphics(figB, fB, 'Resolution', 300);
fprintf('[OK] Fig B → %s\n', fB);

%% ════════════════════════════════════════════════════════════════════════
%  REPORTE TXT PARA EL MANUSCRITO (Conducta y Baseline)
% ════════════════════════════════════════════════════════════════════════
f_rep_cond = fullfile(DIR_OUT, 'Reporte_Conducta.txt');
fid = fopen(f_rep_cond, 'w');

% Recalcular tests para extraer los estadísticos exactos (Z o W) requeridos por las revistas
[p_base_ies, ~, stat_base_ies] = ranksum(CAS.ies_nc, CTR.ies_nc);
[p_cas_ies, ~, stat_cas_ies]  = signrank(CAS.ies_nc, CAS.ies_ch);
[p_ctr_ies, ~, stat_ctr_ies]  = signrank(CTR.ies_nc, CTR.ies_ch);

% Calcular Desviaciones Estándar (SD) para el texto
sd_cas_nc = std(CAS.ies_nc); sd_cas_ch = std(CAS.ies_ch);
sd_ctr_nc = std(CTR.ies_nc); sd_ctr_ch = std(CTR.ies_ch);

fprintf(fid, '================================================================\n');
fprintf(fid, '  CONDUCTA Y BASELINE — Resultados paper (S1_Conducta.m)\n');
fprintf(fid, '  Generado: %s\n', datestr(now));
fprintf(fid, '  N = %d Casos | %d Controles\n', nCas, nCtr);
fprintf(fid, '================================================================\n\n');

fprintf(fid, '── 1. EQUIVALENCIA BASAL (No-Chew / Block 1) ───────────────────\n');
fprintf(fid, '  IES Casos (Nc):     %.1f ± %.1f ms\n', mean(CAS.ies_nc), sd_cas_nc);
fprintf(fid, '  IES Controles (Nc): %.1f ± %.1f ms\n', mean(CTR.ies_nc), sd_ctr_nc);
if isfield(stat_base_ies, 'zval')
    fprintf(fid, '  Mann-Whitney U: p = %.4f (Z = %.3f)\n\n', p_base_ies, stat_base_ies.zval);
else
    fprintf(fid, '  Mann-Whitney U: p = %.4f\n\n', p_base_ies);
end

fprintf(fid, '── 2. EFECTO DE LA MASTICACIÓN (Wilcoxon Signed-Rank) ──────────\n');
fprintf(fid, '  [GRUPO CASOS]\n');
fprintf(fid, '  IES Nc (Reposo): %.1f ± %.1f ms\n', mean(CAS.ies_nc), sd_cas_nc);
fprintf(fid, '  IES Ch (Mastic): %.1f ± %.1f ms\n', mean(CAS.ies_ch), sd_cas_ch);
if isfield(stat_cas_ies, 'zval')
    fprintf(fid, '  Wilcoxon: p = %.4f (Z = %.3f)\n', p_cas_ies, stat_cas_ies.zval);
else
    fprintf(fid, '  Wilcoxon: p = %.4f (W = %.1f)\n', p_cas_ies, stat_cas_ies.signedrank);
end
fprintf(fid, '  Cohen''s d (pareado): %.3f\n\n', d_cas_ies);

fprintf(fid, '  [GRUPO CONTROL]\n');
fprintf(fid, '  IES Block 1: %.1f ± %.1f ms\n', mean(CTR.ies_nc), sd_ctr_nc);
fprintf(fid, '  IES Block 2: %.1f ± %.1f ms\n', mean(CTR.ies_ch), sd_ctr_ch);
if isfield(stat_ctr_ies, 'zval')
    fprintf(fid, '  Wilcoxon: p = %.4f (Z = %.3f)\n\n', p_ctr_ies, stat_ctr_ies.zval);
else
    fprintf(fid, '  Wilcoxon: p = %.4f (W = %.1f)\n\n', p_ctr_ies, stat_ctr_ies.signedrank);
end

fprintf(fid, '── 3. REDUCCIÓN PORCENTUAL (Descomposición ΔIES) ───────────────\n');
fprintf(fid, '  Casos Reducción IES:   %.1f%% (RT contrib: %.1f%%, ACC contrib: %.1f%%)\n', ...
    mean(CAS.delta_ies_pct), mean(CAS.rt_contrib), mean(CAS.acc_contrib));
fprintf(fid, '  Control Reducción IES: %.1f%% (RT contrib: %.1f%%, ACC contrib: %.1f%%)\n', ...
    mean(CTR.delta_ies_pct), mean(CTR.rt_contrib), mean(CTR.acc_contrib));

fclose(fid);
fprintf('\n>>> [ÉXITO] Reporte de conducta guardado en: %s\n', f_rep_cond);
%% ════════════════════════════════════════════════════════════════════════
%  FUNCIONES LOCALES
% ════════════════════════════════════════════════════════════════════════

function patch_bar(ax, xpos, ybot, ytop, bw, clr)
%PATCH_BAR  Dibuja un rectángulo (barra) entre ybot e ytop.
    hw = bw/2;
    xs = [xpos-hw, xpos+hw, xpos+hw, xpos-hw, xpos-hw];
    ys = [ybot,    ybot,    ytop,    ytop,    ybot   ];
    fill(ax, xs, ys, clr, 'EdgeColor','none');
end

function s = p2stars(p)
%P2STARS  Devuelve etiqueta de significancia.
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'n.s';
    end
end
