%% ============================================================================
%  S3b_FOOOF_PeriodicSpectra.m  —  Espectros FOOOF (3 figuras)
%  ----------------------------------------------------------------------------
%  Carga FOOOF_Workspace.mat y genera:
%
%    Fig 1. PSD Broadband (pre-FOOOF) con fit aperiódico superpuesto
%           → muestra el aplanamiento de pendiente en Casos durante Chew
%           → REQUIERE que S4 haya guardado GR.*.PSD_Nc/Ch y GR.*.AP_Nc/Ch
%
%    Fig 2. Delta Periódico (Ch - Nc) con estadística de interacción
%
%    Fig 3. Espectros Periódicos por Bloque (No Chew | Chew)
%  ============================================================================
clear; clc; close all;

%% -- 1. RUTAS Y CARGA DEL WORKSPACE ------------------------------------------
P       = S0_paths();
dir_out = P.fig03;
f_ws    = P.file_fooof_ws;

if ~exist(dir_out, 'dir'); mkdir(dir_out); end

fprintf('>>> Cargando FOOOF_Workspace.mat...\n');
W  = load(f_ws);
GR = W.GR;
f  = GR.Cases.f(:);

%% -- 2. ESTÉTICA CONSISTENTE -------------------------------------------------
% Paleta exacta del paper (tomada de S1_IES / S2_behavior)
c_case      = [0.15 0.55 0.30];   % Verde oscuro Cases  (líneas)
c_case_fill = [0.60 0.85 0.70];   % Verde claro  Cases  (sombras SE)
c_ctrl      = [0.18 0.42 0.78];   % Azul oscuro  Controls (líneas)
c_ctrl_fill = [0.65 0.80 0.95];   % Azul claro   Controls (sombras SE)

% Dimensiones estándar del paper
figW  = 460;
figH  = 520;

% Condición se codifica por estilo de línea, NO por color
% (color = grupo, línea = condición, grosor = grupo)

set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 11, ...
       'DefaultFigureColor', 'w');
DPI = 300;

bnd    = {[4 7], [8 12], [13 30]};
bnd_nm = {'\theta', '\alpha', '\beta'};

%% ============================================================================
%  FIGURA 1: PSD BROADBAND (pre-FOOOF) + componente aperiódico
%  Requiere: GR.*.PSD_Nc, GR.*.PSD_Ch, GR.*.AP_Nc, GR.*.AP_Ch
%  Si aún no están en el workspace, se salta con advertencia.
% ============================================================================
has_raw = isfield(GR.Cases, 'PSD_Nc') && isfield(GR.Cases, 'PSD_Ch');

if ~has_raw
    warning(['Fig 1 OMITIDA: GR no contiene PSD_Nc/PSD_Ch.\n' ...
             'Agregá esos campos en S4_FOOOF y regenerá el workspace.']);
else
    fprintf('>>> Generando Figura 1: PSD Broadband + Pendiente Aperiódica...\n');

    % --- Promedios y SE (real() elimina parte imaginaria residual Python→MATLAB) ---
    psd_cn = real(GR.Cases.PSD_Nc);    psd_cc = real(GR.Cases.PSD_Ch);
    psd_tn = real(GR.Controls.PSD_Nc); psd_tc = real(GR.Controls.PSD_Ch);

    m_cn = nanmean(psd_cn, 2);  se_cn = nanstd(psd_cn, 0, 2) / sqrt(GR.Cases.n);
    m_cc = nanmean(psd_cc, 2);  se_cc = nanstd(psd_cc, 0, 2) / sqrt(GR.Cases.n);
    m_tn = nanmean(psd_tn, 2);  se_tn = nanstd(psd_tn, 0, 2) / sqrt(GR.Controls.n);
    m_tc = nanmean(psd_tc, 2);  se_tc = nanstd(psd_tc, 0, 2) / sqrt(GR.Controls.n);

    % AP fit (opcional — se dibuja si existe)
    has_ap = isfield(GR.Cases, 'AP_Nc');
    if has_ap
        m_ap_cn = nanmean(real(GR.Cases.AP_Nc),    2);
        m_ap_cc = nanmean(real(GR.Cases.AP_Ch),    2);
        m_ap_tn = nanmean(real(GR.Controls.AP_Nc), 2);
        m_ap_tc = nanmean(real(GR.Controls.AP_Ch), 2);
    end

    % Límites Y automáticos con margen
    all_m   = [m_cn; m_cc; m_tn; m_tc];
    ylim_lo = floor(min(all_m) - 1);
    ylim_hi = ceil(max(all_m)  + 1);
    ylims_raw = [ylim_lo ylim_hi];

    fig1 = figure('Name', 'FOOOF_PSD_Raw', 'Position', [80 80 figW figH], 'Visible', 'off');
    ax1  = axes(fig1); hold(ax1, 'on');

    % Fondos de banda (solo para referencia visual)
    for bi = 1:3
        patch(ax1, [bnd{bi}(1) bnd{bi}(2) bnd{bi}(2) bnd{bi}(1)], ...
              [ylims_raw(1) ylims_raw(1) ylims_raw(2) ylims_raw(2)], ...
              [0.94 0.94 0.94], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        text(ax1, mean(bnd{bi}), ylims_raw(2) - 0.6, bnd_nm{bi}, ...
             'HorizontalAlignment', 'center', 'FontSize', 12, ...
             'Color', [.55 .55 .55], 'Interpreter', 'tex');
    end

    % --- Sombras SE (fill claro por grupo, sin alpha tricks) ---
    fill(ax1, [f; flipud(f)], [m_tn+se_tn; flipud(m_tn-se_tn)], ...
         c_ctrl_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax1, [f; flipud(f)], [m_tc+se_tc; flipud(m_tc-se_tc)], ...
         c_ctrl_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax1, [f; flipud(f)], [m_cn+se_cn; flipud(m_cn-se_cn)], ...
         c_case_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax1, [f; flipud(f)], [m_cc+se_cc; flipud(m_cc-se_cc)], ...
         c_case_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % --- Líneas PSD raw ---
    % Color = grupo  |  Estilo = condición  |  Grosor = grupo
    plot(ax1, f, m_tn, '-',  'Color', c_ctrl, 'LineWidth', 1.4, 'DisplayName', 'Controls – No Chew');
    plot(ax1, f, m_tc, '--', 'Color', c_ctrl, 'LineWidth', 1.4, 'DisplayName', 'Controls – Chew');
    plot(ax1, f, m_cn, '-',  'Color', c_case, 'LineWidth', 2.4, 'DisplayName', 'Cases – No Chew');
    plot(ax1, f, m_cc, '--', 'Color', c_case, 'LineWidth', 2.4, 'DisplayName', 'Cases – Chew');

    % --- Fit aperiódico: mismo color, línea punteada, muy fina ---
    if has_ap
        plot(ax1, f, m_ap_tn, ':', 'Color', [c_ctrl 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot(ax1, f, m_ap_tc, ':', 'Color', [c_ctrl 0.55], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot(ax1, f, m_ap_cn, ':', 'Color', [c_case 0.55], 'LineWidth', 1.4, 'HandleVisibility', 'off');
        plot(ax1, f, m_ap_cc, ':', 'Color', [c_case 0.55], 'LineWidth', 1.4, 'HandleVisibility', 'off');
        text(ax1, 0.99, 0.03, 'Dotted: aperiodic fit', 'Units', 'normalized', ...
             'HorizontalAlignment', 'right', 'FontSize', 8, 'Color', [.5 .5 .5]);
    end

    % --- Formato ---
    set(ax1, 'XScale', 'log', 'XLim', [3 35], 'YLim', ylims_raw, ...
        'TickDir', 'out', 'Box', 'off', 'XGrid', 'on', 'YGrid', 'on', ...
        'GridAlpha', 0.12, 'Layer', 'top');
    set(ax1, 'XTick', [4 8 13 20 30], 'XTickLabel', {'4','8','13','20','30'});
    xlabel(ax1, 'Frequency (Hz)', 'FontWeight', 'bold');
    ylabel(ax1, 'Power (dB)', 'FontWeight', 'bold');
    title(ax1, 'Broadband PSD (pre-FOOOF) — aperiodic slope', 'FontSize', 13);
    legend(ax1, 'show', 'Location', 'southwest', 'Box', 'off', 'FontSize', 9);

    exportgraphics(fig1, fullfile(dir_out, 'Fig3_FOOOF_PSD_Raw.png'), 'Resolution', DPI);
    close(fig1);
    fprintf('    → Fig3_FOOOF_PSD_Raw.png guardada.\n');
end

%% ============================================================================
%  FIGURA 2: DELTA PERIÓDICO (Ch - Nc)
% ============================================================================
fprintf('>>> Generando Figura 2: Delta Periódico (Ch - Nc)...\n');
fig2 = figure('Name', 'FOOOF_Delta', 'Position', [80 80 figW figH], 'Visible', 'off');
ax2  = axes(fig2); hold(ax2, 'on');

ylims_delta = [-1.5 2.6];
set(ax2, 'YLim', ylims_delta, 'XLim', [3 35]);

for bi = 1:3
    patch(ax2, [bnd{bi}(1) bnd{bi}(2) bnd{bi}(2) bnd{bi}(1)], ...
          [ylims_delta(1) ylims_delta(1) ylims_delta(2) ylims_delta(2)], ...
          [0.94 0.94 0.94], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(ax2, mean(bnd{bi}), ylims_delta(2)-0.12, bnd_nm{bi}, ...
         'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [.55 .55 .55], ...
         'Interpreter', 'tex');
end

d_cas = GR.Cases.Res_Ch    - GR.Cases.Res_Nc;
d_ctr = GR.Controls.Res_Ch - GR.Controls.Res_Nc;
m_cas = nanmean(d_cas, 2);  se_cas = nanstd(d_cas, 0, 2) / sqrt(GR.Cases.n);
m_ctr = nanmean(d_ctr, 2);  se_ctr = nanstd(d_ctr, 0, 2) / sqrt(GR.Controls.n);

fill(ax2, [f; flipud(f)], [m_ctr+se_ctr; flipud(m_ctr-se_ctr)], ...
     c_ctrl_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax2, f, m_ctr, 'Color', c_ctrl, 'LineWidth', 1.5, 'DisplayName', 'Controls');
fill(ax2, [f; flipud(f)], [m_cas+se_cas; flipud(m_cas-se_cas)], ...
     c_case_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax2, f, m_cas, 'Color', c_case, 'LineWidth', 2.5, 'DisplayName', 'Cases');

y_rows  = [2.22 1.94 1.66];
row_col = {[0 0 0], c_case, c_ctrl};
row_lbl = {'Inter', 'Cases', 'Ctrls'};

for bi = 1:2
    idx  = f >= bnd{bi}(1) & f <= bnd{bi}(2);
    xpos = mean(bnd{bi});
    d1   = nanmean(d_cas(idx,:), 1);
    d2   = nanmean(d_ctr(idx,:), 1);
    [~,p1]=ttest2(d1, d2); [~,p2]=ttest(d1); [~,p3]=ttest(d2);
    pvec = [p1 p2 p3];
    for ri = 1:3
        if pvec(ri) < 0.001; sym2='***'; fs2=14; fw2='bold';
        elseif pvec(ri) < 0.01; sym2='**'; fs2=14; fw2='bold';
        elseif pvec(ri) < 0.05; sym2='*';  fs2=14; fw2='bold';
        else; sym2='n.s.'; fs2=8; fw2='normal'; end
        text(ax2, xpos, y_rows(ri), sym2, 'HorizontalAlignment', 'center', ...
             'FontSize', fs2, 'FontWeight', fw2, 'Color', row_col{ri});
    end
end
for ri = 1:3
    text(ax2, 3.15, y_rows(ri), row_lbl{ri}, 'FontSize', 8, ...
         'Color', row_col{ri}, 'VerticalAlignment', 'middle');
end

yline(ax2, 0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
set(ax2, 'TickDir', 'out', 'Box', 'off', 'XGrid', 'on', 'YGrid', 'on', ...
    'GridAlpha', 0.12, 'Layer', 'top');
xlabel(ax2, 'Frequency (Hz)', 'FontWeight', 'bold');
ylabel(ax2, '\Delta Periodic Power Ch-Nc (dB)', 'FontWeight', 'bold');
title(ax2, 'Chewing-induced oscillatory change (FOOOF residual)', 'FontSize', 13);
legend(ax2, 'show', 'Location', 'northeast', 'Box', 'off');

exportgraphics(fig2, fullfile(dir_out, 'Fig3_FOOOF_DeltaPeriodic.png'), 'Resolution', DPI);
close(fig2);

%% ============================================================================
%  FIGURA 3: ESPECTROS PERIÓDICOS POR BLOQUE (Nc vs Ch)
% ============================================================================
fprintf('>>> Generando Figura 3: Espectros por Bloque (Nc vs Ch)...\n');
fig3 = figure('Name', 'FOOOF_Blocks', 'Position', [100 100 figW*2 figH], 'Visible', 'off');
ylims_block = [-0.2 7];
cond_names  = {'No Chew (Baseline)', 'Chew (Active)'};

for subplot_idx = 1:2
    ax3 = subplot(1, 2, subplot_idx); hold(ax3, 'on');

    set(ax3, 'YLim', ylims_block, 'XLim', [3 35]);
    for bi = 1:3
        patch(ax3, [bnd{bi}(1) bnd{bi}(2) bnd{bi}(2) bnd{bi}(1)], ...
              [ylims_block(1) ylims_block(1) ylims_block(2) ylims_block(2)], ...
              [0.94 0.94 0.94], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        text(ax3, mean(bnd{bi}), ylims_block(2)-0.3, bnd_nm{bi}, ...
             'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [.55 .55 .55], ...
             'Interpreter', 'tex');
    end

    if subplot_idx == 1
        d_cas = GR.Cases.Res_Nc;   d_ctr = GR.Controls.Res_Nc;
    else
        d_cas = GR.Cases.Res_Ch;   d_ctr = GR.Controls.Res_Ch;
    end

    m_cas = nanmean(d_cas, 2);  se_cas = nanstd(d_cas, 0, 2) / sqrt(GR.Cases.n);
    m_ctr = nanmean(d_ctr, 2);  se_ctr = nanstd(d_ctr, 0, 2) / sqrt(GR.Controls.n);

    fill(ax3, [f; flipud(f)], [m_ctr+se_ctr; flipud(m_ctr-se_ctr)], ...
         c_ctrl_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax3, f, m_ctr, 'Color', c_ctrl, 'LineWidth', 1.5, 'DisplayName', 'Controls');
    fill(ax3, [f; flipud(f)], [m_cas+se_cas; flipud(m_cas-se_cas)], ...
         c_case_fill, 'FaceAlpha', 0.50, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax3, f, m_cas, 'Color', c_case, 'LineWidth', 2.5, 'DisplayName', 'Cases');

    y_row_sym = ylims_block(2) - 0.7;
    y_row_txt = ylims_block(2) - 1.1;

    for bi = 1:2
        idx  = f >= bnd{bi}(1) & f <= bnd{bi}(2);
        xpos = mean(bnd{bi});
        d1   = nanmean(d_cas(idx,:), 1);
        d2   = nanmean(d_ctr(idx,:), 1);
        [~, p] = ttest2(d1, d2);
        if p < 0.001; sym='***'; fs=14; fw='bold'; col=[0 0 0];
        elseif p < 0.01; sym='**'; fs=14; fw='bold'; col=[0 0 0];
        elseif p < 0.05; sym='*';  fs=14; fw='bold'; col=[0 0 0];
        else; sym='n.s.'; fs=10; fw='normal'; col=[.5 .5 .5]; end
        text(ax3, xpos, y_row_sym, sym, 'HorizontalAlignment', 'center', ...
             'FontSize', fs, 'FontWeight', fw, 'Color', col);
    end
    text(ax3, 3.5, y_row_txt, 'Cases vs Ctrls', 'FontSize', 9, ...
         'Color', [.4 .4 .4], 'VerticalAlignment', 'middle');

    yline(ax3, 0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    set(ax3, 'TickDir', 'out', 'Box', 'off', 'XGrid', 'on', 'YGrid', 'on', ...
        'GridAlpha', 0.12, 'Layer', 'top');
    xlabel(ax3, 'Frequency (Hz)', 'FontWeight', 'bold');
    if subplot_idx == 1
        ylabel(ax3, 'Periodic Power (FOOOF residual)', 'FontWeight', 'bold');
    end
    title(ax3, cond_names{subplot_idx}, 'FontSize', 13);
    if subplot_idx == 2
        legend(ax3, 'show', 'Location', 'northeast', 'Box', 'off');
    end
end

exportgraphics(fig3, fullfile(dir_out, 'Fig3_FOOOF_Blocks.png'), 'Resolution', DPI);
close(fig3);

%% ============================================================================
%  FIGURA 4: EXPONENT BOXPLOT — Cases NoChew vs Chew (paired)
% ============================================================================
fprintf('>>> Generando Figura 4: Aperiodic Exponent boxplot (Cases)...\n');

exp_cas_nc = GR.Cases.exp_Nc(:);
exp_cas_ch = GR.Cases.exp_Ch(:);
nCas       = GR.Cases.n;

% Colores: gris para baseline (NoChew), verde oscuro para Chew
clrNc  = [0.60 0.60 0.60];   % gris — condición baseline
clrCas = [0.72 0.15 0.15];   % rojo oscuro — condición activa (Chew, Casos)

[p_cas, ~] = signrank(exp_cas_ch, exp_cas_nc);

fig4 = figure('Name','FOOOF_Exponent_Box','Position',[80 80 figW figH],'Visible','off');
ax4  = axes(fig4); hold(ax4,'on');

xNc = 1; xCh = 2;
bw  = 0.35;

draw_box(ax4, exp_cas_nc, xNc, bw, clrNc,  0.55);
draw_box(ax4, exp_cas_ch, xCh, bw, clrCas, 0.55);

% Líneas pareadas individuales
rng(42);
for i = 1:nCas
    line(ax4, [xNc xCh], [exp_cas_nc(i) exp_cas_ch(i)], ...
         'Color', [.55 .55 .55 .25], 'LineWidth', 0.7);
end

% Scatter individual
jit = 0.10;
scatter(ax4, xNc + randn(nCas,1)*jit, exp_cas_nc, 30, clrNc,  ...
        'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
scatter(ax4, xCh + randn(nCas,1)*jit, exp_cas_ch, 30, clrCas, ...
        'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');

% Sin errorbar de media — estructura idéntica al Panel A (Baseline Comparison)
sem_nc = std(exp_cas_nc)/sqrt(nCas);
sem_ch = std(exp_cas_ch)/sqrt(nCas);

% Bracket de significancia
ymax = max([exp_cas_nc; exp_cas_ch]) * 1.02;
yb   = ymax + 0.04;
line(ax4, [xNc xNc xCh xCh], [yb-0.02 yb yb yb-0.02], 'Color','k','LineWidth',1.2);
text(ax4, mean([xNc xCh]), yb+0.025, p2s(p_cas), ...
     'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

set(ax4, 'XTick',[xNc xCh], 'XTickLabel',{'No Chew','Chew'}, ...
         'XLim',[0.5 2.5], 'YGrid','on', 'GridAlpha',0.20, ...
         'GridColor',[.6 .6 .6], 'Box','off', 'TickDir','out', 'FontSize',11);
ylabel(ax4, 'Aperiodic Exponent', 'FontWeight','bold');

exportgraphics(fig4, fullfile(dir_out,'Fig3_FOOOF_Exponent_Box.png'), 'Resolution', DPI);
close(fig4);

%% ============================================================================
%  FIGURA 5: PANEL COMPUESTO — PSD Broadband + Exponent Boxplot
%  → listo para Figure 2A del paper
% ============================================================================
if has_raw
    fprintf('>>> Generando Figura 5: Panel compuesto PSD + Exponent...\n');

    fig5 = figure('Name','FOOOF_PanelA','Position',[80 80 figW*2 figH],'Visible','off');

    % ── Panel izquierdo: PSD broadband (reutiliza datos de Fig 1) ──────────
    ax5L = subplot(1,2,1); hold(ax5L,'on');

    psd_cn = real(GR.Cases.PSD_Nc);    psd_cc = real(GR.Cases.PSD_Ch);
    psd_tn = real(GR.Controls.PSD_Nc); psd_tc = real(GR.Controls.PSD_Ch);
    m_cn = nanmean(psd_cn,2); se_cn = nanstd(psd_cn,0,2)/sqrt(GR.Cases.n);
    m_cc = nanmean(psd_cc,2); se_cc = nanstd(psd_cc,0,2)/sqrt(GR.Cases.n);
    m_tn = nanmean(psd_tn,2); se_tn = nanstd(psd_tn,0,2)/sqrt(GR.Controls.n);
    m_tc = nanmean(psd_tc,2); se_tc = nanstd(psd_tc,0,2)/sqrt(GR.Controls.n);

    all_m     = [m_cn; m_cc; m_tn; m_tc];
    ylims_raw = [floor(min(all_m)-1) ceil(max(all_m)+1)];

    for bi = 1:3
        patch(ax5L,[bnd{bi}(1) bnd{bi}(2) bnd{bi}(2) bnd{bi}(1)], ...
              [ylims_raw(1) ylims_raw(1) ylims_raw(2) ylims_raw(2)], ...
              [0.94 0.94 0.94],'EdgeColor','none','HandleVisibility','off');
        text(ax5L,mean(bnd{bi}),ylims_raw(2)-0.6,bnd_nm{bi}, ...
             'HorizontalAlignment','center','FontSize',11, ...
             'Color',[.55 .55 .55],'Interpreter','tex');
    end

    fill(ax5L,[f;flipud(f)],[m_tn+se_tn;flipud(m_tn-se_tn)],c_ctrl_fill,'FaceAlpha',0.50,'EdgeColor','none','HandleVisibility','off');
    fill(ax5L,[f;flipud(f)],[m_tc+se_tc;flipud(m_tc-se_tc)],c_ctrl_fill,'FaceAlpha',0.50,'EdgeColor','none','HandleVisibility','off');
    fill(ax5L,[f;flipud(f)],[m_cn+se_cn;flipud(m_cn-se_cn)],c_case_fill,'FaceAlpha',0.50,'EdgeColor','none','HandleVisibility','off');
    fill(ax5L,[f;flipud(f)],[m_cc+se_cc;flipud(m_cc-se_cc)],c_case_fill,'FaceAlpha',0.50,'EdgeColor','none','HandleVisibility','off');

    plot(ax5L,f,m_tn,'-', 'Color',c_ctrl,'LineWidth',1.4,'DisplayName','Controls – No Chew');
    plot(ax5L,f,m_tc,'--','Color',c_ctrl,'LineWidth',1.4,'DisplayName','Controls – Chew');
    plot(ax5L,f,m_cn,'-', 'Color',c_case,'LineWidth',2.4,'DisplayName','Cases – No Chew');
    plot(ax5L,f,m_cc,'--','Color',c_case,'LineWidth',2.4,'DisplayName','Cases – Chew');

    if isfield(GR.Cases,'AP_Nc')
        m_ap_cn = nanmean(real(GR.Cases.AP_Nc),2);
        m_ap_cc = nanmean(real(GR.Cases.AP_Ch),2);
        m_ap_tn = nanmean(real(GR.Controls.AP_Nc),2);
        m_ap_tc = nanmean(real(GR.Controls.AP_Ch),2);
        plot(ax5L,f,m_ap_tn,':','Color',[c_ctrl 0.50],'LineWidth',1.0,'HandleVisibility','off');
        plot(ax5L,f,m_ap_tc,':','Color',[c_ctrl 0.50],'LineWidth',1.0,'HandleVisibility','off');
        plot(ax5L,f,m_ap_cn,':','Color',[c_case 0.50],'LineWidth',1.4,'HandleVisibility','off');
        plot(ax5L,f,m_ap_cc,':','Color',[c_case 0.50],'LineWidth',1.4,'HandleVisibility','off');
        text(ax5L,0.99,0.03,'Dotted: aperiodic fit','Units','normalized', ...
             'HorizontalAlignment','right','FontSize',8,'Color',[.5 .5 .5]);
    end

    set(ax5L,'XScale','log','XLim',[3 35],'YLim',ylims_raw, ...
        'XTick',[4 8 13 20 30],'XTickLabel',{'4','8','13','20','30'}, ...
        'TickDir','out','Box','off','XGrid','on','YGrid','on','GridAlpha',0.12,'Layer','top');
    xlabel(ax5L,'Frequency (Hz)','FontWeight','bold');
    ylabel(ax5L,'Power (dB)','FontWeight','bold');
    title(ax5L,'Broadband PSD (pre-FOOOF)','FontSize',12);
    legend(ax5L,'show','Location','southwest','Box','off','FontSize',8);

    % ── Panel derecho: Exponent boxplot ───────────────────────────────────
    ax5R = subplot(1,2,2); hold(ax5R,'on');

    draw_box(ax5R, exp_cas_nc, xNc, bw, clrNc,  0.55);
    draw_box(ax5R, exp_cas_ch, xCh, bw, clrCas, 0.55);

    rng(42);
    for i = 1:nCas
        line(ax5R,[xNc xCh],[exp_cas_nc(i) exp_cas_ch(i)], ...
             'Color',[.55 .55 .55 .25],'LineWidth',0.7);
    end
    scatter(ax5R, xNc+randn(nCas,1)*jit, exp_cas_nc, 30, clrNc,  'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
    scatter(ax5R, xCh+randn(nCas,1)*jit, exp_cas_ch, 30, clrCas, 'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
    

    line(ax5R,[xNc xNc xCh xCh],[yb-0.02 yb yb yb-0.02],'Color','k','LineWidth',1.2);
    text(ax5R,mean([xNc xCh]),yb+0.025,p2s(p_cas), ...
         'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

    set(ax5R,'XTick',[xNc xCh],'XTickLabel',{'No Chew','Chew'}, ...
             'XLim',[0.5 2.5],'YGrid','on','GridAlpha',0.20, ...
             'GridColor',[.6 .6 .6],'Box','off','TickDir','out','FontSize',11);
    ylabel(ax5R,'Aperiodic Exponent','FontWeight','bold');

    exportgraphics(fig5, fullfile(dir_out,'Fig2A_FOOOF_Panel.png'), 'Resolution', DPI);
    close(fig5);
    fprintf('    → Fig2A_FOOOF_Panel.png guardada.\n');
end

fprintf('✓ Script completado. Outputs en:\n  %s\n', dir_out);

%% ── FUNCIONES AUXILIARES ────────────────────────────────────────────────────
function draw_box(ax, data, xc, bw, clr, alpha)
    % Boxplot estándar: IQR box + mediana + bigotes con capuchón
    % Misma estructura visual que Panel A (Baseline Comparison)
    q       = quantile(data, [0.25 0.50 0.75]);
    iqr_val = q(3) - q(1);
    wlo     = max(data(data >= q(1) - 1.5*iqr_val));
    whi     = min(data(data <= q(3) + 1.5*iqr_val));
    hw      = bw / 2;
    cap     = hw * 0.55;   % ancho del capuchón del bigote

    % Caja IQR con relleno
    fill(ax, [xc-hw xc+hw xc+hw xc-hw xc-hw], ...
             [q(1)  q(1)  q(3)  q(3)  q(1)], ...
         clr, 'FaceAlpha', alpha, 'EdgeColor', clr*0.65, 'LineWidth', 1.4);

    % Línea de mediana
    line(ax, [xc-hw xc+hw], [q(2) q(2)], 'Color', 'w', 'LineWidth', 2.2);

    % Bigote inferior con capuchón
    line(ax, [xc xc],         [wlo q(1)],       'Color', clr*0.65, 'LineWidth', 1.2);
    line(ax, [xc-cap xc+cap], [wlo wlo],         'Color', clr*0.65, 'LineWidth', 1.2);

    % Bigote superior con capuchón
    line(ax, [xc xc],         [q(3) whi],        'Color', clr*0.65, 'LineWidth', 1.2);
    line(ax, [xc-cap xc+cap], [whi whi],         'Color', clr*0.65, 'LineWidth', 1.2);
end

function s = p2s(p)
    if     p < 0.001; s = '***';
    elseif p < 0.01;  s = '**';
    elseif p < 0.05;  s = '*';
    else;             s = 'n.s.';
    end
end