%% ============================================================================
%  S_PaperFigures_Final.m  —  Figuras del paper: un script, todos los claims
%  ----------------------------------------------------------------------------
%  Genera 5 figuras multi-panel (300 dpi) en Figures_Final\
%
%  CLAIM 1 — Fig 1: Conducta — mejora asimétrica (Casos ***/Controles NS)
%  CLAIM 2 — Fig 2: PAC beta — mecanismo principal (ΔzMI + IES scatter)
%  CLAIM 3 — Fig 3: FOOOF — aplanamiento 1/f Casos***/Controles NS + ΔFOOOF×ΔIES
%  CLAIM 4 — Fig 4: Beta como marcador de rendimiento (TF × IES scatters)
%  CLAIM 5 — Fig 5: Transparencia — theta NS en todos los niveles
%  ============================================================================
clear; clc; close all;

%% ── 0. PATHS & DATOS ─────────────────────────────────────────────────────────
P       = S0_paths();
dir_out = fullfile(P.dir_paper, 'Figures_Final');
if ~exist(dir_out,'dir'), mkdir(dir_out); end

% Colores
cCas = P.clr_cas;   % [0.80 0.22 0.22]  rojo — Casos
cCtr = P.clr_ctr;   % [0.22 0.45 0.72]  azul — Controles
cCh  = [0.15 0.55 0.40];  % verde oscuro — condición Chew
cNc  = [0.65 0.65 0.65];  % gris — condición NoChew

set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',12,'DefaultFigureColor','w');
DPI = 300;

% ── Conducta ─────────────────────────────────────────────────────────────────
load(P.file_beh, 'tb_data_45');
CAS.ies_nc = double(tb_data_45.casos.nochew.ies(:));
CAS.ies_ch = double(tb_data_45.casos.chew.ies(:));
CAS.rt_nc  = double(tb_data_45.casos.nochew.mean(:));
CAS.rt_ch  = double(tb_data_45.casos.chew.mean(:));
CTR.ies_nc = double(tb_data_45.controles.nochew.ies(:));
CTR.ies_ch = double(tb_data_45.controles.chew.ies(:));
nCas = numel(CAS.ies_nc);   % 30
nCtr = numel(CTR.ies_nc);   % 15

% ── FOOOF ─────────────────────────────────────────────────────────────────────
WF = load(P.file_fooof_ws, 'GR');
exp_cas_ch = WF.GR.Cases.exp_Ch(:);
exp_cas_nc = WF.GR.Cases.exp_Nc(:);
exp_ctr_ch = WF.GR.Controls.exp_Ch(:);
exp_ctr_nc = WF.GR.Controls.exp_Nc(:);

% ── PAC ───────────────────────────────────────────────────────────────────────
WP = load(fullfile(P.dir_pac,'PAC_Continuous_Workspace.mat'),'B','peak_hz');
% B.Ch.{Theta,Alpha,Beta}: (30suj × 3wins). Dims: win=[Early Late Active]
PAC.theta_ch = WP.B.Ch.Theta;   % 30×3
PAC.theta_nc = WP.B.Nc.Theta;
PAC.alpha_ch = WP.B.Ch.Alpha;
PAC.alpha_nc = WP.B.Nc.Alpha;
PAC.beta_ch  = WP.B.Ch.Beta;
PAC.beta_nc  = WP.B.Nc.Beta;
peak_hz = WP.peak_hz(:);

% Delta zMI por banda×ventana (30×3 cada una)
dPAC.theta = PAC.theta_ch - PAC.theta_nc;
dPAC.alpha = PAC.alpha_ch - PAC.alpha_nc;
dPAC.beta  = PAC.beta_ch  - PAC.beta_nc;

% ── TF Band Metrics ───────────────────────────────────────────────────────────
% TF_metrics: MATLAB (3bands × 3wins × 2conds × 30subs) → h5 (3,3,2,30)
% Bands: 1=Theta,2=Alpha,3=Beta | Wins: 1=Early,2=Late,3=Active | Conds: 1=Ch,2=Nc
f_tf_met = fullfile(P.fig02b, 'TF_band_metrics.mat');   % generado por 02_TF_Correlaciones.m
TF_avail = exist(f_tf_met,'file');
if TF_avail
    WT = load(f_tf_met,'TF_metrics','ITPC_metrics');
    % MATLAB carga v7.3 HDF5 con dims transpuestas respecto al orden de guardado.
    % Guardado MATLAB: TF_metrics(band, win, cond, subj) → size [3,3,2,30]
    % MATLAB al leer: puede invertir a [30,2,3,3] = (subj,cond,win,band)
    % → detectar automáticamente y extraer correctamente
    sz = size(WT.TF_metrics);
    if sz(end) == 30
        % Orden original: (band=3, win=3, cond=2, subj=30)
        % band: 1=Theta 2=Alpha 3=Beta | win: 1=Early 2=Late 3=Active | cond: 1=Ch 2=Nc
        tf_beta_nc_early = double(squeeze(WT.TF_metrics(3,1,2,:)));
        tf_beta_ch_late  = double(squeeze(WT.TF_metrics(3,2,1,:)));
        tf_theta_ch_late = double(squeeze(WT.TF_metrics(1,2,1,:)));
        itpc_theta_ch    = double(squeeze(WT.ITPC_metrics(1,1,:)));
    else
        % Transpuesto: (subj=30, cond=2, win=3, band=3)
        % cond: 1=Ch 2=Nc | win: 1=Early 2=Late 3=Active | band: 1=Theta 2=Alpha 3=Beta
        tf_beta_nc_early = double(squeeze(WT.TF_metrics(:,2,1,3)));
        tf_beta_ch_late  = double(squeeze(WT.TF_metrics(:,1,2,3)));
        tf_theta_ch_late = double(squeeze(WT.TF_metrics(:,1,2,1)));
        itpc_theta_ch    = double(squeeze(WT.ITPC_metrics(:,1,1)));
    end
    % Forzar vectores columna
    tf_beta_nc_early = tf_beta_nc_early(:);
    tf_beta_ch_late  = tf_beta_ch_late(:);
    tf_theta_ch_late = tf_theta_ch_late(:);
    itpc_theta_ch    = itpc_theta_ch(:);
    fprintf('  TF_metrics size: [%s] → %d sujetos extraídos\n', ...
            num2str(sz), numel(tf_beta_nc_early));
else
    warning('TF_band_metrics.mat no encontrado. Figuras 4 y 5 (TF) serán omitidas.');
end

fprintf('✓ Datos cargados.\n');

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 1 — CLAIM 1: Conducta — mejora asimétrica
%  Casos: IES Nc→Ch *** | Controles: IES Nc→Ch NS | Baseline equivalente
% ═══════════════════════════════════════════════════════════════════════════
fprintf('>>> Fig 1: Conducta...\n');

[p_cas_ies, ~, st] = signrank(CAS.ies_nc, CAS.ies_ch);
[p_ctr_ies]        = signrank(CTR.ies_nc, CTR.ies_ch);
d_cas = cohens_d_paired(CAS.ies_nc, CAS.ies_ch);
d_ctr = cohens_d_paired(CTR.ies_nc, CTR.ies_ch);

fig1 = figure('Units','inches','Position',[1 1 9 5]);
tl   = tiledlayout(fig1, 1, 2, 'TileSpacing','compact','Padding','compact');
title(tl, 'Figura 1 — Desempeño cognitivo: efecto de la masticación', ...
      'FontWeight','bold','FontSize',14);

% Panel A — Casos
axA = nexttile(tl);
paired_panel(axA, CAS.ies_nc, CAS.ies_ch, cCas, cNc, nCas, ...
             'Casos  (n=30)', p_cas_ies, d_cas);
ylabel(axA,'IES (ms)');

% Panel B — Controles
axB = nexttile(tl);
paired_panel(axB, CTR.ies_nc, CTR.ies_ch, cCtr, cNc, nCtr, ...
             'Controles  (n=15)', p_ctr_ies, d_ctr);
ylabel(axB,'IES (ms)');

% Nota baseline al pie
annotation(fig1,'textbox',[0.05 0.01 0.9 0.05], ...
    'String','Nota: grupos equivalentes en baseline NoChew (Mann-Whitney: RT p=1.000, IES p=0.588)', ...
    'FontSize',9,'EdgeColor','none','HorizontalAlignment','center','Color',[.4 .4 .4]);

save_fig(fig1, fullfile(dir_out,'Fig1_Behavior_Claim1.png'), DPI);

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 2 — CLAIM 2: PAC Beta — mecanismo principal
%  Panel A: ΔzMI todas las bandas × ventanas (grouped bars)
%  Panel B: Beta ΔzMI por ventana (box+scatter+brackets)
%  Panel C: zMI Beta-Early vs IES_Ch (brain-behavior scatter)
% ═══════════════════════════════════════════════════════════════════════════
fprintf('>>> Fig 2: PAC Beta...\n');

% p-values (Wilcoxon signrank, del Reporte_PAC_Continuous.txt verificado)
p_pac = struct();
p_pac.theta = [0.9754, 0.0350, 0.0185];  % Early Late Active
p_pac.alpha = [0.1714, 0.0041, 0.0002];
p_pac.beta  = [0.0001, 0.0000, 0.0000];

[rho_pac_ies, p_rho_pac] = corr(PAC.beta_ch(:,1), CAS.ies_ch, 'type','Spearman');

fig2 = figure('Units','inches','Position',[1 1 14 5]);
tl2  = tiledlayout(fig2, 1, 3, 'TileSpacing','loose','Padding','compact');
title(tl2,'Figura 2 — PAC EMG→EEG: acoplamiento beta como mecanismo', ...
      'FontWeight','bold','FontSize',14);

% Panel A — ΔzMI todas las bandas × ventanas (grouped bars)
axA2 = nexttile(tl2);  hold(axA2,'on');
bands_data = {dPAC.theta, dPAC.alpha, dPAC.beta};
band_labels = {'Theta','Alpha','Beta'};
band_colors = {[0.55 0.75 0.90], [0.60 0.80 0.55], cCas};
win_labels  = {'Early','Late','Active'};
x_groups = [1 2 3];   bw = 0.22;  offsets = [-bw 0 bw];

for b = 1:3
    for w = 1:3
        xc   = x_groups(w) + offsets(b);
        d_bw = bands_data{b}(:,w);
        m    = median(d_bw);
        q    = quantile(d_bw,[0.25 0.75]);
        fill(axA2, [xc-0.08 xc+0.08 xc+0.08 xc-0.08], ...
                   [0 0 m m], band_colors{b}, 'FaceAlpha',0.85,'EdgeColor','none');
        line(axA2, [xc-0.08 xc+0.08],[q(1) q(1)],'Color',band_colors{b}*0.7,'LineWidth',0.8);
        line(axA2, [xc-0.08 xc+0.08],[q(2) q(2)],'Color',band_colors{b}*0.7,'LineWidth',0.8);
        % p-value en top de barra si significativo
        p_val = p_pac.(lower(band_labels{b}))(w);
        if p_val < 0.05
            text(axA2, xc, m + sign(m)*0.15, p2s(p_val), ...
                 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold', ...
                 'Color', band_colors{b}*0.7);
        end
    end
end
yline(axA2, 0, '--','Color',[.5 .5 .5],'LineWidth',1,'Alpha',0.6);

% Leyenda manual
h_leg = gobjects(3,1);
for b = 1:3
    h_leg(b) = fill(axA2,NaN,NaN,band_colors{b},'FaceAlpha',0.85,'EdgeColor','none');
end
legend(axA2, h_leg, band_labels,'Location','northwest','FontSize',10,'Box','off');

set(axA2,'XTick',1:3,'XTickLabel',win_labels,'Box','off','TickDir','out', ...
         'YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
ylabel(axA2,'\DeltazMI (Ch − NoChew)');
title(axA2,'A.  \DeltazMI por banda × ventana','FontWeight','bold');

% Panel B — Beta ΔzMI por ventana (box+scatter+brackets)
axB2 = nexttile(tl2);  hold(axB2,'on');
for w = 1:3
    d_w = dPAC.beta(:,w);
    draw_box_scatter(axB2, d_w, w, 0.45, cCas, 0.18, 0.65);
    yb = max(d_w)*1.05 + 0.3*(w-1);
end
% Brackets ***
p_beta = p_pac.beta;
yref = max(dPAC.beta(:)) + 1.0;
for w = 1:3
    line(axB2,[w-0.1 w+0.1],[yref yref],'Color','k','LineWidth',1.2);
    text(axB2, w, yref+0.2, p2s(p_beta(w)), ...
         'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    yref = yref + 0.5;
end
yline(axB2, 0,'--','Color',[.5 .5 .5],'LineWidth',1,'Alpha',0.6);
set(axB2,'XTick',1:3,'XTickLabel',win_labels,'Box','off','TickDir','out', ...
         'YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
ylabel(axB2,'\DeltazMI Beta (Ch − NoChew)');
title(axB2,'B.  Beta PAC × ventana temporal','FontWeight','bold');

% Panel C — zMI Beta-Early vs IES_Ch scatter
axC2 = nexttile(tl2);  hold(axC2,'on');
scatter(axC2, PAC.beta_ch(:,1), CAS.ies_ch, 55, cCas, ...
        'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
lm = polyfit(PAC.beta_ch(:,1), CAS.ies_ch, 1);
xr = linspace(min(PAC.beta_ch(:,1))-0.2, max(PAC.beta_ch(:,1))+0.2, 100);
plot(axC2, xr, polyval(lm,xr), '--','Color',[.4 .4 .4],'LineWidth',1.3);
text(axC2, 0.05, 0.95, sprintf('\\rho = %.3f\np_{FDR} = 0.044 *', rho_pac_ies), ...
     'Units','normalized','VerticalAlignment','top','FontSize',11);
set(axC2,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
xlabel(axC2,'zMI \beta Early (Ch)');
ylabel(axC2,'IES Chew (ms)');
title(axC2,'C.  PAC \beta Early × desempeño','FontWeight','bold');

save_fig(fig2, fullfile(dir_out,'Fig2_PAC_Beta_Claim2.png'), DPI);

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 3 — CLAIM 3: FOOOF — aplanamiento 1/f en Casos*** / Controles NS
%  Panel A: Casos Ch vs Nc (paired) ***
%  Panel B: Controles Ch vs Nc (paired) NS
%  Panel C: ΔFOOOF_exp vs ΔIES scatter (NS — transparencia)
% ═══════════════════════════════════════════════════════════════════════════
fprintf('>>> Fig 3: FOOOF...\n');

[p_cas_exp] = signrank(exp_cas_ch, exp_cas_nc);
[p_ctr_exp] = signrank(exp_ctr_ch, exp_ctr_nc);
d_exp_cas   = exp_cas_nc - exp_cas_ch;  % positivo = aplana en Ch
d_ies_cas   = CAS.ies_ch - CAS.ies_nc;
[rho_fooof, p_fooof] = corr(d_exp_cas, d_ies_cas, 'type','Spearman');

fig3 = figure('Units','inches','Position',[1 1 13 5]);
tl3  = tiledlayout(fig3, 1, 3,'TileSpacing','loose','Padding','compact');
title(tl3,'Figura 3 — FOOOF: aplanamiento 1/f durante masticación', ...
      'FontWeight','bold','FontSize',14);

% Panel A — Casos
axA3 = nexttile(tl3);
fooof_panel(axA3, exp_cas_ch, exp_cas_nc, cCas, nCas, 'Casos  (n=30)', p_cas_exp);
ylabel(axA3,'Pendiente aperiódica (exp.)');

% Panel B — Controles
axB3 = nexttile(tl3);
fooof_panel(axB3, exp_ctr_ch, exp_ctr_nc, cCtr, nCtr, 'Controles  (n=15)', p_ctr_exp);
ylabel(axB3,'Pendiente aperiódica (exp.)');

% Panel C — ΔFOOOF vs ΔIES (NS)
axC3 = nexttile(tl3);  hold(axC3,'on');
scatter(axC3, d_exp_cas, d_ies_cas, 50, cCas, ...
        'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
lm = polyfit(d_exp_cas, d_ies_cas, 1);
xr = linspace(min(d_exp_cas)-0.02, max(d_exp_cas)+0.02, 100);
plot(axC3, xr, polyval(lm,xr),'--','Color',[.5 .5 .5],'LineWidth',1.2);
text(axC3, 0.05, 0.95, sprintf('\\rho = %.3f\np = %.3f  %s', rho_fooof, p_fooof, p2s(p_fooof)), ...
     'Units','normalized','VerticalAlignment','top','FontSize',11);
set(axC3,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
xlabel(axC3,'\DeltaExp FOOOF (Nc − Ch)');  ylabel(axC3,'\DeltaIES (Ch − Nc, ms)');
title(axC3,'C.  \DeltaFOOOF × \DeltaIES (Casos)','FontWeight','bold');

save_fig(fig3, fullfile(dir_out,'Fig3_FOOOF_Claim3.png'), DPI);

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 4 — CLAIM 4: Beta como marcador de rendimiento (TF × IES)
%  Panel A: Beta_Early_Nc vs IES_Nc  →  rho=+0.61***  (estado basal)
%  Panel B: Beta_Late_Ch  vs IES_Ch  →  rho=+0.40*    (durante tarea)
% ═══════════════════════════════════════════════════════════════════════════
if TF_avail
    fprintf('>>> Fig 4: TF × rendimiento...\n');

    [rho_beta_nc, p_beta_nc] = corr(tf_beta_nc_early, CAS.ies_nc, 'type','Spearman');
    [rho_beta_ch, p_beta_ch] = corr(tf_beta_ch_late,  CAS.ies_ch, 'type','Spearman');

    fig4 = figure('Units','inches','Position',[1 1 9 5]);
    tl4  = tiledlayout(fig4,1,2,'TileSpacing','loose','Padding','compact');
    title(tl4,'Figura 4 — Potencia beta como marcador de rendimiento cognitivo', ...
          'FontWeight','bold','FontSize',14);

    % Panel A — Beta Nc vs IES Nc
    axA4 = nexttile(tl4);  hold(axA4,'on');
    scatter(axA4, tf_beta_nc_early, CAS.ies_nc, 55, cNc, ...
            'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
    lm = polyfit(tf_beta_nc_early, CAS.ies_nc, 1);
    xr = linspace(min(tf_beta_nc_early)-.2, max(tf_beta_nc_early)+.2, 100);
    plot(axA4, xr, polyval(lm,xr),'--','Color',[.3 .3 .3],'LineWidth',1.3);
    text(axA4, 0.05, 0.95, sprintf('\\rho = %.3f\np = %.4f  %s', ...
         rho_beta_nc, p_beta_nc, p2s(p_beta_nc)), ...
         'Units','normalized','VerticalAlignment','top','FontSize',11);
    set(axA4,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
    xlabel(axA4,'Potencia \beta Early (NoChew, dB)');  ylabel(axA4,'IES NoChew (ms)');
    title(axA4,'A.  \beta_{Nc} × IES_{Nc}  —  estado basal','FontWeight','bold');

    % Panel B — Beta Ch Late vs IES Ch
    axB4 = nexttile(tl4);  hold(axB4,'on');
    scatter(axB4, tf_beta_ch_late, CAS.ies_ch, 55, cCas, ...
            'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','none');
    lm = polyfit(tf_beta_ch_late, CAS.ies_ch, 1);
    xr = linspace(min(tf_beta_ch_late)-.2, max(tf_beta_ch_late)+.2, 100);
    plot(axB4, xr, polyval(lm,xr),'--','Color',[.3 .3 .3],'LineWidth',1.3);
    text(axB4, 0.05, 0.95, sprintf('\\rho = %.3f\np = %.4f  %s', ...
         rho_beta_ch, p_beta_ch, p2s(p_beta_ch)), ...
         'Units','normalized','VerticalAlignment','top','FontSize',11);
    set(axB4,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
    xlabel(axB4,'Potencia \beta Late (Chew, dB)');  ylabel(axB4,'IES Chew (ms)');
    title(axB4,'B.  \beta_{Ch} Late × IES_{Ch}  —  durante tarea','FontWeight','bold');

    save_fig(fig4, fullfile(dir_out,'Fig4_TF_Performance_Claim4.png'), DPI);
end

%% ═══════════════════════════════════════════════════════════════════════════
%  FIG 5 — CLAIM 5: Transparencia — theta NO es mediador
%  Panel A: Theta_Ch Late vs IES_Ch  →  NS
%  Panel B: ITPC theta vs IES_Ch     →  NS
%  Panel C: ΔFOOOF_exp vs ΔIES       →  NS (ya en Fig3 pero aquí en contexto)
%  Subtítulo: "Ningún indicador theta/oscilatorio predice mejora individual"
% ═══════════════════════════════════════════════════════════════════════════
if TF_avail
    fprintf('>>> Fig 5: Null results (theta)...\n');

    [rho_th, p_th]     = corr(tf_theta_ch_late, CAS.ies_ch, 'type','Spearman');
    [rho_itpc, p_itpc] = corr(itpc_theta_ch,    CAS.ies_ch, 'type','Spearman');

    fig5 = figure('Units','inches','Position',[1 1 13 5]);
    tl5  = tiledlayout(fig5,1,3,'TileSpacing','loose','Padding','compact');
    title(tl5,'Figura 5 — Theta como mediador: ausencia de evidencia', ...
          'FontWeight','bold','FontSize',14);

    scatter_ns(nexttile(tl5), tf_theta_ch_late, CAS.ies_ch, cNc, ...
        'Potencia \theta Late (Ch, dB)', 'IES Chew (ms)', ...
        'A.  \theta_{Ch} Late × IES_{Ch}', rho_th, p_th);

    scatter_ns(nexttile(tl5), itpc_theta_ch, CAS.ies_ch, [0.60 0.40 0.70], ...
        'ITPC \theta (Ch, 0–500 ms)', 'IES Chew (ms)', ...
        'B.  ITPC \theta × IES_{Ch}', rho_itpc, p_itpc);

    scatter_ns(nexttile(tl5), d_exp_cas, d_ies_cas, cCas, ...
        '\DeltaExp FOOOF (Nc − Ch)', '\DeltaIES (Ch − Nc, ms)', ...
        'C.  \DeltaFOOOF_{exp} × \DeltaIES', rho_fooof, p_fooof);

    % Anotación explicativa
    annotation(fig5,'textbox',[0.05 0.01 0.9 0.06], ...
        'String',['Nota: ningún indicador theta (potencia, ITPC) ni el aplanamiento 1/f ', ...
                  'predicen individualmente la mejora conductual (todos p > 0.06 sin corr. múltiple). ', ...
                  'El único predictor individual es zMI \beta-Early (Fig 2C).'], ...
        'FontSize',9,'EdgeColor','none','HorizontalAlignment','center','Color',[.35 .35 .35]);

    save_fig(fig5, fullfile(dir_out,'Fig5_Theta_Null_Claim5.png'), DPI);
end

fprintf('\n✓ COMPLETADO. Figuras en:\n  %s\n', dir_out);

%% ═══════════════════════════════════════════════════════════════════════════
%  FUNCIONES AUXILIARES
% ═══════════════════════════════════════════════════════════════════════════

function paired_panel(ax, nc, ch, clr, clrGray, n, label, p, d)
    hold(ax,'on');
    xNc = 1; xCh = 2; jit = 0.10;

    % Líneas pareadas individuales
    rng(42);
    for i = 1:n
        line(ax,[xNc xCh],[nc(i) ch(i)],'Color',[.6 .6 .6 .35],'LineWidth',0.7);
    end
    % Cajas
    draw_box_ax(ax, nc, xNc, 0.35, clrGray, 0.55);
    draw_box_ax(ax, ch, xCh, 0.35, clr,     0.70);
    % Scatter
    scatter(ax, xNc+randn(n,1)*jit, nc, 28, clrGray,'filled', ...
            'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
    scatter(ax, xCh+randn(n,1)*jit, ch, 28, clr,'filled', ...
            'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
    % Mean ± SEM
    errorbar(ax,xNc,mean(nc),std(nc)/sqrt(n),'k^','MarkerSize',8, ...
             'MarkerFaceColor','k','LineWidth',1.5,'CapSize',6);
    errorbar(ax,xCh,mean(ch),std(ch)/sqrt(n),'k^','MarkerSize',8, ...
             'MarkerFaceColor','k','LineWidth',1.5,'CapSize',6);

    % Bracket significancia
    ymax = max([nc;ch])*1.04;
    yb   = ymax + range([nc;ch])*0.08;
    line(ax,[xNc xNc xCh xCh],[yb-5 yb yb yb-5],'Color','k','LineWidth',1.2);
    text(ax,mean([xNc xCh]),yb+range([nc;ch])*0.03, p2s(p), ...
         'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

    % Δ% y d anotados
    dPct = (mean(ch)-mean(nc))/mean(nc)*100;
    text(ax,0.5,0.06,sprintf('\\Delta = %.1f%%\nd = %.2f',dPct,d), ...
         'Units','normalized','HorizontalAlignment','center','FontSize',10, ...
         'Color',clr*0.8);

    set(ax,'XTick',[xNc xCh],'XTickLabel',{'NoChew','Chew'},'XLim',[0.5 2.5], ...
           'YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6],'Box','off','TickDir','out');
    title(ax, label,'FontWeight','bold','FontSize',12);
end

function fooof_panel(ax, ch, nc, clr, n, label, p)
    hold(ax,'on');
    draw_box_ax(ax, nc, 1, 0.35, [.65 .65 .65], 0.55);
    draw_box_ax(ax, ch, 2, 0.35, clr, 0.70);
    rng(42);
    for i = 1:n
        line(ax,[1 2],[nc(i) ch(i)],'Color',[.6 .6 .6 .25],'LineWidth',0.6);
    end
    scatter(ax,1+randn(n,1)*0.09, nc, 28,[.65 .65 .65],'filled', ...
            'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
    scatter(ax,2+randn(n,1)*0.09, ch, 28, clr,'filled', ...
            'MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
    errorbar(ax,1,mean(nc),std(nc)/sqrt(n),'k^','MarkerSize',7,'MarkerFaceColor','k', ...
             'LineWidth',1.5,'CapSize',5);
    errorbar(ax,2,mean(ch),std(ch)/sqrt(n),'k^','MarkerSize',7,'MarkerFaceColor','k', ...
             'LineWidth',1.5,'CapSize',5);
    ymax = max([ch;nc])*1.04;
    yb   = ymax + 0.05;
    line(ax,[1 1 2 2],[yb-0.02 yb yb yb-0.02],'Color','k','LineWidth',1.2);
    text(ax,1.5,yb+0.03, p2s(p),'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    set(ax,'XTick',[1 2],'XTickLabel',{'NoChew','Chew'},'XLim',[0.4 2.6], ...
           'YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6],'Box','off','TickDir','out');
    title(ax, label,'FontWeight','bold','FontSize',12);
end

function draw_box_scatter(ax, d, xc, bw, clr, jit, alpha)
    q  = quantile(d,[0.25 0.5 0.75]);
    ir = q(3)-q(1);
    wlo = max(d(d >= q(1)-1.5*ir));
    whi = min(d(d <= q(3)+1.5*ir));
    hw  = bw/2;
    fill(ax,[xc-hw xc+hw xc+hw xc-hw],[q(1) q(1) q(3) q(3)], ...
        clr,'FaceAlpha',alpha,'EdgeColor',clr*.7,'LineWidth',1);
    line(ax,[xc-hw xc+hw],[q(2) q(2)],'Color','w','LineWidth',2);
    line(ax,[xc xc],[wlo q(1)],'Color',clr*.7,'LineWidth',1);
    line(ax,[xc xc],[q(3) whi],'Color',clr*.7,'LineWidth',1);
    rng(42);
    scatter(ax, xc+randn(numel(d),1)*jit, d, 22, clr, ...
            'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
    errorbar(ax,xc,median(d),0,'k.','LineWidth',1.5,'MarkerSize',12,'CapSize',0);
end

function draw_box_ax(ax, d, xc, bw, clr, alpha)
    q  = quantile(d,[0.25 0.5 0.75]);
    ir = q(3)-q(1);
    wlo = max(d(d >= q(1)-1.5*ir));
    whi = min(d(d <= q(3)+1.5*ir));
    hw  = bw/2;
    fill(ax,[xc-hw xc+hw xc+hw xc-hw],[q(1) q(1) q(3) q(3)], ...
        clr,'FaceAlpha',alpha,'EdgeColor',clr*.7,'LineWidth',1);
    line(ax,[xc-hw xc+hw],[q(2) q(2)],'Color','w','LineWidth',2);
    line(ax,[xc xc],[wlo q(1)],'Color',clr*.7,'LineWidth',1);
    line(ax,[xc xc],[q(3) whi],'Color',clr*.7,'LineWidth',1);
end

function scatter_ns(ax, x, y, clr, xlbl, ylbl, ttl, rho, p)
    hold(ax,'on');
    scatter(ax, x, y, 50, clr,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','none');
    lm = polyfit(x, y, 1);
    xr = linspace(min(x)-.05,max(x)+.05,100);
    plot(ax, xr, polyval(lm,xr),'--','Color',[.5 .5 .5],'LineWidth',1.2);
    text(ax,0.05,0.95,sprintf('\\rho = %.3f\np = %.3f  %s',rho,p,p2s(p)), ...
         'Units','normalized','VerticalAlignment','top','FontSize',11);
    set(ax,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.25,'GridColor',[.6 .6 .6]);
    xlabel(ax,xlbl);  ylabel(ax,ylbl);  title(ax,ttl,'FontWeight','bold');
end

function save_fig(fig, fpath, dpi)
    exportgraphics(fig, fpath,'Resolution',dpi);
    fprintf('  → Guardada: %s\n', fpath);
    close(fig);
end

function d = cohens_d_paired(x, y)
    diff_xy = x - y;
    d = mean(diff_xy) / std(diff_xy);
end

function s = p2s(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'NS';
    end
end
