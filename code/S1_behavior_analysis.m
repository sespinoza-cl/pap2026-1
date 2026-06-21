%% S1_behavior_analysis.m
% Análisis conductual completo desde data_beh_tb.mat.
% Equivalente MATLAB del script Python S1_behavior_figs.py.
%
% Outputs:
%   figures/ — S1a..S1h (boxplots, deltas, descomposición IES)
%   stats/   — S1_behavior_results.mat  (deltas por sujeto para correlaciones)
%
% Uso posterior: correlacionar ies_delta / rt_delta con PAC, TF, chewing freq.

clear; clc;

ROOT    = 'C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1';
BEH_RAW = fullfile(ROOT, 'Analysis_V1_Final', 'data', 'computed', 'data_beh_tb.mat');
DIR_FIG = fullfile(ROOT, 'Analysis_V1_Final', 'outputs', 'figures');
DIR_STA = fullfile(ROOT, 'Analysis_V1_Final', 'outputs', 'stats');
if ~exist(DIR_FIG,'dir'), mkdir(DIR_FIG); end
if ~exist(DIR_STA,'dir'), mkdir(DIR_STA); end

%% ── 1. CARGAR DATOS ──────────────────────────────────────────────────────
tb = load(BEH_RAW);
d  = tb.tb_data;

rt_cas_nc  = d.casos.nochew.med(:);
rt_cas_ch  = d.casos.chew.med(:);
rt_ctr_nc  = d.controles.nochew.med(:);
rt_ctr_ch  = d.controles.chew.med(:);

ies_cas_nc = d.casos.nochew.ies(:);
ies_cas_ch = d.casos.chew.ies(:);
ies_ctr_nc = d.controles.nochew.ies(:);
ies_ctr_ch = d.controles.chew.ies(:);

acc_cas_nc = d.casos.nochew.acc(:);
acc_cas_ch = d.casos.chew.acc(:);
acc_ctr_nc = d.controles.nochew.acc(:);
acc_ctr_ch = d.controles.chew.acc(:);

rt_mean_cas_nc = d.casos.nochew.mean(:);
rt_mean_cas_ch = d.casos.chew.mean(:);
rt_mean_ctr_nc = d.controles.nochew.mean(:);
rt_mean_ctr_ch = d.controles.chew.mean(:);

n_cas = numel(rt_cas_nc);
n_ctr = numel(rt_ctr_nc);
fprintf('N = %d Cases  |  %d Controls\n', n_cas, n_ctr);

%% ── 2. DELTAS POR SUJETO ─────────────────────────────────────────────────
ies_delta_cas = ies_cas_ch - ies_cas_nc;
ies_delta_ctr = ies_ctr_ch - ies_ctr_nc;
rt_delta_cas  = rt_cas_ch  - rt_cas_nc;
rt_delta_ctr  = rt_ctr_ch  - rt_ctr_nc;
acc_delta_cas = acc_cas_ch - acc_cas_nc;
acc_delta_ctr = acc_ctr_ch - acc_ctr_nc;

ies_delta_pct_cas = ies_delta_cas ./ ies_cas_nc * 100;
ies_delta_pct_ctr = ies_delta_ctr ./ ies_ctr_nc * 100;
rt_delta_pct_cas  = rt_delta_cas  ./ rt_cas_nc  * 100;
rt_delta_pct_ctr  = rt_delta_ctr  ./ rt_ctr_nc  * 100;

%% ── 3. DESCOMPOSICIÓN IES ────────────────────────────────────────────────
% IES ≈ mean_RT / acc  →  ΔIES = RT_ch/acc_ch − RT_nc/acc_nc
% RT  contrib_i = (RT_ch − RT_nc) / acc_nc
% ACC contrib_i = RT_ch * (1/acc_ch − 1/acc_nc)
rt_contrib_cas  = (rt_mean_cas_ch - rt_mean_cas_nc) ./ acc_cas_nc;
acc_contrib_cas = rt_mean_cas_ch .* (1./acc_cas_ch - 1./acc_cas_nc);
rt_pct_cas      = rt_contrib_cas  ./ ies_cas_nc * 100;
acc_pct_cas     = acc_contrib_cas ./ ies_cas_nc * 100;

rt_contrib_ctr  = (rt_mean_ctr_ch - rt_mean_ctr_nc) ./ acc_ctr_nc;
acc_contrib_ctr = rt_mean_ctr_ch .* (1./acc_ctr_ch - 1./acc_ctr_nc);
rt_pct_ctr      = rt_contrib_ctr  ./ ies_ctr_nc * 100;
acc_pct_ctr     = acc_contrib_ctr ./ ies_ctr_nc * 100;

%% ── 4. ESTADÍSTICAS ──────────────────────────────────────────────────────
fprintf('\n=== ESTADÍSTICAS ===\n');

% Equivalencia baseline
[p_ies_b1] = ranksum(ies_cas_nc, ies_ctr_nc);
[p_rt_b1]  = ranksum(rt_cas_nc,  rt_ctr_nc);
fprintf('Baseline IES Cases vs Controls : p = %.3f\n', p_ies_b1);
fprintf('Baseline RT  Cases vs Controls : p = %.3f\n', p_rt_b1);

% Efecto masticatorio within (Wilcoxon signed-rank vs 0)
[p_ies_cas_w] = signrank(ies_delta_cas);
[p_rt_cas_w]  = signrank(rt_delta_cas);
[p_ies_ctr_w] = signrank(ies_delta_ctr);
[p_rt_ctr_w]  = signrank(rt_delta_ctr);
fprintf('Wilcoxon IES Cases vs 0        : p = %.4f\n', p_ies_cas_w);
fprintf('Wilcoxon RT  Cases vs 0        : p = %.4f\n', p_rt_cas_w);
fprintf('Wilcoxon IES Controls vs 0     : p = %.4f\n', p_ies_ctr_w);
fprintf('Wilcoxon RT  Controls vs 0     : p = %.4f\n', p_rt_ctr_w);

% Between-group delta
[p_ies_delta_mw] = ranksum(ies_delta_cas, ies_delta_ctr);
[p_rt_delta_mw]  = ranksum(rt_delta_cas,  rt_delta_ctr);
fprintf('Mann-Whitney ΔIES Cases vs Ctr : p = %.3f\n', p_ies_delta_mw);
fprintf('Mann-Whitney ΔRT  Cases vs Ctr : p = %.3f\n', p_rt_delta_mw);

% LME Group × Block
subj_id = [(1:n_cas)'; (1:n_cas)'; (n_cas+1:n_cas+n_ctr)'; (n_cas+1:n_cas+n_ctr)'];
grp_lbl = [repmat({'Cases'},    n_cas,1); repmat({'Cases'},    n_cas,1); ...
           repmat({'Controls'}, n_ctr,1); repmat({'Controls'}, n_ctr,1)];
blk_lbl = [ones(n_cas,1); 2*ones(n_cas,1); ones(n_ctr,1); 2*ones(n_ctr,1)];

T = table([rt_cas_nc; rt_cas_ch; rt_ctr_nc; rt_ctr_ch], ...
          [ies_cas_nc; ies_cas_ch; ies_ctr_nc; ies_ctr_ch], ...
          categorical(grp_lbl), categorical(blk_lbl), subj_id, ...
          'VariableNames', {'RT','IES','Group','Block','Subject'});

lme_rt  = fitlme(T, 'RT  ~ Group * Block + (1|Subject)');
lme_ies = fitlme(T, 'IES ~ Group * Block + (1|Subject)');

an_rt  = anova(lme_rt);
an_ies = anova(lme_ies);

idx_rt  = find(contains(an_rt.Term,  'Group') & contains(an_rt.Term,  'Block'));
idx_ies = find(contains(an_ies.Term, 'Group') & contains(an_ies.Term, 'Block'));
p_lme_rt  = an_rt.pValue(idx_rt);
p_lme_ies = an_ies.pValue(idx_ies);

fprintf('LME Group×Block RT  : F(%d,%d)=%.2f  p = %.4f\n', ...
    an_rt.DF1(idx_rt),  an_rt.DF2(idx_rt),  an_rt.FStat(idx_rt),  p_lme_rt);
fprintf('LME Group×Block IES : F(%d,%d)=%.2f  p = %.4f\n', ...
    an_ies.DF1(idx_ies), an_ies.DF2(idx_ies), an_ies.FStat(idx_ies), p_lme_ies);

% Descomposición IES
[p_tot_cas] = signrank(ies_delta_pct_cas);
[p_tot_ctr] = signrank(ies_delta_pct_ctr);
fprintf('\n=== IES DECOMPOSITION ===\n');
fprintf('Cases     RT=%.1f%%  ACC=%.1f%%  Total=%.1f%%  p=%.4f\n', ...
    mean(rt_pct_cas),  mean(acc_pct_cas),  mean(ies_delta_pct_cas), p_tot_cas);
fprintf('Controls  RT=%.1f%%  ACC=%.1f%%  Total=%.1f%%  p=%.4f\n', ...
    mean(rt_pct_ctr),  mean(acc_pct_ctr),  mean(ies_delta_pct_ctr), p_tot_ctr);

%% ── 5. GUARDAR RESULTADOS POR SUJETO ─────────────────────────────────────
beh_results.cas.ies_nc         = ies_cas_nc;
beh_results.cas.ies_ch         = ies_cas_ch;
beh_results.cas.ies_delta      = ies_delta_cas;
beh_results.cas.ies_delta_pct  = ies_delta_pct_cas;
beh_results.cas.rt_nc          = rt_cas_nc;
beh_results.cas.rt_ch          = rt_cas_ch;
beh_results.cas.rt_delta       = rt_delta_cas;
beh_results.cas.rt_delta_pct   = rt_delta_pct_cas;
beh_results.cas.acc_nc         = acc_cas_nc;
beh_results.cas.acc_ch         = acc_cas_ch;
beh_results.cas.acc_delta      = acc_delta_cas;

beh_results.ctr.ies_nc         = ies_ctr_nc;
beh_results.ctr.ies_ch         = ies_ctr_ch;
beh_results.ctr.ies_delta      = ies_delta_ctr;
beh_results.ctr.ies_delta_pct  = ies_delta_pct_ctr;
beh_results.ctr.rt_nc          = rt_ctr_nc;
beh_results.ctr.rt_ch          = rt_ctr_ch;
beh_results.ctr.rt_delta       = rt_delta_ctr;
beh_results.ctr.rt_delta_pct   = rt_delta_pct_ctr;
beh_results.ctr.acc_nc         = acc_ctr_nc;
beh_results.ctr.acc_ch         = acc_ctr_ch;
beh_results.ctr.acc_delta      = acc_delta_ctr;

beh_results.stats.p_ies_b1       = p_ies_b1;
beh_results.stats.p_rt_b1        = p_rt_b1;
beh_results.stats.p_ies_cas_w    = p_ies_cas_w;
beh_results.stats.p_rt_cas_w     = p_rt_cas_w;
beh_results.stats.p_ies_ctr_w    = p_ies_ctr_w;
beh_results.stats.p_rt_ctr_w     = p_rt_ctr_w;
beh_results.stats.p_lme_rt       = p_lme_rt;
beh_results.stats.p_lme_ies      = p_lme_ies;
beh_results.stats.p_ies_delta_mw = p_ies_delta_mw;
beh_results.stats.p_rt_delta_mw  = p_rt_delta_mw;

save(fullfile(DIR_STA, 'S1_behavior_results.mat'), 'beh_results');
fprintf('\nResultados guardados: S1_behavior_results.mat\n');

%% ── 6. FIGURAS ───────────────────────────────────────────────────────────
COL_CAS = [0.000, 0.620, 0.451];   % Okabe-Ito bluish-green  (Cases)
COL_CTR = [0.000, 0.447, 0.698];   % Okabe-Ito blue          (Controls)
COL_NC  = [0.533, 0.533, 0.533];   % gris neutro              (NoChew)

% ── FIG A: RT Baseline (B1) Cases vs Controls ─────────────────────────────
fh = figure('Units','inches','Position',[1 1 4 5],'Color','w');
ax = axes(fh); hold(ax,'on');
for g = 1:2
    if g==1, dat=rt_cas_nc; xc=1; col=COL_CAS; else, dat=rt_ctr_nc; xc=2; col=COL_CTR; end
    rng(g); jit = (rand(numel(dat),1)-0.5)*0.18;
    scatter(ax, xc+jit, dat, 24, 'filled', 'MarkerFaceColor',col, ...
        'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
    draw_box(ax, xc, dat, col, 0.22);
end
ym = max([rt_cas_nc; rt_ctr_nc])*1.25;
plot(ax,[1 2],[ym*0.90 ym*0.90],'-k','LineWidth',0.9);
text(ax,1.5, ym*0.95, sprintf('p = %.3f', p_rt_b1), ...
    'HorizontalAlignment','center','FontSize',9);
ax.XTick=[1 2]; ax.XTickLabel={'Cases','Controls'};
ax.XLim=[0.5 2.5]; ax.YLim=[0 ym];
ylabel(ax,'Median RT — NoChew block (ms)','FontSize',10);
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh, fullfile(DIR_FIG,'S1a_RT_baseline'), '-dpng','-r300'); close(fh);
fprintf('Guardada: S1a_RT_baseline.png  |  p=%.3f\n', p_rt_b1);

% ── FIG B: IES Cases (NoChew vs Chew) ────────────────────────────────────
fh = figure('Units','inches','Position',[1 1 4 5],'Color','w');
ax = axes(fh); hold(ax,'on');
draw_paired(ax, 1, 2, ies_cas_nc, ies_cas_ch, COL_NC, COL_CAS, 10);
ym = max([ies_cas_nc; ies_cas_ch])*1.22;
add_sig_bracket(ax, 1, 2, ym*0.90, ym*0.95, pstar(p_ies_cas_w));
ax.XTick=[1 2]; ax.XTickLabel={'NoChew','Chew'};
ax.XLim=[0.5 2.5]; ax.YLim=[0 ym];
ylabel(ax,'IES (ms)','FontSize',10);
title(ax,'Cases','FontSize',11,'FontWeight','bold');
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh,fullfile(DIR_FIG,'S1b_IES_Cases'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1b_IES_Cases.png  |  %s (p=%.4f)\n', pstar(p_ies_cas_w), p_ies_cas_w);

% ── FIG C: IES Controls ───────────────────────────────────────────────────
fh = figure('Units','inches','Position',[1 1 4 5],'Color','w');
ax = axes(fh); hold(ax,'on');
draw_paired(ax, 1, 2, ies_ctr_nc, ies_ctr_ch, COL_NC, COL_CTR, 20);
ym = max([ies_ctr_nc; ies_ctr_ch])*1.22;
add_sig_bracket(ax, 1, 2, ym*0.90, ym*0.95, pstar(p_ies_ctr_w));
ax.XTick=[1 2]; ax.XTickLabel={'NoChew','Chew'};
ax.XLim=[0.5 2.5]; ax.YLim=[0 ym];
ylabel(ax,'IES (ms)','FontSize',10);
title(ax,'Controls','FontSize',11,'FontWeight','bold');
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh,fullfile(DIR_FIG,'S1c_IES_Controls'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1c_IES_Controls.png  |  %s (p=%.4f)\n', pstar(p_ies_ctr_w), p_ies_ctr_w);

% ── FIG D: RT Cases ───────────────────────────────────────────────────────
fh = figure('Units','inches','Position',[1 1 4 5],'Color','w');
ax = axes(fh); hold(ax,'on');
draw_paired(ax, 1, 2, rt_cas_nc, rt_cas_ch, COL_NC, COL_CAS, 30);
ym = max([rt_cas_nc; rt_cas_ch])*1.22;
add_sig_bracket(ax, 1, 2, ym*0.90, ym*0.95, pstar(p_rt_cas_w));
ax.XTick=[1 2]; ax.XTickLabel={'NoChew','Chew'};
ax.XLim=[0.5 2.5]; ax.YLim=[0 ym];
ylabel(ax,'Median RT (ms)','FontSize',10);
title(ax,'Cases','FontSize',11,'FontWeight','bold');
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh,fullfile(DIR_FIG,'S1d_RT_Cases'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1d_RT_Cases.png  |  %s (p=%.4f)\n', pstar(p_rt_cas_w), p_rt_cas_w);

% ── FIG E: RT Controls ────────────────────────────────────────────────────
fh = figure('Units','inches','Position',[1 1 4 5],'Color','w');
ax = axes(fh); hold(ax,'on');
draw_paired(ax, 1, 2, rt_ctr_nc, rt_ctr_ch, COL_NC, COL_CTR, 40);
ym = max([rt_ctr_nc; rt_ctr_ch])*1.22;
add_sig_bracket(ax, 1, 2, ym*0.90, ym*0.95, pstar(p_rt_ctr_w));
ax.XTick=[1 2]; ax.XTickLabel={'NoChew','Chew'};
ax.XLim=[0.5 2.5]; ax.YLim=[0 ym];
ylabel(ax,'Median RT (ms)','FontSize',10);
title(ax,'Controls','FontSize',11,'FontWeight','bold');
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh,fullfile(DIR_FIG,'S1e_RT_Controls'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1e_RT_Controls.png  |  %s (p=%.4f)\n', pstar(p_rt_ctr_w), p_rt_ctr_w);

% ── FIG F: Δ IES (Cases y Controls) ──────────────────────────────────────
fh = figure('Units','inches','Position',[1 1 5 5],'Color','w');
ax = axes(fh); hold(ax,'on');
all_d = [ies_delta_cas; ies_delta_ctr];
ylo   = min(all_d)*1.30;
yhi   = abs(ylo)*0.30;                    % espacio positivo para brackets
plot(ax, [0.5 2.5], [0 0], '--', 'Color', [0.55 0.55 0.55], 'LineWidth',0.8);
for g = 1:2
    if g==1, dat=ies_delta_cas; xc=1; col=COL_CAS; pw=p_ies_cas_w;
    else,    dat=ies_delta_ctr; xc=2; col=COL_CTR; pw=p_ies_ctr_w; end
    rng(g*50); jit=(rand(numel(dat),1)-0.5)*0.18;
    scatter(ax, xc+jit, dat, 24, 'filled', 'MarkerFaceColor',col, ...
        'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
    draw_box(ax, xc, dat, col, 0.22);
    % estrella vs 0 (encima del whisker superior)
    whi = max(dat(dat <= prctile(dat,75) + 1.5*(prctile(dat,75)-prctile(dat,25))));
    text(ax, xc, whi + abs(ylo)*0.05, pstar(pw), ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
end
% bracket LME
bk_y = yhi * 0.70;
plot(ax,[1 2],[bk_y bk_y],'-k','LineWidth',0.9);
plot(ax,[1 1],[bk_y*0.85 bk_y],'-k','LineWidth',0.9);
plot(ax,[2 2],[bk_y*0.85 bk_y],'-k','LineWidth',0.9);
text(ax,1.5, bk_y*1.08, sprintf('LME p = %.3f', p_lme_ies), ...
    'HorizontalAlignment','center','FontSize',9);
ax.XTick=[1 2]; ax.XTickLabel={'Cases','Controls'};
ax.XLim=[0.5 2.5]; ax.YLim=[ylo yhi];
ylabel(ax,'Δ IES  (Chew − NoChew, ms)','FontSize',10);
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh,fullfile(DIR_FIG,'S1f_IES_delta'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1f_IES_delta.png  |  Cases %s  Ctrl %s  LME p=%.4f\n', ...
    pstar(p_ies_cas_w), pstar(p_ies_ctr_w), p_lme_ies);

% ── FIG G: Δ RT ───────────────────────────────────────────────────────────
fh = figure('Units','inches','Position',[1 1 5 5],'Color','w');
ax = axes(fh); hold(ax,'on');
all_d2 = [rt_delta_cas; rt_delta_ctr];
ylo2   = min(all_d2)*1.30;
yhi2   = abs(ylo2)*0.30;
plot(ax, [0.5 2.5], [0 0], '--', 'Color', [0.55 0.55 0.55], 'LineWidth',0.8);
for g = 1:2
    if g==1, dat=rt_delta_cas; xc=1; col=COL_CAS; pw=p_rt_cas_w;
    else,    dat=rt_delta_ctr; xc=2; col=COL_CTR; pw=p_rt_ctr_w; end
    rng(g*60); jit=(rand(numel(dat),1)-0.5)*0.18;
    scatter(ax, xc+jit, dat, 24, 'filled', 'MarkerFaceColor',col, ...
        'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
    draw_box(ax, xc, dat, col, 0.22);
    whi = max(dat(dat <= prctile(dat,75) + 1.5*(prctile(dat,75)-prctile(dat,25))));
    text(ax, xc, whi + abs(ylo2)*0.05, pstar(pw), ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
end
bk_y2 = yhi2 * 0.70;
plot(ax,[1 2],[bk_y2 bk_y2],'-k','LineWidth',0.9);
plot(ax,[1 1],[bk_y2*0.85 bk_y2],'-k','LineWidth',0.9);
plot(ax,[2 2],[bk_y2*0.85 bk_y2],'-k','LineWidth',0.9);
text(ax,1.5, bk_y2*1.08, sprintf('LME p = %.3f', p_lme_rt), ...
    'HorizontalAlignment','center','FontSize',9);
ax.XTick=[1 2]; ax.XTickLabel={'Cases','Controls'};
ax.XLim=[0.5 2.5]; ax.YLim=[ylo2 yhi2];
ylabel(ax,'Δ RT  (Chew − NoChew, ms)','FontSize',10);
ax.Box='off'; ax.TickDir='out'; ax.FontSize=9;
print(fh,fullfile(DIR_FIG,'S1g_RT_delta'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1g_RT_delta.png  |  Cases %s  Ctrl %s  LME p=%.4f\n', ...
    pstar(p_rt_cas_w), pstar(p_rt_ctr_w), p_lme_rt);

% ── FIG H: IES Shift Decomposition ───────────────────────────────────────
% Para descomposición IES: RT=color saturado del grupo, ACC=versión clara
COLS_RT  = [COL_CTR; COL_CAS];              % [Controls; Cases]
COLS_ACC = [0.48, 0.73, 0.86; ...           % azul claro  (Controls)
            0.44, 0.80, 0.69];              % verde claro (Cases)
rt_m  = [mean(rt_pct_ctr),  mean(rt_pct_cas)];
acc_m = [mean(acc_pct_ctr), mean(acc_pct_cas)];
tot_m = [mean(ies_delta_pct_ctr), mean(ies_delta_pct_cas)];

fh = figure('Units','inches','Position',[1 1 5 5.5],'Color','w');
ax = axes(fh); hold(ax,'on'); ax.Color = 'w';
bw = 0.30;
for i = 1:2
    xc = i - 1;
    % Barra RT (de 0 a rt_m) — color del grupo
    patch(ax, xc + [-bw bw bw -bw -bw], [0 0 rt_m(i) rt_m(i) 0], ...
          COLS_RT(i,:), 'EdgeColor','none');
    % Barra ACC (de rt_m a rt_m+acc_m, apilada) — versión clara del grupo
    bot = rt_m(i);
    top = rt_m(i) + acc_m(i);
    patch(ax, xc + [-bw bw bw -bw -bw], [bot bot top top bot], ...
          COLS_ACC(i,:), 'EdgeColor','none');
    % Diamante = total ΔIES
    plot(ax, xc, tot_m(i), 'd', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', ...
        'MarkerSize', 9);
    % Significancia debajo del diamante  (i=1→Controls, i=2→Cases)
    if i == 1, p_tot_i = p_tot_ctr; else, p_tot_i = p_tot_cas; end
    offset = abs(min(tot_m))*0.10;
    text(ax, xc, tot_m(i) - offset, pstar(p_tot_i), ...
        'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
end
% Línea cero
plot(ax, [-0.5 1.5], [0 0], '-k', 'LineWidth',0.9);
% Grid Y
ax.YGrid = 'on'; ax.GridLineStyle = '--'; ax.GridAlpha = 0.25;
% Leyenda
patch(ax,[NaN NaN NaN NaN],[NaN NaN NaN NaN], COLS_RT(1,:),  'DisplayName','RT Contribution');
patch(ax,[NaN NaN NaN NaN],[NaN NaN NaN NaN], COLS_ACC(1,:), 'DisplayName','ACC Contribution');
plot(ax, NaN, NaN, 'd', 'MarkerFaceColor','k','MarkerEdgeColor','k', ...
    'DisplayName','Total Δ IES (%)');
legend(ax,'Location','southwest','Box','off','FontSize',9);

ymin_d = min([rt_m + acc_m, tot_m]) * 1.30;
ax.XTick = [0 1]; ax.XTickLabel = {'Controls','Cases'};
ax.XLim = [-0.5 1.5]; ax.YLim = [ymin_d, 1];
ylabel(ax,'Δ IES (%)  (Chew − NoChew)','FontSize',11,'FontWeight','bold');
title(ax,'IES Shift Decomposition','FontSize',13,'FontWeight','bold');
ax.Box='off'; ax.TickDir='out'; ax.FontSize=10;
set(ax,'FontWeight','bold');
print(fh,fullfile(DIR_FIG,'S1h_IES_shift_decomposition'),'-dpng','-r300'); close(fh);
fprintf('Guardada: S1h_IES_shift_decomposition.png\n');
fprintf('\n=== FIN ===\n');

%% ── FUNCIONES LOCALES ────────────────────────────────────────────────────
function draw_box(ax, xc, dat, col, bw)
    q    = prctile(dat, [25 50 75]);
    iqrv = q(3) - q(1);
    wlo  = min(dat(dat >= q(1) - 1.5*iqrv));
    whi  = max(dat(dat <= q(3) + 1.5*iqrv));
    patch(ax, xc+[-bw bw bw -bw -bw], [q(1) q(1) q(3) q(3) q(1)], ...
          col, 'FaceAlpha',0.35, 'EdgeColor',col, 'LineWidth',1.4);
    plot(ax, [xc-bw xc+bw], [q(2) q(2)], '-w', 'LineWidth',2.0);
    plot(ax, [xc xc], [wlo q(1)], '-', 'Color',col, 'LineWidth',0.9);
    plot(ax, [xc xc], [q(3) whi], '-', 'Color',col, 'LineWidth',0.9);
    plot(ax, xc, mean(dat), '^k', 'MarkerSize',5, 'MarkerFaceColor','k');
end

function draw_paired(ax, x1, x2, dat1, dat2, col1, col2, seed)
    rng(seed);
    n   = numel(dat1);
    jit = (rand(n,1) - 0.5) * 0.18;
    for i = 1:n
        plot(ax, [x1+jit(i) x2+jit(i)], [dat1(i) dat2(i)], ...
             '-', 'Color', [0.7 0.7 0.7], 'LineWidth',0.6);
    end
    scatter(ax, x1+jit, dat1, 22, 'filled', 'MarkerFaceColor',col1, ...
        'MarkerFaceAlpha',0.55, 'MarkerEdgeColor','none');
    scatter(ax, x2+jit, dat2, 22, 'filled', 'MarkerFaceColor',col2, ...
        'MarkerFaceAlpha',0.55, 'MarkerEdgeColor','none');
    draw_box(ax, x1, dat1, col1, 0.20);
    draw_box(ax, x2, dat2, col2, 0.20);
end

function add_sig_bracket(ax, x1, x2, y_bar, y_txt, lbl)
    plot(ax, [x1 x2], [y_bar y_bar], '-k', 'LineWidth',0.9);
    plot(ax, [x1 x1], [y_bar*0.96 y_bar], '-k', 'LineWidth',0.9);
    plot(ax, [x2 x2], [y_bar*0.96 y_bar], '-k', 'LineWidth',0.9);
    text(ax, (x1+x2)/2, y_txt, lbl, ...
        'HorizontalAlignment','center', 'FontSize',12, 'FontWeight','bold');
end

function s = pstar(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'n.s.';
    end
end
