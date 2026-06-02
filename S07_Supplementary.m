%% ============================================================================
%  S_Supplementary_Plots.m
%  Supplementary figures + Rayleigh (preferred phase)
%
%  Generates:
%    FigS_A_Alpha_IES.png     — α power vs. IES scatter + CI
%    FigS_A_BetaPAC_IES.png  — β-PAC vs. IES scatter + CI
%    FigS_B_Mediacion.png    — Mediation: path diagram + bootstrap dist.
%    FigS_C_ITPC.png         — ITPC α: raw + partial scatter
%    FigS_D_Heatmaps.png     — EEG-PAC: ρ(zMI,IES) and median zMI heatmaps
%    FigS_E_Rayleigh.png     — Rayleigh: preferred phase by band × window
%
%  Sebastian, 2026
%% ============================================================================
clear; clc; close all;
fprintf('S_Supplementary_Plots — %s\n\n', datestr(now));

%% ── PATHS (desde S0_paths) ───────────────────────────────────────────────────
P         = S0_paths();
f_beh     = P.file_beh;
f_pac_emg = fullfile(P.fig04,  'PAC_4Groups_Workspace.mat');   % generado por 06_PAC.m
f_pac_eeg = fullfile(P.fig04,  'PAC_EEG_Workspace.mat');      % generado por S06b_PAC_EEG.m (θ-fase × amp)
f_tf      = fullfile(P.fig02b, 'TF_band_metrics.mat');         % generado por 02_TF_Correlaciones.m
dir_out   = P.dir_supp;
if ~exist(dir_out,'dir'); mkdir(dir_out); end

DPI = 300;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',11,'DefaultFigureColor','w');

% Colour palette
c_cas  = [0.18 0.55 0.80];   % blue — cases
c_gray = [0.55 0.55 0.55];
c_red  = [0.82 0.22 0.22];
c_grn  = [0.20 0.65 0.40];

NBOOT = 2000;   % bootstrap iterations (shared)
rng(2026);      % reproducibility

%% ── LOAD DATA ────────────────────────────────────────────────────────────────
% ── Behavioural table ─────────────────────────────────────────────────────────
beh    = load(f_beh); fn = fieldnames(beh); D = beh.(fn{1});
BEH_Ch = D.casos.chew;
BEH_Nc = D.casos.nochew;

beh_ids = string(BEH_Ch.Participantes);
beh_ids = arrayfun(@(s) extractBefore(s + "_X","_X"), beh_ids);
need_e3 = ~startsWith(beh_ids,"E3");
beh_ids(need_e3) = "E3" + beh_ids(need_e3);

get_ies = @(ids) arrayfun(@(s) ...
    ternary_val(any(beh_ids==s), double(BEH_Ch.ies_m(find(beh_ids==s,1))), NaN), ids);
get_ies_nc = @(ids) arrayfun(@(s) ...
    ternary_val(any(beh_ids==s), double(BEH_Nc.ies_m(find(beh_ids==s,1))), NaN), ids);

% ── β-PAC workspace (D: drive — same file as Fig2 Panel C) ──────────────────
W_pac      = load(f_pac_emg, 'B', 'casos');
casos_pac  = string(W_pac.casos(:));
beta_early = W_pac.B.Casos.Ch.Beta(:, 1);   % β-Early [nPAC × 1]
IES_pac    = get_ies(casos_pac);             % IES aligned to PAC order
IES_nc_pac = get_ies_nc(casos_pac);
fprintf('  β-PAC: %d sujetos | IES matched: %d\n', numel(casos_pac), sum(~isnan(IES_pac)));

% ── TF workspace ─────────────────────────────────────────────────────────────
tf         = load(f_tf);
TFM        = tf.TF_metrics;         % [nTF 2 3 3]
casos_tf   = string(tf.casos_ids(:));
alpha_act_raw = squeeze(TFM(:,1,2,3));   % α-Mid [nTF × 1]
IES_tf     = get_ies(casos_tf);          % IES aligned to TF order
fprintf('  α-TF : %d sujetos | IES matched: %d\n', numel(casos_tf), sum(~isnan(IES_tf)));

% ── For individual-figure variables: use their own aligned IES ────────────────
% (IES_ch / alpha_act / beta_early used in scatter panels, mediation, etc.)
alpha_act = alpha_act_raw;
IES_ch    = IES_pac;    % default: PAC order (used for most figures)
IES_nc    = IES_nc_pac;
N         = numel(casos_pac);

% ── ITPC & EEG-PAC workspace ─────────────────────────────────────────────────
ITPCM      = tf.ITPC_metrics;          % [nTF 1 3]
alpha_itpc_raw = squeeze(ITPCM(:,1,2));

have_eeg = exist(f_pac_eeg,'file')==2;   % EEG-EEG PAC opcional (comodulograma θ-fase × amp)
if have_eeg
    E        = load(f_pac_eeg);
    zMI_eeg  = E.zMI_ch;          % [30 24 3]
    ies_eeg  = double(E.IES_Ch(:));
    pref_ch  = E.pref_ch;
    pref_nc  = E.pref_nc;
else
    fprintf(2, ['[SKIP] EEG-EEG PAC workspace no encontrado → se omiten FigS_D ' ...
        '(heatmaps) y FigS_E (Rayleigh theta-fase). Steiger/mediacion/ITPC corren igual.\n  %s\n'], f_pac_eeg);
    zMI_eeg = []; ies_eeg = []; pref_ch = []; pref_nc = [];
end

% ── Three-way intersection for Steiger (α, β-PAC, IES must be same subjects) ─
common_ids = intersect(casos_pac, casos_tf);
nC = numel(common_ids);
alpha_act_3 = nan(nC,1);
beta_early_3 = nan(nC,1);
IES_3        = nan(nC,1);
for j = 1:nC
    s = common_ids(j);
    i_tf  = find(casos_tf  == s, 1);
    i_pac = find(casos_pac == s, 1);
    if ~isempty(i_tf);  alpha_act_3(j)  = alpha_act_raw(i_tf);  end
    if ~isempty(i_pac); beta_early_3(j) = beta_early(i_pac);    end
    ix_b = find(beh_ids == s, 1);
    if ~isempty(ix_b);  IES_3(j) = double(BEH_Ch.ies_m(ix_b)); end
end
fprintf('  Steiger common N = %d (α∩β-PAC∩IES)\n\n', sum(~isnan(IES_3) & ~isnan(alpha_act_3) & ~isnan(beta_early_3)));

% alpha_itpc — aligned to TF subject order
alpha_itpc = arrayfun(@(s) ...
    ternary_val(any(casos_tf==s), alpha_itpc_raw(find(casos_tf==s,1)), NaN), casos_pac);

amp_freqs = 8:4:100;           % 24 frequencies
wins_lbl  = {'Early','Late','Mid'};
amp_bands = struct('name', {'LowAlpha','HighAlpha','Beta','LowGamma','HighGamma','BroadGamma'}, ...
                   'range',{[8 13],[14 20],[20 30],[30 50],[50 80],[30 100]});
n_ab      = numel(amp_bands);
BAND_CLRS = lines(n_ab);

%% ── PRE-COMPUTE CORRELATIONS & STEIGER ───────────────────────────────────────
% r12 & r13 computed on the SAME subjects (3-way intersection) with SAME IES
v3   = ~isnan(alpha_act_3) & ~isnan(beta_early_3) & ~isnan(IES_3);
N_st = sum(v3);

[r12, p12] = corr(alpha_act_3(v3),  IES_3(v3),         'Type','Spearman');
[r13, p13] = corr(beta_early_3(v3), IES_3(v3),         'Type','Spearman');
[r23, p23] = corr(alpha_act_3(v3),  beta_early_3(v3),  'Type','Spearman');
[Z_st, p_st] = steiger_test(r12, r13, r23, N_st);

fprintf('  r12 (α, IES)    = %+.3f   p = %.4f  %s\n', r12, p12, sig_sym(p12));
fprintf('  r13 (β-PAC,IES) = %+.3f   p = %.4f  %s\n', r13, p13, sig_sym(p13));
fprintf('  r23 (α, β-PAC)  = %+.3f   p = %.4f  %s\n', r23, p23, sig_sym(p23));
fprintf('  Steiger Z = %+.3f   p = %.4f  %s\n\n', Z_st, p_st, sig_sym(p_st));

%% ════════════════════════════════════════════════════════════════════════════
%  FIG S_A_1 — α power vs. IES  (with bootstrap CI)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('Generating FigS_A_Alpha_IES.png...\n');

% α-Active (TF order) × IES_tf — both aligned to TF subject list
v_tf  = ~isnan(alpha_act_raw) & ~isnan(IES_tf);
x1 = alpha_act_raw(v_tf);
y1 = IES_tf(v_tf);
n1 = numel(x1);

[rho_a1, p_a1] = corr(x1, y1, 'Type','Spearman');

% Regression line
Pf1  = polyfit(x1, y1, 1);
xl1  = linspace(min(x1), max(x1), 200);
yy1  = polyval(Pf1, xl1);

% Bootstrap CI: regression band 
rb1   = nan(NBOOT,1);
pb1   = nan(NBOOT, numel(xl1));
for k = 1:NBOOT
    idx_b     = randi(n1, n1, 1);
    rb1(k)    = corr(x1(idx_b), y1(idx_b), 'Type','Spearman');
    pb1(k,:)  = polyval(polyfit(x1(idx_b), y1(idx_b), 1), xl1);
end
ci1_lo = prctile(pb1, 2.5);
ci1_hi = prctile(pb1, 97.5);

fig1 = figure('Color','w','Units','inches','Position',[0 0 4 4]);
ax1  = axes('Parent',fig1); hold(ax1,'on');

fill(ax1, [xl1 fliplr(xl1)], [ci1_lo fliplr(ci1_hi)], c_cas, ...
    'FaceAlpha',0.15, 'EdgeColor','none');
plot(ax1, xl1, yy1, '-', 'Color',c_cas, 'LineWidth',1.8);
scatter(ax1, x1, y1, 50, c_cas, 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerEdgeColor','w', 'LineWidth',0.6);

% Texto modificado: solo rho y p-valor
text(ax1, 0.05, 0.95, ...
    sprintf('\\rho = %+.2f\np = %.3f%s', rho_a1, p_a1, sig_sym(p_a1)), ...
    'Units','norm', 'VerticalAlignment','top', 'FontSize',10, ...
    'FontWeight','bold', 'Color', c_cas*0.8);

xlabel(ax1, '\alpha_{Mid} power (dB)', 'FontSize',11);
ylabel(ax1, 'IES (ms)',                   'FontSize',11);
title(ax1,  '\alpha power vs. IES',       'FontWeight','bold','FontSize',12);
box(ax1,'off'); set(ax1,'TickDir','out','FontSize',10);

exportgraphics(fig1, fullfile(dir_out,'FigS_A_Alpha_IES.png'), 'Resolution',DPI);
close(fig1);

%% ════════════════════════════════════════════════════════════════════════════
%  FIG S_A_2 — β-PAC vs. IES  (with bootstrap CI)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('Generating FigS_A_BetaPAC_IES.png...\n');

% β-PAC × IES_pac — both aligned to PAC subject order
v_be = ~isnan(beta_early) & ~isnan(IES_pac);
x2 = beta_early(v_be);
y2 = IES_pac(v_be);
n2 = numel(x2);
[rho_b2, p_b2] = corr(x2, y2, 'Type','Spearman');

Pf2 = polyfit(x2, y2, 1);
xl2 = linspace(min(x2), max(x2), 200);
yy2 = polyval(Pf2, xl2);

rb2  = nan(NBOOT,1);
pb2  = nan(NBOOT, numel(xl2));
for k = 1:NBOOT
    idx_b    = randi(n2, n2, 1);
    rb2(k)   = corr(x2(idx_b), y2(idx_b), 'Type','Spearman');
    pb2(k,:) = polyval(polyfit(x2(idx_b), y2(idx_b), 1), xl2);
end
ci2_lo  = prctile(pb2, 2.5);
ci2_hi  = prctile(pb2, 97.5);
ci_rho2 = prctile(rb2, [2.5 97.5]);

fig2 = figure('Color','w','Units','inches','Position',[0 0 4 4]);
ax2  = axes('Parent',fig2); hold(ax2,'on');

fill(ax2, [xl2 fliplr(xl2)], [ci2_lo fliplr(ci2_hi)], c_gray, ...
    'FaceAlpha',0.15, 'EdgeColor','none');
plot(ax2, xl2, yy2, '-', 'Color',c_gray, 'LineWidth',1.8);
scatter(ax2, x2, y2, 50, c_gray, 'filled', 'MarkerFaceAlpha',0.75, ...
    'MarkerEdgeColor','w', 'LineWidth',0.6);

text(ax2, 0.05, 0.95, ...
    sprintf('\\rho = %+.2f\n95%% CI [%+.2f, %+.2f]\np = %.3f%s\nN = %d', ...
    rho_b2, ci_rho2(1), ci_rho2(2), p_b2, sig_sym(p_b2), n2), ...
    'Units','norm', 'VerticalAlignment','top', 'FontSize',10, ...
    'FontWeight','bold', 'Color', [0.35 0.35 0.35]);

xlabel(ax2, '\beta-PAC_{Early} (zMI)', 'FontSize',11);
ylabel(ax2, 'IES (ms)',                'FontSize',11);
title(ax2,  '\beta-PAC vs. IES',       'FontWeight','bold','FontSize',12);
box(ax2,'off'); set(ax2,'TickDir','out','FontSize',10);

exportgraphics(fig2, fullfile(dir_out,'FigS_A_BetaPAC_IES.png'), 'Resolution',DPI);
close(fig2);

%% ════════════════════════════════════════════════════════════════════════════
%  FIG B — Mediation: bootstrap + path diagram
%% ════════════════════════════════════════════════════════════════════════════
fprintf('Generating FigS_B (Mediation bootstrap)...\n');

% Use 3-way aligned vectors (same subjects, same IES) for valid mediation
X = beta_early_3; M = alpha_act_3; Y = IES_3;
v = ~isnan(X) & ~isnan(M) & ~isnan(Y);
X = X(v); M = M(v); Y = Y(v); n = sum(v);

a   = ols_slope(X,M);
bv  = ols_slopes([M X], Y); b = bv(1); cp = bv(2);
c   = ols_slope(X,Y);
ab  = a*b;

ab_boot = nan(NBOOT,1);
for k = 1:NBOOT
    idx = randsample(n,n,true);
    bv_k      = ols_slopes([M(idx) X(idx)], Y(idx));
    ab_boot(k) = ols_slope(X(idx),M(idx)) * bv_k(1);
end
ci_lo = prctile(ab_boot,2.5);
ci_hi = prctile(ab_boot,97.5);
sig_boot = (ci_lo>0 && ci_hi>0) || (ci_lo<0 && ci_hi<0);

figB = figure('Units','inches','Position',[0 0 10 4.5]);
tl2  = tiledlayout(figB, 1, 2, 'TileSpacing','compact','Padding','compact');

% Left panel: bootstrap distribution of ab
ax3 = nexttile; hold(ax3,'on');
histogram(ax3, ab_boot, 60, 'FaceColor',c_cas,'EdgeColor','none', ...
    'FaceAlpha',0.7,'Normalization','probability');
xline(ax3, 0,    '--', 'Color','k',   'LineWidth',1.5);
xline(ax3, ab,   '-',  'Color',c_red, 'LineWidth',2.5, ...
    'Label',sprintf('ab=%.3f',ab), 'LabelHorizontalAlignment','left');
xline(ax3, ci_lo, ':', 'Color',c_red, 'LineWidth',1.8);
xline(ax3, ci_hi, ':', 'Color',c_red, 'LineWidth',1.8, ...
    'Label',sprintf('95%% CI [%.2f, %.2f]',ci_lo,ci_hi));
xlabel(ax3, 'Indirect effect ab (bootstrap)', 'FontSize',11);
ylabel(ax3, 'Probability',                    'FontSize',11);
title(ax3,  'Bootstrap ab (5000 iter)',        'FontWeight','bold');
text(ax3, 0.97, 0.95, ternary(sig_boot,'✓ CI excludes 0','✗ CI crosses 0'), ...
    'Units','norm','HorizontalAlignment','right','VerticalAlignment','top','FontSize',11);
box(ax3,'off'); set(ax3,'TickDir','out');

% Right panel: path diagram
ax4 = nexttile;
axis(ax4,'off'); hold(ax4,'on');
set(ax4,'XLim',[0 10],'YLim',[0 6]);

node_X   = [1.5  5.0  8.5];
node_Y   = [3.0  5.2  3.0];
node_lbl = {'\beta_{PAC-Early}', '\alpha_{Mid}', 'IES_{Ch}'};
node_sub = {'(X)', '(M)', '(Y)'};
node_c   = {c_gray, c_cas, c_red};
for ni = 1:3
    rectangle(ax4, 'Position',[node_X(ni)-1.1 node_Y(ni)-0.55 2.2 1.1], ...
        'Curvature',0.3, 'FaceColor',node_c{ni}, 'EdgeColor','none');
    text(ax4, node_X(ni), node_Y(ni)+0.15, node_lbl{ni}, ...
        'HorizontalAlignment','center','Color','w','FontWeight','bold','FontSize',10);
    text(ax4, node_X(ni), node_Y(ni)-0.20, node_sub{ni}, ...
        'HorizontalAlignment','center','Color','w','FontSize',9);
end
annotate_arrow(ax4, [node_X(1)+1.1 node_X(2)-1.1], [node_Y(1)+0.3 node_Y(2)-0.3], sprintf('a=%.4f',a),   'k');
annotate_arrow(ax4, [node_X(2)+1.1 node_X(3)-1.1], [node_Y(2)-0.3 node_Y(3)+0.3], sprintf('b=%.2f',b),   'k');
annotate_arrow(ax4, [node_X(1)+1.1 node_X(3)-1.1], [node_Y(1)-0.1 node_Y(3)-0.1], sprintf('c''=%.2f',cp), c_gray);
text(ax4, 5, 1.5, sprintf('ab = %.4f\n95%% CI [%.3f, %.3f]\n%s', ab, ci_lo, ci_hi, ...
    ternary(sig_boot,'✓ Significant','✗ n.s.')), ...
    'HorizontalAlignment','center','FontSize',10,'Color',c_red,'FontWeight','bold');
title(ax4, 'Path diagram (OLS)', 'FontWeight','bold');

exportgraphics(figB, fullfile(dir_out,'FigS_B_Mediacion.png'), 'Resolution',DPI);
close(figB);

%% ════════════════════════════════════════════════════════════════════════════
%  FIG C — ITPC α: correlations (raw and partial)
%% ════════════════════════════════════════════════════════════════════════════
fprintf('Generating FigS_C (ITPC alpha)...\n');

[r_i, p_i] = corr(alpha_itpc, IES_ch, 'Type','Spearman');
rp          = partial_spearman(alpha_itpc, IES_ch, alpha_act);
df_p        = N-3;
t_p         = rp * sqrt(df_p / max(1e-10, 1-rp^2));
p_p         = 2*(1-tcdf(abs(t_p), df_p));

figC = figure('Units','inches','Position',[0 0 10 4.5]);
tl3  = tiledlayout(figC, 1, 2, 'TileSpacing','compact','Padding','compact');

% Left: raw scatter
ax5 = nexttile; hold(ax5,'on');
scatter(ax5, alpha_itpc, IES_ch, 60, c_cas, 'filled', 'MarkerFaceAlpha',0.8);
lmi = polyfit(alpha_itpc, IES_ch, 1);
xli = linspace(min(alpha_itpc)-0.005, max(alpha_itpc)+0.005, 80);
plot(ax5, xli, polyval(lmi,xli), '--', 'Color',c_gray, 'LineWidth',1.5);
text(ax5, 0.05, 0.95, sprintf('\\rho = %+.3f\np = %.3f %s', r_i, p_i, sig_sym(p_i)), ...
    'Units','norm','VerticalAlignment','top','FontSize',11);
xlabel(ax5, '\alpha-ITPC (Chew)',  'FontSize',12);
ylabel(ax5, 'IES Chew (ms)',       'FontSize',12);
title(ax5,  '\alpha-ITPC vs. IES (raw)', 'FontWeight','bold');
box(ax5,'off'); set(ax5,'TickDir','out');

% Right: partial scatter (residuals controlling α power)
ax6 = nexttile; hold(ax6,'on');
res_itpc = alpha_itpc - polyval(polyfit(alpha_act,alpha_itpc,1), alpha_act);
res_ies  = IES_ch     - polyval(polyfit(alpha_act,IES_ch,1),    alpha_act);
scatter(ax6, res_itpc, res_ies, 60, c_cas, 'filled', 'MarkerFaceAlpha',0.8);
lmr = polyfit(res_itpc, res_ies, 1);
xlr = linspace(min(res_itpc)-0.002, max(res_itpc)+0.002, 80);
plot(ax6, xlr, polyval(lmr,xlr), '--', 'Color',c_gray, 'LineWidth',1.5);
xline(ax6, 0, '-', 'Color',[.8 .8 .8], 'LineWidth',0.8);
yline(ax6, 0, '-', 'Color',[.8 .8 .8], 'LineWidth',0.8);
text(ax6, 0.05, 0.95, sprintf('\\rho_{partial} = %+.3f\np = %.3f %s', rp, p_p, sig_sym(p_p)), ...
    'Units','norm','VerticalAlignment','top','FontSize',11);
xlabel(ax6, 'Residual \alpha-ITPC | \alpha-power', 'FontSize',12);
ylabel(ax6, 'Residual IES | \alpha-power',         'FontSize',12);
title(ax6,  '\alpha-ITPC vs. IES (partial)',        'FontWeight','bold');
box(ax6,'off'); set(ax6,'TickDir','out');

exportgraphics(figC, fullfile(dir_out,'FigS_C_ITPC.png'), 'Resolution',DPI);
close(figC);

%% ════════════════════════════════════════════════════════════════════════════
%  FIG D — EEG-PAC heatmaps (median zMI and ρ vs IES)   [requiere EEG-EEG PAC]
%% ════════════════════════════════════════════════════════════════════════════
if have_eeg
fprintf('Generating FigS_D (heatmaps EEG-PAC)...\n');

rho_mat = nan(n_ab,3); p_mat = nan(n_ab,3); med_mat = nan(n_ab,3);
for ab = 1:n_ab
    fi = amp_freqs >= amp_bands(ab).range(1) & amp_freqs <= amp_bands(ab).range(2);
    if ~any(fi); continue; end
    for w = 1:3
        zw = squeeze(mean(zMI_eeg(:,fi,w), 2, 'omitnan'));
        med_mat(ab,w) = median(zw,'omitnan');
        vv = ~isnan(zw) & ~isnan(ies_eeg);
        if sum(vv) >= 5
            [rho_mat(ab,w), p_mat(ab,w)] = corr(double(zw(vv)), ies_eeg(vv), 'Type','Spearman');
        end
    end
end
band_names = {amp_bands.name};

figD = figure('Units','inches','Position',[0 0 11 4.5]);
tl4  = tiledlayout(figD, 1, 2, 'TileSpacing','compact','Padding','compact');

% Heatmap: median zMI
ax7 = nexttile;
imagesc(ax7, med_mat);
clim_m = max(abs(med_mat(:))); clim(ax7, [-clim_m clim_m]);
colormap(ax7, div_cmap(256,[0.8 0.2 0.2],[0.2 0.5 0.8]));
cb7 = colorbar(ax7); cb7.Label.String = 'median zMI';
set(ax7,'YTick',1:n_ab,'YTickLabel',band_names, ...
    'XTick',1:3,'XTickLabel',wins_lbl,'TickDir','out','Box','on');
title(ax7, 'Median zMI_{Ch} (\theta-phase \times band)', 'FontWeight','bold');
for ab = 1:n_ab; for w = 1:3
    [pw,~] = signrank(double(squeeze(mean(zMI_eeg(:, ...
        amp_freqs>=amp_bands(ab).range(1)&amp_freqs<=amp_bands(ab).range(2), w), 2,'omitnan'))));
    if pw < 0.05
        text(ax7, w, ab, sig_sym(pw), 'HorizontalAlignment','center', ...
            'FontSize',14,'FontWeight','bold','Color','w');
    end
end; end

% Heatmap: ρ vs IES
ax8 = nexttile;
imagesc(ax8, rho_mat);
clim_r = max(abs(rho_mat(~isnan(rho_mat))));
if isempty(clim_r) || clim_r == 0; clim_r = 0.5; end
clim(ax8, [-clim_r clim_r]);
colormap(ax8, div_cmap(256,[0.8 0.2 0.2],[0.2 0.5 0.8]));
cb8 = colorbar(ax8); cb8.Label.String = '\rho (IES_{Ch})';
set(ax8,'YTick',1:n_ab,'YTickLabel',band_names, ...
    'XTick',1:3,'XTickLabel',wins_lbl,'TickDir','out','Box','on');
title(ax8, '\rho (zMI_{Ch}, IES_{Ch}) — Spearman', 'FontWeight','bold');
for ab = 1:n_ab; for w = 1:3
    if ~isnan(p_mat(ab,w)) && p_mat(ab,w) < 0.05
        text(ax8, w, ab, sprintf('%.2f\n%s', rho_mat(ab,w), sig_sym(p_mat(ab,w))), ...
            'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color','w');
    else
        text(ax8, w, ab, sprintf('%.2f', rho_mat(ab,w)), ...
            'HorizontalAlignment','center','FontSize',8,'Color','k');
    end
end; end

exportgraphics(figD, fullfile(dir_out,'FigS_D_Heatmaps.png'), 'Resolution',DPI);
close(figD);

%% ════════════════════════════════════════════════════════════════════════════
%  FIG E — Rayleigh: preferred phase by band × window
%% ════════════════════════════════════════════════════════════════════════════
fprintf('Generating FigS_E (Rayleigh preferred phases)...\n');

sz_p       = size(pref_ch);
n_wins_rayl = 3;
band_tex_names = {'Low-\alpha (8-13 Hz)','High-\alpha (14-20 Hz)','\beta (20-30 Hz)', ...
                  'Low-\gamma (30-50 Hz)','High-\gamma (50-80 Hz)','Broad-\gamma (30-100 Hz)'};

figE = figure('Color','w','Units','inches','Position',[0 0 10 12]);
tl5  = tiledlayout(figE, n_ab, n_wins_rayl, 'TileSpacing','compact','Padding','compact');
title(tl5, 'Phase alignment to frontocentral \theta-valley', ...
    'FontWeight','bold','FontSize',14);

for ab = 1:n_ab
    fi = amp_freqs >= amp_bands(ab).range(1) & amp_freqs <= amp_bands(ab).range(2);
    for w = 1:n_wins_rayl
        ax = nexttile;
        angles = extract_pref_angles(pref_ch, sz_p, N, fi, w, numel(amp_freqs));
        if isempty(angles) || all(isnan(angles))
            text(0.5, 0.5, 'N/A', 'Units','norm','HorizontalAlignment','center');
            axis(ax,'off'); continue;
        end
        angles_clean = angles(~isnan(angles));
        [p_ray, R_ray, mu_ray] = rayleigh_test(angles_clean);

        n_bins  = 16;
        edges   = linspace(-pi, pi, n_bins+1);
        cnts    = histcounts(angles_clean, edges);
        hold(ax,'on');
        for bi = 1:n_bins
            if cnts(bi) == 0; continue; end
            theta_fill = linspace(edges(bi), edges(bi+1), 20);
            r_fill     = (cnts(bi) / numel(angles_clean)) * 1.5;
            fill(ax, [0, r_fill*cos(theta_fill), 0], [0, r_fill*sin(theta_fill), 0], ...
                BAND_CLRS(ab,:), 'FaceAlpha',0.85, 'EdgeColor','none');
        end

        th_circ = linspace(-pi, pi, 100);
        plot(ax, 0.5*cos(th_circ), 0.5*sin(th_circ), '-', 'Color',[0.85 0.85 0.85],'LineWidth',0.8);
        plot(ax, [-0.5 0.5], [0 0], ':','Color',[0.8 0.8 0.8],'LineWidth',0.8);
        plot(ax, [0 0], [-0.5 0.5], ':','Color',[0.8 0.8 0.8],'LineWidth',0.8);

        r_mean_plot = min(R_ray, 0.5);
        quiver(ax, 0, 0, r_mean_plot*cos(mu_ray), r_mean_plot*sin(mu_ray), ...
            0, 'Color','k','LineWidth',1.8,'MaxHeadSize',0.8);
        axis(ax,'equal','off');
        xlim(ax, [-0.65 0.65]); ylim(ax, [-0.65 0.65]);

        if ab == 1
            title(ax, wins_lbl{w}, 'FontSize',12,'FontWeight','bold');
        end
        if w == 1
            text(ax, -0.75, 0, band_tex_names{ab}, ...
                'HorizontalAlignment','right','VerticalAlignment','middle', ...
                'FontSize',11,'FontWeight','bold');
        end
        text(ax,  0.55, 0, '0',   'HorizontalAlignment','left',  'FontSize',9,'Color',[0.3 0.3 0.3]);
        text(ax, -0.55, 0, '\pi', 'HorizontalAlignment','right', 'FontSize',9,'Color',[0.3 0.3 0.3]);

        stat_color  = ternary(p_ray < 0.05, 'k', [0.4 0.4 0.4]);
        stat_weight = ternary(p_ray < 0.05, 'bold', 'normal');
        text(ax, 0, -0.65, sprintf('R = %.2f | p = %.3f%s', R_ray, p_ray, ...
            ternary(p_ray<0.05,' *','')), ...
            'HorizontalAlignment','center','FontSize',9, ...
            'FontWeight',stat_weight,'Color',stat_color);
    end
end

exportgraphics(figE, fullfile(dir_out,'FigS_E_Rayleigh_Nature.png'), 'Resolution',DPI);
close(figE);

%% ════════════════════════════════════════════════════════════════════════════
%  RAYLEIGH — 3 stand-alone highlight plots
%% ════════════════════════════════════════════════════════════════════════════
% [band_idx, window_idx]  bands: 1=LowAlpha(8-13), 3=Beta(20-30), 6=BroadGamma
% windows: 1=Early, 2=Late, 3=Mid
stars_to_plot     = [1 3; 2 1; 3 3];
idx_main          = [1, 3, 6];
band_tex_names_main = {'Low-\alpha (8-13 Hz)','\beta (20-30 Hz)','Broad-\gamma (30-100 Hz)'};

for i = 1:size(stars_to_plot,1)
    ab_i     = stars_to_plot(i,1);
    w_i      = stars_to_plot(i,2);
    idx_orig = idx_main(ab_i);
    fi       = amp_freqs >= amp_bands(idx_orig).range(1) & amp_freqs <= amp_bands(idx_orig).range(2);

    figE_ind = figure('Color','w','Units','inches','Position',[0 0 4 4]);
    ax = axes('Parent',figE_ind);

    angles       = extract_pref_angles(pref_ch, sz_p, N, fi, w_i, numel(amp_freqs));
    angles_clean = angles(~isnan(angles));
    [p_ray, R_ray, mu_ray] = rayleigh_test(angles_clean);

    n_bins = 16; edges = linspace(-pi, pi, n_bins+1);
    cnts   = histcounts(angles_clean, edges);
    hold(ax,'on');
    for bi = 1:n_bins
        if cnts(bi) == 0; continue; end
        theta_fill = linspace(edges(bi), edges(bi+1), 20);
        r_fill     = (cnts(bi)/numel(angles_clean)) * 1.5;
        fill(ax, [0, r_fill*cos(theta_fill), 0], [0, r_fill*sin(theta_fill), 0], ...
            BAND_CLRS(idx_orig,:), 'FaceAlpha',0.8,'EdgeColor','none');
    end
    th_circ = linspace(-pi, pi, 100);
    plot(ax, 0.5*cos(th_circ), 0.5*sin(th_circ), '-','Color',[0.7 0.7 0.7],'LineWidth',0.8);
    plot(ax, [-0.5 0.5],[0 0], '-','Color',[0.8 0.8 0.8],'LineWidth',0.8);
    plot(ax, [0 0],[-0.5 0.5], '-','Color',[0.8 0.8 0.8],'LineWidth',0.8);

    r_mean_plot = min(R_ray, 0.48);
    quiver(ax, 0, 0, r_mean_plot*cos(mu_ray), r_mean_plot*sin(mu_ray), ...
        0, 'Color','k','LineWidth',2.5,'MaxHeadSize',0.8);
    axis(ax,'equal','off');
    xlim(ax,[-0.7 0.7]); ylim(ax,[-0.6 0.6]);

    if     p_ray < 0.001; sig_txt = '***';
    elseif p_ray < 0.01;  sig_txt = '**';
    elseif p_ray < 0.05;  sig_txt = '*';
    else;                 sig_txt = 'n.s.';
    end
    text(ax, 0, 0.55, sig_txt, 'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
    title(ax, sprintf('%s (%s)', band_tex_names_main{ab_i}, wins_lbl{w_i}), ...
        'FontSize',12,'FontWeight','bold');
    text(ax,  0.55, 0, '0',   'HorizontalAlignment','left',  'FontSize',10,'Color',[0.3 0.3 0.3]);
    text(ax, -0.55, 0, '\pi', 'HorizontalAlignment','right', 'FontSize',10,'Color',[0.3 0.3 0.3]);

    fname = sprintf('Fig2_Rayleigh_%s_%s.png', amp_bands(idx_orig).name, wins_lbl{w_i});
    exportgraphics(figE_ind, fullfile(dir_out,fname), 'Resolution',300);
    close(figE_ind);
end
end   % if have_eeg  (FigS_D + FigS_E + Rayleigh individuales requieren EEG-EEG PAC)

%% ── TEXT REPORT ──────────────────────────────────────────────────────────────
f_rep = fullfile(dir_out, 'Reporte_Supplementary.txt');
fid   = fopen(f_rep,'w');
W_rep = @(varargin) fprintf(fid, varargin{:});

W_rep('================================================================\n');
W_rep('  SUPPLEMENTARY ANALYSES — Report (S_Supplementary_Plots.m)\n');
W_rep('  Generated: %s\n', datestr(now));
W_rep('  N = %d Cases | ROI: F1+FC1\n', N);
W_rep('================================================================\n\n');

% A: Steiger
W_rep('── A: STEIGER Z — formal double dissociation ───────────────────\n');
W_rep('  r12 = rho(alpha_Mid,  IES_Ch)    = %+.4f   p = %.4f  %s\n', r12, p12, sig_sym(p12));
W_rep('  r13 = rho(beta_PAC_Ear, IES_Ch)    = %+.4f   p = %.4f  %s\n', r13, p13, sig_sym(p13));
W_rep('  r23 = rho(alpha_Mid,  beta_PAC)  = %+.4f   p = %.4f  %s\n', r23, p23, sig_sym(p23));
W_rep('  N   = %d\n', N);
W_rep('  Z   = %+.4f   p = %.4f  %s\n', Z_st, p_st, sig_sym(p_st));
W_rep('  r23 interp: %+.3f (%s) — predictors covary; despite shared variance,\n', r23, sig_sym(p23));
W_rep('  directional effects on IES are formally opposite (Steiger **).\n\n');

% B: Mediation
W_rep('── B: BOOTSTRAPPED MEDIATION (%d iterations) ────────────────────\n', NBOOT);
W_rep('  Model: beta_PAC_Early (X) -> alpha_Mid (M) -> IES_Ch (Y)\n');
W_rep('  a  (X->M):          beta = %+.4f\n', a);
W_rep('  b  (M->Y|X):        beta = %+.4f\n', b);
W_rep('  c  (total X->Y):    beta = %+.4f\n', c);
W_rep('  c'' (direct X->Y):  beta = %+.4f\n', cp);
W_rep('  ab (indirect):      beta = %+.4f   95%% CI [%.4f, %.4f]\n', ab, ci_lo, ci_hi);
if sig_boot
    W_rep('  -> SIGNIFICANT MEDIATION (CI does not cross zero)\n\n');
else
    W_rep('  -> Mediation NOT significant (CI crosses zero)\n\n');
end

% C: ITPC
W_rep('── C: ITPC alpha — predictor of IES ───────────────────────────\n');
W_rep('  rho(alpha_ITPC, IES_Ch)                   = %+.4f   p = %.4f  %s\n', r_i, p_i, sig_sym(p_i));
W_rep('  rho_partial(alpha_ITPC, IES | alpha_power) = %+.4f\n', rp);
W_rep('  t(%d) = %.4f   p = %.4f  %s\n', df_p, t_p, p_p, sig_sym(p_p));
W_rep('  -> Null result: mechanism is power suppression, not phase.\n\n');

if have_eeg
% D: EEG-PAC heatmaps
W_rep('── D: EEG-EEG PAC — zMI_Ch vs IES_Ch by band x window ─────────\n');
W_rep('  %-12s  %-8s  %-8s  %-8s  %-8s  %-8s  %-8s\n', ...
      'Band','med_E','med_L','med_M','rho_E','rho_L','rho_M');
W_rep('  %s\n', repmat('-',1,68));
for ab = 1:n_ab
    fi2      = amp_freqs >= amp_bands(ab).range(1) & amp_freqs <= amp_bands(ab).range(2);
    row_med  = nan(1,3); row_rho = nan(1,3); row_p = nan(1,3);
    for w = 1:3
        zw = squeeze(mean(zMI_eeg(:,fi2,w),2,'omitnan'));
        row_med(w) = median(zw,'omitnan');
        vv = ~isnan(zw) & ~isnan(ies_eeg);
        if sum(vv) >= 5
            [row_rho(w), row_p(w)] = corr(double(zw(vv)), ies_eeg(vv), 'Type','Spearman');
        end
    end
    W_rep('  %-12s  %+.3f   %+.3f   %+.3f   %+.3f%s  %+.3f%s  %+.3f%s\n', ...
        amp_bands(ab).name, row_med(1), row_med(2), row_med(3), ...
        row_rho(1), ternary(~isnan(row_p(1))&&row_p(1)<0.05,'*',' '), ...
        row_rho(2), ternary(~isnan(row_p(2))&&row_p(2)<0.05,'*',' '), ...
        row_rho(3), ternary(~isnan(row_p(3))&&row_p(3)<0.05,'*',' '));
end
W_rep('\n');

% E: Rayleigh
W_rep('── E: RAYLEIGH — preferred phase theta x band x window ─────────\n');
W_rep('  %-12s  %-22s  %-22s  %-22s\n','Band','Early','Late','Mid');
W_rep('  %s\n', repmat('-',1,82));
for ab = 1:n_ab
    fi2 = amp_freqs >= amp_bands(ab).range(1) & amp_freqs <= amp_bands(ab).range(2);
    row = sprintf('  %-12s', amp_bands(ab).name);
    for w = 1:3
        ang = extract_pref_angles(pref_ch, size(pref_ch), N, fi2, w, numel(amp_freqs));
        ang = ang(~isnan(ang));
        if numel(ang) >= 3
            [p_r, R_r, mu_r] = rayleigh_test(ang);
            row = [row sprintf('  R=%.3f p=%.3f%-3s mu=%.2f', R_r, p_r, ...
                ternary(p_r<0.05,' * ','   '), mu_r)]; %#ok<AGROW>
        else
            row = [row sprintf('  %-22s','N/A')]; %#ok<AGROW>
        end
    end
    W_rep('%s\n', row);
end
W_rep('\n  mu in radians (0 = theta peak, pi = theta trough)\n');
W_rep('  Interpretation: Early and Mid -> consistent preferred phase (~pi)\n');
W_rep('  Late -> uniform distribution (n.s. across all bands)\n');
W_rep('  Mechanism: theta organises broadband excitability, not selective gamma nesting\n\n');
end   % if have_eeg (reporte D/E)
W_rep('================================================================\n  END OF REPORT\n================================================================\n');
fclose(fid);

%% ── SUMMARY ──────────────────────────────────────────────────────────────────
fprintf('\nOutputs saved to:\n  %s\n\n', dir_out);
fprintf('  FigS_A_Alpha_IES.png\n');
fprintf('  FigS_A_BetaPAC_IES.png\n');
fprintf('  FigS_B_Mediacion.png\n');
fprintf('  FigS_C_ITPC.png\n');
fprintf('  FigS_D_Heatmaps.png\n');
fprintf('  FigS_E_Rayleigh_Nature.png\n');
fprintf('  Reporte_Supplementary.txt\n\n');

%% ════════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════════
function [p, R, mu] = rayleigh_test(angles)
    angles = angles(~isnan(angles(:)));
    n = numel(angles);
    if n < 3; p=1; R=0; mu=0; return; end
    C  = mean(cos(angles));
    S  = mean(sin(angles));
    R  = sqrt(C^2 + S^2);
    z  = n * R^2;
    p  = exp(-z) * (1 + (2*z - z^2)/(4*n) ...
         - (24*z - 132*z^2 + 76*z^3 - 9*z^4)/(288*n^2));
    p  = max(0, min(1, real(p)));
    mu = atan2(S, C);
end

function angles = extract_pref_angles(pref_ch, sz_p, N, fi, w, n_freqs)
    angles = nan(N,1);
    if numel(sz_p) == 3
        if sz_p(1)==N && sz_p(2)==n_freqs && sz_p(3)>=w
            slice = squeeze(pref_ch(:, fi, w));
            if isvector(slice)
                angles = double(slice(:));
            else
                angles = circ_mean_mat(double(slice));
            end
        elseif sz_p(3)==N && sz_p(1)==n_freqs && sz_p(2)>=w
            slice = squeeze(pref_ch(fi, w, :));
            if isvector(slice)
                angles = double(slice(:));
            else
                angles = circ_mean_mat(double(slice'));
            end
        else
            try
                angles = double(squeeze(mean(pref_ch(:, fi, w), 2,'omitnan')));
            catch
                angles = nan(N,1);
            end
        end
    elseif numel(sz_p) == 2
        if sz_p(1)==N && sz_p(2)>=w
            angles = double(pref_ch(:,w));
        elseif sz_p(2)==N && sz_p(1)>=w
            angles = double(pref_ch(w,:)');
        end
    elseif numel(pref_ch) == N
        angles = double(pref_ch(:));
    end
end

function mu = circ_mean_mat(A)
    mu = atan2(mean(sin(A),2,'omitnan'), mean(cos(A),2,'omitnan'));
end

function [Z, p] = steiger_test(r12, r13, r23, n)
    z12   = atanh(r12); z13 = atanh(r13);
    r_bar = (r12^2 + r13^2) / 2;
    f     = (1-r23) / (2*(1-r_bar));
    h     = (1-f*r_bar) / (1-r_bar);
    Z     = (z12-z13) * sqrt((n-3) / (2*(1-r23)*h));
    p     = 2*(1-normcdf(abs(Z)));
end

function b = ols_slope(X,Y)
    X1 = [ones(numel(X),1) X(:)]; bv = X1\Y(:); b = bv(2);
end

function b = ols_slopes(Xm,Y)
    X1 = [ones(size(Xm,1),1) Xm]; bv = X1\Y(:); b = bv(2:end);
end

function rp = partial_spearman(X,Y,Z)
    rx = tiedrank(X(:)); rx = rx-mean(rx);
    ry = tiedrank(Y(:)); ry = ry-mean(ry);
    rz = tiedrank(Z(:)); rz = rz-mean(rz);
    rxy = dot(rx,ry)/(norm(rx)*norm(ry));
    rxz = dot(rx,rz)/(norm(rx)*norm(rz));
    ryz = dot(ry,rz)/(norm(ry)*norm(rz));
    rp  = (rxy-rxz*ryz) / (sqrt(max(0,1-rxz^2))*sqrt(max(0,1-ryz^2)));
end

function s = sig_sym(p)
    if     p < 0.001; s = '***';
    elseif p < 0.01;  s = '**';
    elseif p < 0.05;  s = '*';
    elseif p < 0.10;  s = '(+)';
    else;             s = 'NS';
    end
end

function s = ternary(c,a,b)
    if c; s=a; else; s=b; end
end

function v = ternary_val(c,a,b)
    if c; v=a; else; v=b; end
end

function cmap = div_cmap(N,cold,hot)
    half = ceil(N/2); w = [1 1 1];
    c1   = interp1([0;1],[cold;w],linspace(0,1,half)');
    c2   = interp1([0;1],[w;hot], linspace(0,1,half)');
    cmap = [c1; c2(2:end,:)];
end

function annotate_arrow(ax,x,y,lbl,clr)
    dx = x(2)-x(1); dy = y(2)-y(1);
    quiver(ax,x(1),y(1),dx*0.85,dy*0.85,0,'Color',clr,'LineWidth',1.8,'MaxHeadSize',0.35);
    text(ax,mean(x),mean(y)+0.22,lbl,'FontSize',9,'Color',clr,'FontWeight','bold', ...
        'HorizontalAlignment','center');
end
