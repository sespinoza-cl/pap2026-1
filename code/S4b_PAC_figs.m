% S4b_PAC_figs.m  —  Figuras PAC para el paper (P2V1)
%
% Genera 2 figuras desde v1_S4_PAC_stats.mat:
%   Fig 1: S4b_PAC_zMI_theta.png
%           zMI theta vs null (boxplot+scatter), n=15/31 sig, p<0.0001
%   Fig 2: S4b_PAC_theta_specificity.png  [2-panel]
%           Panel A: Delta MI (Late-Base) x 3 bandas
%           Panel B: MI_late absoluto x 3 bandas (Friedman)
%
% Fuente: Analysis_V1_Final/data/computed/v1_S4_PAC_stats.mat
% Guardar en: Analysis_V1_Final/outputs/figures_ok/
%             Paper_plots/Final_Figures/

run(fullfile(fileparts(mfilename('fullpath')), 'S0_config.m'));

%% ── Cargar datos ─────────────────────────────────────────────────────────
mat_path = fullfile(DATA_FOR_PLOTS, 'v1_S4_PAC_stats.mat');
if ~exist(mat_path,'file')
    error('No encontrado: %s', mat_path);
end
S = load(mat_path);

zMI    = S.zMI_all(:);          % 31x1
MI_late = S.MI_late;            % 3x31
MI_base = S.MI_base;            % 3x31
MI_early= S.MI_early;           % 3x31
delta   = MI_late - MI_base;    % 3x31

bands   = {'θ (4–7 Hz)','α (8–12 Hz)','β (13–20 Hz)'};
clrs    = {[0.84 0.15 0.16], [0.12 0.47 0.71], [0.17 0.63 0.17]};

% Directorios de salida
DIR_FIG  = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                    'outputs','figures_ok');
DIR_PAPER = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
                     'Paper_plots','Final_Figures');

%% ── Stats ────────────────────────────────────────────────────────────────
% zMI vs 0
[~, p_zmi]   = signrank(zMI);
k_sig        = sum(zMI > 1.96);
p_binom      = 1 - binocdf(k_sig-1, sum(~isnan(zMI)), 0.05);

% Delta MI por banda (Wilcoxon vs 0)
p_delta = nan(3,1);
for b = 1:3
    [~, p_delta(b)] = signrank(delta(b,:)');
end

% Friedman sobre MI_late
[p_fried, ~] = friedman(MI_late', 1, 'off');

% Post-hoc Wilcoxon pairwise (sin corrección, para indicar)
pairs = [1 2; 1 3; 2 3];
p_ph  = nan(3,1);
for k = 1:3
    [~, p_ph(k)] = ranksum(MI_late(pairs(k,1),:)', MI_late(pairs(k,2),:)');
end
% p_ph(1): theta vs alpha = 0.165
% p_ph(2): theta vs beta  = 0.010
% p_ph(3): alpha vs beta  = 0.047

fprintf('\n=== PAC Stats ===\n');
fprintf('zMI vs 0: Wilcoxon p=%.4e  | %d/31 sig | binom p=%.4e\n', ...
        p_zmi, k_sig, p_binom);
fprintf('Delta MI (Late-Base) Wilcoxon p: theta=%.3f  alpha=%.3f  beta=%.3f\n', ...
        p_delta(1), p_delta(2), p_delta(3));
fprintf('MI_late Friedman p=%.3f\n', p_fried);
fprintf('Post-hoc theta/alpha=%.3f  theta/beta=%.3f  alpha/beta=%.3f\n', ...
        p_ph(1), p_ph(2), p_ph(3));

%% ════════════════════════════════════════════════════════════════════════
%  FIGURA 1: zMI theta vs null
%  ════════════════════════════════════════════════════════════════════════
fig1 = figure('Units','centimeters','Position',[2 2 10 12],...
              'Color','w','Renderer','painters');

ax = axes(fig1);
hold(ax,'on');

YMAX = 22;
z_plot = min(zMI, YMAX);   % cap para visualización
n_clip = sum(zMI > YMAX);

% Scatter con jitter
rng(42);
jit = (rand(numel(z_plot),1)-0.5)*0.24;
scatter(ax, 1+jit, z_plot, 28, [0.84 0.15 0.16], 'filled',...
        'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor','none');

% Boxplot manual
draw_box(ax, z_plot, 1, 0.30, [0.84 0.15 0.16], 0.35);

% Líneas de referencia
yline(ax, 0,    '--k', 'LineWidth', 0.8, 'Alpha', 0.5);
yline(ax, 1.96, ':',  'Color',[0.5 0.5 0.5], 'LineWidth', 0.8);
text(ax, 1.22, 1.96, 'z = 1.96', 'VerticalAlignment','middle',...
     'FontSize', 8.5, 'Color', [0.55 0.55 0.55]);

% Stats text
text(ax, 0.97, 0.97, sprintf('vs null: p < 0.0001\n%d/31 sig. (binom. p < 0.0001)', k_sig),...
     'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
     'FontSize', 8.5, 'Color', [0.2 0.2 0.2]);

% Nota outliers
if n_clip > 0
    text(ax, 1.15, YMAX-0.5, sprintf('(%d subjects z>%d)', n_clip, YMAX),...
         'FontSize', 7.5, 'Color', [0.6 0.6 0.6], 'FontAngle','italic');
end

ylim(ax, [-5 YMAX+1.5]);
xlim(ax, [0.55 1.55]);
set(ax,'XTick',1,'XTickLabel',{'\theta-PAC  (4–7 Hz)'},'TickDir','out',...
       'Box','off','FontSize',11,'XAxis.Visible','off');
ylabel(ax, 'z-score (MI vs. shuffled null)','FontSize',11);
ax.XAxis.Visible = 'off';
ax.Box = 'off';

print_fig(fig1, DIR_FIG,   'S4b_PAC_zMI_theta');
print_fig(fig1, DIR_PAPER, 'S4b_PAC_zMI_theta');
fprintf('Fig1 guardada.\n');

%% ════════════════════════════════════════════════════════════════════════
%  FIGURA 2: 2-panel theta specificity
%  ════════════════════════════════════════════════════════════════════════
fig2 = figure('Units','centimeters','Position',[2 16 22 11],...
              'Color','w','Renderer','painters');

xs = [0 1 2];
bw = 0.38;

% ── Panel A: Delta MI (Late-Base) x banda ────────────────────────────────
ax1 = subplot(1,2,1,'Parent',fig2);
hold(ax1,'on');

YMIN_A = -148; YMAX_A = 130;
delta5  = delta * 1e5;   % unidades x10^-5

for b = 1:3
    d_v = delta5(b,:)';
    d_c = max(min(d_v, YMAX_A), YMIN_A);   % clip
    rng(b);
    jit = (rand(numel(d_c),1)-0.5)*0.18;
    scatter(ax1, xs(b)+jit, d_c, 18, clrs{b}, 'filled',...
            'MarkerFaceAlpha',0.22,'MarkerEdgeColor','none');
    draw_box(ax1, d_c, xs(b), bw, clrs{b}, 0.40);
end

yline(ax1, 0, '--k', 'LineWidth', 0.8, 'Alpha', 0.50);

% Marcadores de significancia (fila superior fija)
for b = 1:3
    lbl = p2s(p_delta(b));
    bold_w = 'normal';
    if contains(lbl,'*'), bold_w='bold'; end
    if p_delta(b) < 0.05
        % Beta: significativo, con flecha apuntando abajo
        text(ax1, xs(b), YMAX_A-10, lbl, 'HorizontalAlignment','center',...
             'VerticalAlignment','top','FontSize',11,'FontWeight','bold',...
             'Color',clrs{b});
        annotation(fig2,'arrow',...
            ax2fig(ax1,[xs(b) xs(b)],[YMAX_A-18 mean(delta5(b,:))]),...
            'Color',clrs{b},'LineWidth',1.2,'HeadStyle','vback3','HeadLength',6,'HeadWidth',5);
    else
        text(ax1, xs(b), YMAX_A-8, lbl, 'HorizontalAlignment','center',...
             'VerticalAlignment','top','FontSize',9.5,'FontWeight',bold_w,...
             'Color',clrs{b});
    end
end

n_cl = sum(delta5(3,:) < YMIN_A);
if n_cl>0
    text(ax1, xs(3), YMIN_A+3, sprintf('(%d outlier no mostrado)',n_cl),...
         'HorizontalAlignment','center','VerticalAlignment','bottom',...
         'FontSize',7,'Color',[0.6 0.6 0.6],'FontAngle','italic');
end

xlim(ax1,[-0.6 2.6]); ylim(ax1,[YMIN_A YMAX_A]);
set(ax1,'XTick',xs,'XTickLabel',bands,'TickDir','out','Box','off','FontSize',10);
ylabel(ax1,'MI change: Late \minus Base  (\times10^{-5})','FontSize',10);
text(ax1,-0.12,1.04,'A','Units','normalized','FontSize',14,'FontWeight','bold');

% ── Panel B: MI_late absoluto x banda ────────────────────────────────────
ax2 = subplot(1,2,2,'Parent',fig2);
hold(ax2,'on');

mil3 = MI_late * 1e3;   % unidades x10^-3

for b = 1:3
    rng(b+10);
    jit = (rand(31,1)-0.5)*0.18;
    scatter(ax2, xs(b)+jit, mil3(b,:)', 18, clrs{b}, 'filled',...
            'MarkerFaceAlpha',0.22,'MarkerEdgeColor','none');
    draw_box(ax2, mil3(b,:)', xs(b), bw, clrs{b}, 0.40);
end

% Brackets post-hoc
ym = max(prctile(mil3(:),98));
arm_y = 0.005*ym;
% theta vs beta (p=0.010, *)
draw_bracket(ax2, xs(1), xs(3), ym*1.08, arm_y, '*');
% alpha vs beta (p=0.047, *)
draw_bracket(ax2, xs(2), xs(3), ym*1.25, arm_y, '*');
text(ax2, 0.97, 0.97, sprintf('Friedman p = %.3f', p_fried),...
     'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top',...
     'FontSize', 8.5, 'Color', [0.35 0.35 0.35]);

xlim(ax2,[-0.6 2.6]); ylim(ax2,[0 ym*1.44]);
set(ax2,'XTick',xs,'XTickLabel',bands,'TickDir','out','Box','off','FontSize',10);
ylabel(ax2,'MI in WIN\_LATE  (\times10^{-3})','FontSize',10);
text(ax2,-0.12,1.04,'B','Units','normalized','FontSize',14,'FontWeight','bold');

print_fig(fig2, DIR_FIG,   'S4b_PAC_theta_specificity');
print_fig(fig2, DIR_PAPER, 'S4b_PAC_theta_specificity');
fprintf('Fig2 guardada.\n');

%% ── Funciones locales ────────────────────────────────────────────────────
function draw_box(ax, data, xc, bw, clr, alp)
    if nargin<6, alp=0.40; end
    v  = sort(data(~isnan(data)));
    q1 = prctile(v,25); q2 = prctile(v,50); q3 = prctile(v,75);
    iqr_v = q3-q1;
    low_f = q1-1.5*iqr_v; hi_f = q3+1.5*iqr_v;
    wlo = min(v(v>=low_f)); if isempty(wlo), wlo=q1; end
    whi = max(v(v<=hi_f));  if isempty(whi), whi=q3; end
    out_v= v(v<low_f | v>hi_f);
    % Box
    fill(ax, xc+[-bw/2 bw/2 bw/2 -bw/2 -bw/2], ...
             [q1 q1 q3 q3 q1], clr, 'FaceAlpha',alp,'EdgeColor',clr,...
             'LineWidth',1.6,'HandleVisibility','off');
    % Median
    plot(ax, xc+[-bw/2 bw/2], [q2 q2], 'Color',clr,'LineWidth',2.2,...
         'HandleVisibility','off');
    % Whiskers
    plot(ax, [xc xc], [wlo q1], 'Color',clr,'LineWidth',1.1,'HandleVisibility','off');
    plot(ax, [xc xc], [q3 whi], 'Color',clr,'LineWidth',1.1,'HandleVisibility','off');
    % Outliers
    if ~isempty(out_v)
        scatter(ax, repmat(xc,numel(out_v),1), out_v, 22, clr,...
                'filled','MarkerFaceAlpha',0.65,'HandleVisibility','off');
    end
end

function draw_bracket(ax, x1, x2, y, arm, lbl)
    plot(ax, [x1 x1 x2 x2], [y-arm y y y-arm], 'k-','LineWidth',0.9,...
         'HandleVisibility','off');
    text(ax, (x1+x2)/2, y+arm*0.5, lbl, 'HorizontalAlignment','center',...
         'VerticalAlignment','bottom','FontSize',10,'FontWeight','bold',...
         'HandleVisibility','off');
end

function s = p2s(p)
    if p < 0.001,      s = '***';
    elseif p < 0.01,   s = '**';
    elseif p < 0.05,   s = '*';
    else,              s = 'n.s.';
    end
end

function xy_fig = ax2fig(ax, xd, yd)
    % Convierte coordenadas de datos a fracciones de figura
    axpos = ax.Position;
    xl = xlim(ax); yl = ylim(ax);
    xf = axpos(1) + axpos(3)*(xd-xl(1))/(xl(2)-xl(1));
    yf = axpos(2) + axpos(4)*(yd-yl(1))/(yl(2)-yl(1));
    xy_fig = [xf; yf];
end

function print_fig(fig, dir_out, fname)
    if ~exist(dir_out,'dir'), mkdir(dir_out); end
    fpath = fullfile(dir_out, [fname '.png']);
    exportgraphics(fig, fpath, 'Resolution', 300, 'BackgroundColor','white');
    fprintf('  → %s\n', fpath);
end
