%% ============================================================================
%  P_PAC_Panels.m — Paneles individuales de PAC (Figura 3 del paper)
%  ----------------------------------------------------------------------------
%  LEE el workspace ya calculado (no recalcula nada): rápido y re-ejecutable.
%    Entrada : outputs/Figure04_PAC/PAC_4Groups_Workspace.mat  (de S06_PAC.m)
%              data/data_beh_tb_45.mat                          (IES, ID-matched)
%    Salida  : outputs/figures/  (PNG individuales para ensamblar en Inkscape)
%
%  Paneles:
%    P_Fig3A_PAC_Heatmap.png        — ΔzMI (Casos Ch-Nc) banda×ventana + signif.
%    P_Fig3B_BetaEarly_vs_IES.png   — scatter β-PAC Early vs IES (Spearman)
%    P_Fig3D_Violin_BetaActive.png  — disociación β-Active en las 4 celdas Grupo×Cond
%
%  (Los scatters β-PAC×IES y α×IES "Nature" también los genera S07 como FigS_A_*.)
%  ============================================================================
clear; clc; close all;

P = S0_paths();
DIR_FIG = fullfile(P.dir_out, 'figures');
if ~exist(DIR_FIG,'dir'); mkdir(DIR_FIG); end
DPI = 300;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',11,'DefaultFigureColor','w');

c_case = [0.80 0.22 0.22];   c_ctrl = [0.22 0.45 0.72];
c_beta = [0.40 0.15 0.55];

%% ── Cargar workspace PAC + conducta (ID-matched) ───────────────────────────
f_ws = fullfile(P.fig04, 'PAC_4Groups_Workspace.mat');
assert(exist(f_ws,'file')==2, 'Falta %s — corre S06_PAC.m primero', f_ws);
W = load(f_ws);
B = W.B;  casos = W.casos;
BANDS = P.bands;  WINS = P.wins;          % {'Theta','Alpha','Beta'} ; {'Early','Late','Active'}
nB = numel(BANDS); nW = numel(WINS);

beh   = load(P.file_beh, 'tb_data_45');
BEH   = beh.tb_data_45.casos.chew;
IES   = nan(numel(casos),1);
for i = 1:numel(casos)
    ix = find(strcmp(string(BEH.Participantes), casos{i}), 1);
    if ~isempty(ix); IES(i) = double(BEH.ies_m(ix)); end
end

%% ── PANEL A: Heatmap ΔzMI (Casos Ch-Nc) ────────────────────────────────────
dz = nan(nB,nW); pm = nan(nB,nW);
for b = 1:nB
    for w = 1:nW
        d = B.Casos.Ch.(BANDS{b})(:,w) - B.Casos.Nc.(BANDS{b})(:,w);
        d = d(~isnan(d));
        dz(b,w) = median(d);
        pm(b,w) = signrank(d);
    end
end

% Workspace cols: [1=Early, 2=Late, 3=Mid] — reorder to chronological display [Early, Mid, Late]
disp_ord  = [1 3 2];                        % Early(1), Mid(3), Late(2)
dz_d  = dz(:, disp_ord);                   % [nB x 3] display order
pm_d  = pm(:, disp_ord);
WINS_D = WINS(disp_ord);                    % {'Early','Mid','Late'}

figA = figure('Position',[60 60 500 520],'Visible','off');
axA  = axes(figA,'Position',[0.20 0.12 0.62 0.76]);
imagesc(axA, dz_d');                        % filas = ventanas display, cols = bandas
cmax = max(abs(dz_d(:)))*1.05;
clim(axA, [-cmax cmax]);
colormap(axA, redblue(256));
cb = colorbar(axA); cb.Label.String = 'Median \DeltazMI (Chew - No-Chew)';
set(axA,'XTick',1:nB,'XTickLabel',{'\theta','\alpha','\beta'}, ...
        'YTick',1:nW,'YTickLabel',WINS_D,'TickLabelInterpreter','tex', ...
        'TickDir','out','Box','on','FontSize',12);
xlabel(axA,'Frequency band','FontWeight','bold');
ylabel(axA,'Time window','FontWeight','bold');
title(axA,'Cases: \DeltazMI (Chew - No-Chew)','FontSize',12);
for b = 1:nB
    for w = 1:nW
        txtc = 'k'; if abs(dz_d(b,w)) > cmax*0.55; txtc = 'w'; end
        text(axA, b, w, sprintf('%.2f\n%s', dz_d(b,w), p2stars(pm_d(b,w))), ...
            'HorizontalAlignment','center','FontSize',11,'FontWeight','bold','Color',txtc);
    end
end
exportgraphics(figA, fullfile(DIR_FIG,'P_Fig3A_PAC_Heatmap.png'),'Resolution',DPI);
close(figA); fprintf('-> P_Fig3A_PAC_Heatmap.png\n');

%% ── PANEL B: β-PAC Early vs IES (Spearman) ─────────────────────────────────
x = B.Casos.Ch.Beta(:,1);   y = IES;
v = ~isnan(x) & ~isnan(y);
[rho, prho] = corr(x(v), y(v), 'Type','Spearman');

figB = figure('Position',[60 60 500 460],'Visible','off');
axB = axes(figB); hold(axB,'on');
scatter(axB, x(v), y(v), 55, c_beta, 'filled','MarkerFaceAlpha',0.75,'MarkerEdgeColor','w');
pf = polyfit(x(v), y(v), 1);
xx = linspace(min(x(v)), max(x(v)), 100);
plot(axB, xx, polyval(pf,xx), '-', 'Color', c_beta, 'LineWidth', 2);
text(axB, 0.04, 0.96, sprintf('\\rho = %+.3f\np_{raw} = %.3f', rho, prho), ...
     'Units','normalized','VerticalAlignment','top','FontSize',12,'FontWeight','bold','Color',c_beta);
set(axB,'Box','off','TickDir','out','YGrid','on','GridAlpha',0.15);
xlabel(axB,'\beta-PAC zMI (Early, Chew)','FontWeight','bold','Interpreter','tex');
ylabel(axB,'IES (ms) - Chew','FontWeight','bold');
title(axB,'\beta-PAC predicts efficiency (uncorrected)','FontSize',12);
exportgraphics(figB, fullfile(DIR_FIG,'P_Fig3B_BetaEarly_vs_IES.png'),'Resolution',DPI);
close(figB); fprintf('-> P_Fig3B_BetaEarly_vs_IES.png  (rho=%+.3f, p_raw=%.3f)\n', rho, prho);

%% ── PANEL D: Violín β-Active en las 4 celdas Grupo×Condición ────────────────
cells = { B.Casos.Ch.Beta(:,3),     B.Casos.Nc.Beta(:,3), ...
          B.Controles.Ch.Beta(:,3), B.Controles.Nc.Beta(:,3) };
labs  = {'Cases/Ch','Cases/Nc','Ctrl/Ch','Ctrl/Nc'};
clrs  = [c_case; c_case*0.6+0.4; c_ctrl; c_ctrl*0.6+0.4];

figD = figure('Position',[60 60 560 460],'Visible','off');
axD = axes(figD); hold(axD,'on');
for k = 1:4
    d = cells{k}; d = d(~isnan(d));
    % violín simple (kernel density) + box
    [f,xi] = ksdensity(d);
    f = f/max(f)*0.35;
    fill(axD, [k+f, fliplr(k-f)], [xi, fliplr(xi)], clrs(k,:), ...
        'FaceAlpha',0.45,'EdgeColor',clrs(k,:)*0.7,'LineWidth',1);
    q = quantile(d,[.25 .5 .75]);
    line(axD,[k-0.12 k+0.12],[q(2) q(2)],'Color','k','LineWidth',2);
    line(axD,[k k],[q(1) q(3)],'Color','k','LineWidth',1.2);
    scatter(axD, k+(rand(numel(d),1)-.5)*0.10, d, 14, clrs(k,:)*0.7, 'filled','MarkerFaceAlpha',0.5);
    % signrank vs 0
    text(axD, k, max(d)+0.4, p2stars(signrank(d)), 'HorizontalAlignment','center','FontWeight','bold','FontSize',13);
end
yline(axD,0,'k--','Color',[.5 .5 .5]);
set(axD,'XTick',1:4,'XTickLabel',labs,'XLim',[0.4 4.6],'Box','off','TickDir','out','YGrid','on','GridAlpha',0.15);
ylabel(axD,'\beta-PAC zMI (Mid)','FontWeight','bold','Interpreter','tex');
title(axD,'Triple dissociation: \beta-PAC only in Cases/Chew','FontSize',12);
exportgraphics(figD, fullfile(DIR_FIG,'P_Fig3D_Violin_BetaMid.png'),'Resolution',DPI);
close(figD); fprintf('-> P_Fig3D_Violin_BetaMid.png\n');

fprintf('\n✓ P_PAC_Panels.m listo. Paneles en: %s\n', DIR_FIG);

%% ── Funciones locales ──────────────────────────────────────────────────────
function s = p2stars(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'n.s.';
    end
end

function cm = redblue(n)
    bot=[0.18 0.42 0.78]; mid=[1 1 1]; top=[0.78 0.15 0.15];
    n2=floor(n/2);
    cm=[interp1([0 1],[bot;mid],linspace(0,1,n2));
        interp1([0 1],[mid;top],linspace(0,1,n-n2))];
end
