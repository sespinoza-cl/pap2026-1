%% ============================================================================
%  P_Rayleigh_Panels.m — Polar plots de fase preferida (acoplamiento masticatorio)
%  ----------------------------------------------------------------------------
%  REEMPLAZA el viejo Rayleigh θ-γ broadband (cortado: γ artefactual + bandas
%  inconsistentes). Aquí: fase preferida de la amplitud EEG respecto al ciclo
%  masticatorio (EMG), en las BANDAS CANÓNICAS θ/α/β × ventanas Early/Late/Active,
%  condición Cases/Chew. Estilo rosa con barras de distribución + vector resultante.
%
%  LEE: outputs/Figure04_PAC/PAC_4Groups_Workspace.mat  (R.Casos.Ch.pref_phase)
%  ESCRIBE: outputs/figures/P_Rayleigh_<Band>_<Win>.png   (9 paneles individuales)
%  ============================================================================
clear; clc; close all;

P = S0_paths();
DIR_FIG = fullfile(P.dir_out, 'figures');
if ~exist(DIR_FIG,'dir'); mkdir(DIR_FIG); end
DPI = 300;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',10,'DefaultFigureColor','w');

% Paleta: hue por banda (θ verde, α azul, β violeta); tono por ventana temporal
band_base = [0.16 0.52 0.30;    % theta  — verde
             0.13 0.42 0.70;    % alpha  — azul
             0.50 0.20 0.55];   % beta   — violeta

%% ── Cargar workspace PAC ────────────────────────────────────────────────────
f_ws = fullfile(P.fig04,'PAC_4Groups_Workspace.mat');
assert(exist(f_ws,'file')==2,'Falta %s — corre S06_PAC.m primero', f_ws);
W = load(f_ws);
pref = W.R.Casos.Ch.pref_phase;          % [nS × nF × nW]
amp_freqs = W.P.amp_freqs;               % 4:2:60
BANDS = P.bands; BHZ = P.bands_hz; WINS = P.wins;   % θ/α/β ; Early/Late/Active
band_nm = {'\theta','\alpha','\beta'};

%% ── Un panel polar por banda × ventana ─────────────────────────────────────
for bi = 1:numel(BANDS)
    bmask = amp_freqs >= BHZ{bi}(1) & amp_freqs <= BHZ{bi}(2);
    for wi = 1:numel(WINS)
        phi = squeeze(mean(pref(:, bmask, wi), 2, 'omitnan'));   % fase por sujeto
        phi = phi(~isnan(phi));
        ray = rayleigh(phi);

        % Color: base de la banda + tono según ventana (cronológico: Early claro → Late oscuro)
        switch WINS{wi}
            case 'Early', sf =  0.50;   % más claro
            case 'Mid',   sf =  0.12;
            case 'Late',  sf = -0.28;   % más oscuro
            otherwise,    sf =  0;
        end
        clr = shade_color(band_base(bi,:), sf);

        fig = figure('Position',[60 60 420 460],'Visible','off');
        ax  = axes(fig); hold(ax,'on'); axis(ax,'equal','off');

        nb = 12; edges = linspace(-pi,pi,nb+1);
        counts = histcounts(phi, edges);
        r_max  = max(counts)*1.25 + 0.5;

        % Sectores (barras de distribución)
        for b = 1:nb
            if counts(b)==0; continue; end
            th = linspace(edges(b),edges(b+1),20);
            patch(ax, [0 counts(b)*sin(th) 0], [0 counts(b)*cos(th) 0], clr, ...
                  'EdgeColor','w','FaceAlpha',0.85,'LineWidth',0.5);
        end
        % Círculo de referencia + ejes
        tt = linspace(0,2*pi,100);
        plot(ax, r_max*sin(tt), r_max*cos(tt), '-','Color',[.75 .75 .75],'LineWidth',0.8);
        plot(ax, [-r_max r_max],[0 0],'-','Color',[.85 .85 .85]);
        plot(ax, [0 0],[-r_max r_max],'-','Color',[.85 .85 .85]);
        % Vector resultante (longitud ∝ R)
        quiver(ax, 0,0, r_max*0.9*ray.R*sin(ray.mu), r_max*0.9*ray.R*cos(ray.mu), 0, ...
               'Color','k','LineWidth',2.4,'MaxHeadSize',0.5);
        % Etiquetas de fase (0 = cresta masticatoria, π = valle)
        text(ax, 0, r_max*1.06, '0','HorizontalAlignment','center','FontSize',9);
        text(ax,-r_max*1.10, 0, '\pi','HorizontalAlignment','center','FontSize',10);

        sig = ''; if ray.p<0.001, sig='***'; elseif ray.p<0.01, sig='**'; elseif ray.p<0.05, sig='*'; end
        title(ax, sprintf('%s — %s\nR=%.3f  p=%.3g %s', band_nm{bi}, WINS{wi}, ray.R, ray.p, sig), ...
              'FontSize',11,'Interpreter','tex');

        fn = sprintf('P_Rayleigh_%s_%s.png', BANDS{bi}, WINS{wi});
        exportgraphics(fig, fullfile(DIR_FIG,fn),'Resolution',DPI);
        close(fig);
        fprintf('-> %s   (R=%.3f, p=%.3g, n=%d)\n', fn, ray.R, ray.p, numel(phi));
    end
end
fprintf('\n✓ P_Rayleigh_Panels.m listo. 9 paneles en: %s\n', DIR_FIG);

%% ── Rayleigh test ──────────────────────────────────────────────────────────
function ray = rayleigh(phi)
    phi = phi(:); n = numel(phi);
    R  = abs(mean(exp(1i*phi)));
    mu = angle(mean(exp(1i*phi)));
    Z  = n*R^2;
    p  = exp(-Z)*(1 + (2*Z-Z^2)/(4*n) - (24*Z-132*Z^2+76*Z^3-9*Z^4)/(288*n^2));
    p  = max(0,min(1,p));
    ray = struct('R',R,'mu',mu,'p',p,'n',n);
end

function c = shade_color(base, f)
    % f>0 aclara (hacia blanco); f<0 oscurece (hacia negro)
    if f >= 0
        c = base + (1-base)*f;
    else
        c = base * (1+f);
    end
    c = min(max(c,0),1);
end
