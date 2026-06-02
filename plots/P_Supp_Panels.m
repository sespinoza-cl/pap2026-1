%% S_Supp_NewFigures.m
%  ============================================================================
%  Script único para TODAS las figuras suplementarias nuevas del Paper 2.
%  Cada sección (%%SN1 – %%SN7) es independiente y puede correrse por separado
%  siempre que la sección %%CONFIG haya corrido primero.
%
%  FIGURAS GENERADAS  (todas individuales, dimensiones uniformes 560×520 px)
%  ─────────────────────────────────────────────────────────────────────────────
%  S-nuevo 1  FigSN1_FOOOF_DeltaPeriodic.png
%             Δespectro periódico (Chew − NoChew) FOOOF residual, con bandas
%             marcadas y significancia bootstrap por bin de frecuencia.
%             → Contenido único: qué componentes periódicos cambian con masticar.
%
%  S-nuevo 2  FigSN2a_FOOOF_Spectra_NoChew.png
%             FigSN2b_FOOOF_Spectra_Chew.png
%             Espectros periódicos individuales (FOOOF residual), mismo eje Y.
%             FigSN2c_FOOOF_tCF_Distribution.png
%             Histograma de la frecuencia central theta (escalera, ambas condiciones).
%             → El revisor preguntó si los picos en bordes de banda son reales.
%
%  S-nuevo 3  FigSN3a_BetaPAC_SpectralSpecificity.png
%             ΔzMI por banda (θ, α, β) en ventana Early — Casos/Chew.
%             FigSN3b_BetaPAC_TemporalDynamics.png
%             Dinámica temporal del β-PAC (3 ventanas, Ch vs Nc).
%             → Anti-artefacto EMG: solo β significativo.
%
%  S-nuevo 4  FigSN4a_Alpha_ERD_TimeCourse.png
%             Curso temporal de α-ERD (8–13 Hz, ROI frontal).
%             FigSN4b_Alpha_ERD_Topography.png
%             Topograma en ventana activa (si EEGLAB disponible).
%             → "Precision gating" y correlación α × IES.
%
%  S-nuevo 5  FigSN5_Mediation_Path.png
%             Path diagram: β-PAC(Early) → α-ERD(Mid) → IES.
%
%  S-nuevo 6  *** DIAGRAMA MANUAL ***
%             Flujograma de validación del pipeline. Crear en Illustrator.
%
%  S-nuevo 7  FigSN7_FOOOF_Independence.png
%             Scatter Δexponente aperiódico vs Δpotencia theta periódica.
%             → Cierra la pregunta "arousal generalizado".
%
%  Salida: DIR_SUPP (ver sección CONFIG)
%  Sebastián — Mayo 2026
%  ============================================================================

clear; clc; close all; rng(2026);
fprintf('=== S_Supp_NewFigures.m  |  %s ===\n\n', datestr(now));

% ============================================================================
% Script: S_Supp_NewFigures.m
% Edición Especial: Formateo de Alta Resolución Estilo Editorial "Nature"
% ============================================================================
% Todas las figuras se generan con dimensiones físicas idénticas en milímetros
% para garantizar consistencia absoluta al unirlas en un collage final.
% Se eliminan decoraciones redundantes y textos en negrita según guías de la revista.
% ============================================================================

clear; clc; close all; rng(2026);
fprintf('=== S_Supp_NewFigures.m  |  Estética Nature Unificada  |  %s ===\n\n', datestr(now));

% ============================================================================
% Script: S_Supp_NewFigures.m
% Edición Especial: Formateo de Alta Resolución Estilo Editorial "Nature"
% ============================================================================
% Todas las figuras se generan con dimensiones físicas idénticas en centímetros
% (8.9 x 8.2 cm) para garantizar consistencia absoluta al unirlas en un collage final.
% Se eliminan decoraciones redundantes y textos en negrita según guías de la revista.
% ============================================================================

clear; clc; close all; rng(2026);
fprintf('=== S_Supp_NewFigures.m  |  Estética Nature Unificada  |  %s ===\n\n', datestr(now));

%% ── CONFIGURACIÓN (desde S0_paths) ─────────────────────────────────────────
S = S0_paths();

% Archivos de entrada (workspaces ya calculados — referenciados/proyecto)
F_FOOOF      = S.file_fooof_ws;                              % FOOOF_Workspace.mat (referenciado)
F_PAC4G      = fullfile(S.fig04,  'PAC_4Groups_Workspace.mat');  % generado por S06_PAC.m
F_BEH        = S.file_beh;
F_TF_CCH     = fullfile(S.dir_tf, 'Group_Cases_Ch_tfraw.mat');   % TF crudo (referenciado)
F_TF_CNC     = fullfile(S.dir_tf, 'Group_Cases_Nc_tfraw.mat');
F_TF_METRICS = fullfile(S.fig02b, 'TF_band_metrics.mat');        % generado por S02
DIR_EPOCH    = S.dir_epoch;
DIR_EEGLAB   = S.eeglab_path;

% Directorio de salida → paneles individuales del proyecto
DIR_SUPP = fullfile(S.dir_out, 'figures');
if ~exist(DIR_SUPP,'dir'); mkdir(DIR_SUPP); end

% Parámetros de Exportación Estilo Nature (Medidas en CENTÍMETROS para Collage Perfecto)
DPI         = 300;
fig_width   = 8.9;   % Ancho estándar de columna simple en Nature (89 mm)
fig_height  = 8.2;   % Alto estándar unificado para paneles del collage (82 mm)

% Paleta de colores original del manuscrito
c_case      = [0.15 0.55 0.30];   % Verde oscuro — Casos
c_case_fill = [0.60 0.85 0.70];   % Verde claro  — SE Casos
c_ctrl      = [0.18 0.42 0.78];   % Azul oscuro  — Controles
c_ctrl_fill = [0.65 0.80 0.95];   % Azul claro   — SE Controles
c_beta      = [0.40 0.15 0.55];   % Violeta      — Banda Beta

% Configuración de tipos de letra y tamaños de imprenta (Arial 7pt)
set(0, 'DefaultAxesFontName', 'Arial', ...
       'DefaultAxesFontSize', 7, ...
       'DefaultTextFontName', 'Arial', ...
       'DefaultTextFontSize', 7, ...
       'DefaultFigureColor', 'w');

bnd    = S.bands_hz;          % θ4-7 / α8-12 / β13-30 (desde config)
bnd_nm = {'\theta','\alpha','\beta'};


%% ============================================================================
%%SN1 ── Delta-espectro periódico FOOOF (S-nuevo 1)
% ============================================================================
fprintf('>>> [SN1]  FOOOF delta-espectro periódico (Ch − Nc)...\n');
if ~isfile(F_FOOOF)
    fprintf('  [SKIP] %s\n\n', F_FOOOF);
else
    W = load(F_FOOOF);  GR = W.(firstfield(W));
    if ~isfield(GR,'Cases') && isfield(GR,'casos'); GR.Cases = GR.casos; end
    f1   = GR.Cases.f(:);
    R_ch1 = real(GR.Cases.Res_Ch);   
    R_nc1 = real(GR.Cases.Res_Nc);
    nCas1 = size(R_ch1,2);
    f_mask1 = f1>=2 & f1<=35;
    f_p1    = f1(f_mask1);
    dR1     = R_ch1(f_mask1,:) - R_nc1(f_mask1,:);   
    mu1  = nanmean(dR1,2);
    se1  = nanstd(dR1,0,2)/sqrt(nCas1);
    
    pbin1 = ones(numel(f_p1),1);
    for fi=1:numel(f_p1)
        v = dR1(fi,:); v = v(~isnan(v));
        if numel(v)>=5
            pbin1(fi) = signrank(v);
        end
    end
    sig1 = pbin1 < 0.05;
    
    fig1 = figure('Name','FOOOF_DeltaPeriodic',...
                  'Units','centimeters','Position',[2 2 fig_width fig_height],'Visible','off');
    ax1  = axes(fig1); hold(ax1,'on');
    
    % Ampliamos el límite superior para dar espacio sin pisar 'beta'
    ylims1_tmp = [min(mu1-se1)-0.08,  max(mu1+se1)+0.28];
    for bi=1:3
        if bnd{bi}(1)<35 && bnd{bi}(2)>2
            patch(ax1,[bnd{bi}(1) bnd{bi}(2) bnd{bi}(2) bnd{bi}(1)],...
                  [ylims1_tmp(1) ylims1_tmp(1) ylims1_tmp(2) ylims1_tmp(2)],...
                  [0.97 0.97 0.97],'EdgeColor','none','HandleVisibility','off');
            % Bajamos el texto de la banda un poco para alejarlo del borde superior
            text(ax1,mean(bnd{bi}),ylims1_tmp(2)-0.10*range(ylims1_tmp),...
                 bnd_nm{bi},'HorizontalAlignment','center','FontSize',8,...
                 'Color',[.50 .50 .50],'Interpreter','tex');
        end
    end
    
    fill(ax1,[f_p1;flipud(f_p1)],[mu1+se1;flipud(mu1-se1)],...
         c_case_fill,'FaceAlpha',.50,'EdgeColor','none','HandleVisibility','off');
    
    % Leyenda limpia: Sin el 'n=30'
    plot(ax1,f_p1,mu1,'-','Color',c_case,'LineWidth',1.2,...
         'DisplayName','Cases Ch-Nc');
    yline(ax1,0,'k-','LineWidth',0.4,'Color',[.4 .4 .4],'HandleVisibility','off');
    
    % Línea continua para la ventana significativa en la base inferior
    if any(sig1)
        ymark1 = ylims1_tmp(1) + 0.02*range(ylims1_tmp); % Ajustado levemente hacia abajo
        sig_y = repmat(ymark1, size(f_p1));
        sig_y(~sig1) = NaN; % Mantiene la línea únicamente donde p < 0.05
        plot(ax1, f_p1, sig_y, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8, ...
             'DisplayName', 'Significant window');
    end
    
    set(ax1,'XLim',[2 35],'YLim',ylims1_tmp,...
        'XTick',[4 7 8 13 15 30],'XTickLabel',{'4','7','8','13','15','30'},...
        'TickDir','out','Box','off','XGrid','on','YGrid','on',...
        'GridAlpha',0.08,'Layer','top','LineWidth',0.5,'TickLength',[0.02 0.02]);
    xlabel(ax1,'Frequency (Hz)');
    
    % Eje Y limpio: Solo dB
    ylabel(ax1,'\Delta Periodic Power (dB)','Interpreter','tex');
    title(ax1,'FOOOF Residual: Chew vs No-Chew (Cases)','FontSize',8,'FontWeight','normal');
    
    % Leyenda reubicada abajo a la izquierda para máxima claridad visual
    legend(ax1,'show','Location','southwest','Box','off','FontSize',6);
    
    exportgraphics(fig1,fullfile(DIR_SUPP,'FigSN1_FOOOF_DeltaPeriodic.png'),'Resolution',DPI);
    close(fig1);
    fprintf('  -> FigSN1_FOOOF_DeltaPeriodic.png\n\n');
end
%% ============================================================================
%%SN2 ── Espectros FOOOF individuales + distribución tCF (S-nuevo 2)
% ============================================================================
fprintf('>>> [SN2]  FOOOF espectros individuales + histograma tCF...\n');
if ~isfile(F_FOOOF)
    fprintf('  [SKIP] %s\n\n', F_FOOOF);
else
    if ~exist('GR','var') || ~isfield(GR,'Cases')
        W = load(F_FOOOF);  GR = W.(firstfield(W));
        if ~isfield(GR,'Cases') && isfield(GR,'casos'); GR.Cases = GR.casos; end
    end
    f2     = GR.Cases.f(:);
    R_ch2  = real(GR.Cases.Res_Ch);
    R_nc2  = real(GR.Cases.Res_Nc);
    tCF_ch2 = double(GR.Cases.tCF_Ch(:));  tCF_ch2(tCF_ch2<1) = NaN;
    tCF_nc2 = double(GR.Cases.tCF_Nc(:));  tCF_nc2(tCF_nc2<1) = NaN;
    nS2     = size(R_ch2,2);
    f_mask2 = f2>=2 & f2<=16;
    f_p2    = f2(f_mask2);
    R_ch2_p = R_ch2(f_mask2,:);
    R_nc2_p = R_nc2(f_mask2,:);
    
    all_r2  = [R_ch2_p(:); R_nc2_p(:)];
    ylo2 = prctile(all_r2, 1) - 0.1;
    yhi2 = prctile(all_r2,99) + 0.3;
    ylim2 = [ylo2 yhi2];
    
    cond2 = {struct('R',R_nc2_p,'tCF',tCF_nc2,'clr',c_ctrl,'clr_f',c_ctrl_fill, ...
                    'lbl','No Chew','fname','FigSN2a_FOOOF_Spectra_NoChew.png'), ...
             struct('R',R_ch2_p, 'tCF',tCF_ch2, 'clr',c_case,'clr_f',c_case_fill, ...
                    'lbl','Chew',   'fname','FigSN2b_FOOOF_Spectra_Chew.png')};
                
    for ci = 1:2
        cd2  = cond2{ci};
        R_u  = cd2.R;  tCF_u = cd2.tCF;
        clr  = cd2.clr;  clr_f = cd2.clr_f;
        
        fig2 = figure('Name',['FOOOF_Spectra_' num2str(ci)],...
                      'Units','centimeters','Position',[2.5 2.5 fig_width fig_height],'Visible','off');
        ax2  = axes(fig2); hold(ax2,'on');
        
        % Trazos individuales (se ocultan de la leyenda para evitar desorden)
        for s = 1:nS2
            is_edge = ~isnan(tCF_u(s)) && (tCF_u(s)<4.5 || tCF_u(s)>6.5);
            ec = [0.85 0.50 0.10 0.40]*is_edge + [clr 0.15]*(~is_edge);
            plot(ax2, f_p2, R_u(:,s), '-', 'Color', ec, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        
        % Sombreado de error y promedio
        m2  = nanmean(R_u,2);  se2 = nanstd(R_u,0,2)/sqrt(nS2);
        fill(ax2,[f_p2;flipud(f_p2)],[m2+se2;flipud(m2-se2)],...
             clr_f,'FaceAlpha',0.50,'EdgeColor','none','HandleVisibility','off');
        
        % 1. Etiqueta para el promedio en la leyenda
        plot(ax2,f_p2,m2,'-','Color',clr,'LineWidth',1.3,'DisplayName','Mean \pm SE');
        
        % Líneas ficticias para explicar los colores en la leyenda ordenadamente debajo de Mean +- SE
        plot(ax2, NaN, NaN, '-', 'Color', clr, 'LineWidth', 1.0, 'DisplayName', 'Individual (central peak)');
        plot(ax2, NaN, NaN, '-', 'Color', [0.85 0.50 0.10], 'LineWidth', 1.0, 'DisplayName', 'Individual (edge peak)');
        
        % Sombreado de la banda Theta
        patch(ax2,[4 7 7 4],[ylim2(1) ylim2(1) ylim2(2) ylim2(2)],...
              [1.00 0.96 0.85],'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
        text(ax2,5.5,ylim2(2)-0.08*range(ylim2),'\theta','FontSize',9,...
             'HorizontalAlignment','center','Color',[0.60 0.40 0.10],'Interpreter','tex');
        
        % Marcas inferiores (ticks) de frecuencia central detectada
        has_pk2 = ~isnan(tCF_u);
        for s = 1:nS2
            if has_pk2(s)
                is_e = tCF_u(s)<4.5 || tCF_u(s)>6.5;
                mk = [0.85 0.50 0.10]*is_e + clr*(~is_e);
                plot(ax2, tCF_u(s), ylim2(1)+0.04*range(ylim2), ...
                     '|', 'Color', mk, 'MarkerSize', 5, 'LineWidth', 1.0, ...
                     'HandleVisibility', 'off');
            end
        end
        
        yline(ax2,0,'k-','LineWidth',0.4,'Color',[.4 .4 .4],'HandleVisibility','off');
        set(ax2,'XLim',[2 16],'YLim',ylim2,...
            'XTick',[4 7 8 13],'XTickLabel',{'4','7','8','13'},...
            'TickDir','out','Box','off','XGrid','on','YGrid','on',...
            'GridAlpha',0.08,'Layer','top','LineWidth',0.5,'TickLength',[0.02 0.02]);
        xlabel(ax2,'Frequency (Hz)');
        
        % Eje Y limpio
        ylabel(ax2,'Periodic Power (dB)');
        
        % Título limpio (sin 'n=30')
        title(ax2,sprintf('FOOOF Residuals: Cases %s', cd2.lbl),'FontSize',8,'FontWeight','normal');
        
        legend(ax2,'show','Location','northeast','Box','off','FontSize',6);
        
        exportgraphics(fig2,fullfile(DIR_SUPP,cd2.fname),'Resolution',DPI);
        close(fig2);
        fprintf('  -> %s\n', cd2.fname);
    end
    
    % ── 2C: Histograma tCF (Estilo escalera limpio)
    fig2c = figure('Name','FOOOF_tCF_Hist',...
                   'Units','centimeters','Position',[3 3 fig_width fig_height],'Visible','off');
    ax2c  = axes(fig2c); hold(ax2c,'on');
    edges2c = 3:0.4:9;
    v_nc2c  = tCF_nc2(~isnan(tCF_nc2));
    v_ch2c  = tCF_ch2(~isnan(tCF_ch2));
    
    % Histogramas con nombres limpios para la leyenda
    histogram(ax2c, v_nc2c, edges2c, 'FaceColor', c_ctrl, 'FaceAlpha', .40,...
              'EdgeColor', c_ctrl*0.70, 'LineWidth', 0.6,...
              'DisplayName', 'No Chew');
    histogram(ax2c, v_ch2c, edges2c, 'DisplayStyle', 'stairs',...
              'EdgeColor', c_case, 'LineWidth', 1.2,...
              'DisplayName', 'Chew');
              
    % Línea ficticia para agregar el significado de las líneas punteadas a la leyenda
    plot(ax2c, NaN, NaN, 'k:', 'LineWidth', 1.0, 'Color', [.5 .5 .5], 'DisplayName', '\theta boundaries');
    
    yl2c = ylim(ax2c);
    % Líneas reales punteadas ocultas de la leyenda para no duplicar
    xline(ax2c,4.0,'k:','LineWidth',1.0,'Color',[.5 .5 .5],'HandleVisibility','off');
    xline(ax2c,7.0,'k:','LineWidth',1.0,'Color',[.5 .5 .5],'HandleVisibility','off');
    
    % Triángulos de las medianas ocultos en la leyenda
    plot(ax2c,median(v_nc2c),yl2c(1)+0.03*range(yl2c),'v','MarkerSize',6,'Color',c_ctrl,...
         'MarkerFaceColor',c_ctrl,'HandleVisibility','off');
    plot(ax2c,median(v_ch2c),yl2c(1)+0.03*range(yl2c),'v','MarkerSize',6,'Color',c_case,...
         'MarkerFaceColor',c_case,'HandleVisibility','off');
    
    both2c = ~isnan(tCF_ch2) & ~isnan(tCF_nc2);
    if sum(both2c)>=5
        [p_cf2,~] = signrank(tCF_ch2(both2c), tCF_nc2(both2c));
        text(ax2c,0.95,0.92,sprintf('Shift p = %.4f',p_cf2),...
             'Units','normalized','HorizontalAlignment','right',...
             'FontSize',7,'Color',[.3 .3 .3]);
    end
    
    set(ax2c,'TickDir','out','Box','off','XGrid','on','YGrid','on',...
        'GridAlpha',0.08,'LineWidth',0.5,'TickLength',[0.02 0.02]);
    xlabel(ax2c,'Theta Peak Frequency (Hz)');
    
    % Eje Y claro que exprese frecuencia de ocurrencia 
    ylabel(ax2c,'Number of subjects');
    title(ax2c,'Theta Central Frequency Peak Distribution','FontSize',8,'FontWeight','normal');
    legend(ax2c,'show','Location','northwest','Box','off','FontSize',6);
    
    exportgraphics(fig2c,fullfile(DIR_SUPP,'FigSN2c_FOOOF_tCF_Distribution.png'),'Resolution',DPI);
    close(fig2c);
    fprintf('  -> FigSN2c_FOOOF_tCF_Distribution.png\n\n');
end


%% ============================================================================
%%SN3 ── Beta-PAC especificidad espectral + dinámica temporal (S-nuevo 3)
% ============================================================================
fprintf('>>> [SN3]  beta-PAC especificidad + dinámica temporal...\n');
if ~isfile(F_PAC4G)
    fprintf('  [SKIP] %s\n\n', F_PAC4G);
else
    WS4 = load(F_PAC4G);
    try
        z_th_E = WS4.B.Casos.Ch.Theta(:,1);
        z_al_E = WS4.B.Casos.Ch.Alpha(:,1);
        z_be_E = WS4.B.Casos.Ch.Beta(:,1);
        % Reorder columns to chronological: Early(col1), Mid/200-700ms(col3), Late(col2)
        % Original workspace: col1=Early(0-300ms), col2=Late(300-900ms), col3=Mid(200-700ms)
        z_be_Ch3 = [WS4.B.Casos.Ch.Beta(:,1), WS4.B.Casos.Ch.Beta(:,3), WS4.B.Casos.Ch.Beta(:,2)];
        z_be_Nc3 = [WS4.B.Casos.Nc.Beta(:,1), WS4.B.Casos.Nc.Beta(:,3), WS4.B.Casos.Nc.Beta(:,2)];
        sn3_ok = true;
    catch ME
        fprintf('  [ERROR] Estructura inesperada: %s\n\n', ME.message);
        sn3_ok = false;
    end
    
    if sn3_ok
        nSp3   = numel(z_th_E);
        bands3 = {z_th_E, z_al_E, z_be_E};
        pvals3 = [signrank(z_th_E) signrank(z_al_E) signrank(z_be_E)];
        win_lbl3     = {'Early','Mid','Late'};   % chronological: 0-300ms, 200-700ms, 300-900ms
        p_be_wins_Ch = arrayfun(@(w) signrank(z_be_Ch3(:,w)), 1:3);
        
        % ── SN3a: Especificidad espectral (Barras compactas y alineación superior)
        fig3a = figure('Name','BetaPAC_Specificity',...
                       'Units','centimeters','Position',[3.5 3.5 fig_width fig_height],'Visible','off');
        ax3a  = axes(fig3a); hold(ax3a,'on');
        clrs3  = {c_ctrl,[0.68 0.68 0.22],c_beta};
        bnd_lbl3 = {'\theta','\alpha','\beta'};
        
        for bi = 1:3
            d3 = bands3{bi};  d3 = d3(~isnan(d3));
            mu3a = mean(d3);  se3a = std(d3)/sqrt(numel(d3));
            bar(ax3a,bi,mu3a,0.50,'FaceColor',clrs3{bi},'EdgeColor','none');
            
            % Barras de error de un solo lado (solo hacia arriba) para coherencia visual con 3b
            errorbar(ax3a,bi,mu3a,zeros(size(mu3a)),se3a,'k','LineStyle','none','LineWidth',0.8,'CapSize',4);
            
            jit3 = (rand(numel(d3),1)-.5)*0.18;
            scatter(ax3a,bi+jit3,d3,8, 'k','filled','MarkerFaceAlpha',.20);
        end
        
        yline(ax3a,0,'k-','LineWidth',0.4,'Color',[.4 .4 .4]);
        
        set(ax3a,'XLim',[0.5 3.5],'XTick',1:3,'XTickLabel',bnd_lbl3,'TickDir','out','Box','off',...
            'XGrid','on','YGrid','on','GridAlpha',0.08,'FontSize',7,...
            'TickLabelInterpreter','tex','LineWidth',0.5,'TickLength',[0.02 0.02]);
            
        yl_3a = ylim(ax3a);
        y_star_3a = yl_3a(2) - 0.02 * range(yl_3a); 
        for bi = 1:3
            sym3 = p2stars(pvals3(bi));
            text(ax3a, bi, y_star_3a, sym3, 'HorizontalAlignment','center',...
                 'VerticalAlignment','top','FontSize',7,'Color','k');
        end
        
        ylabel(ax3a,'\DeltazMI (Chew vs 0)','Interpreter','tex');
        title(ax3a,'\beta-PAC Spectral Specificity','FontSize',8,'FontWeight','normal','Interpreter','tex');
        
        exportgraphics(fig3a,fullfile(DIR_SUPP,'FigSN3a_BetaPAC_SpectralSpecificity.png'),'Resolution',DPI);
        close(fig3a);
        fprintf('  -> FigSN3a_BetaPAC_SpectralSpecificity.png\n');
        
        % ── SN3b: Dinámica temporal β-PAC (Barras agrupadas y error unidireccional)
        fig3b = figure('Name','BetaPAC_TemporalDynamics',...
                       'Units','centimeters','Position',[4 4 fig_width fig_height],'Visible','off');
        ax3b  = axes(fig3b); hold(ax3b,'on');
        mu_ch3b = mean(z_be_Ch3,1,'omitnan');
        se_ch3b = std(z_be_Ch3,0,1,'omitnan')/sqrt(nSp3);
        mu_nc3b = mean(z_be_Nc3,1,'omitnan');
        se_nc3b = std(z_be_Nc3,0,1,'omitnan')/sqrt(nSp3);
        
        b3 = bar(ax3b, 1:3, [mu_nc3b; mu_ch3b]', 'grouped', 'EdgeColor', 'none');
        b3(1).FaceColor = c_ctrl; 
        b3(2).FaceColor = c_case; 
        
        x_nc3 = b3(1).XEndPoints;
        x_ch3 = b3(2).XEndPoints;
        
        % Corrección: zeros(size(...)) en el cuarto parámetro (yneg) fuerza el error solo hacia arriba
        errorbar(ax3b, x_nc3, mu_nc3b, zeros(size(mu_nc3b)), se_nc3b, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
        errorbar(ax3b, x_ch3, mu_ch3b, zeros(size(mu_ch3b)), se_ch3b, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
        
        yline(ax3b,0,'k-','LineWidth',0.4,'Color',[.4 .4 .4],'HandleVisibility','off');
        
        set(ax3b,'XLim',[0.5 3.5],'XTick',1:3,'XTickLabel',win_lbl3,'TickDir','out','Box','off',...
            'XGrid','on','YGrid','on','GridAlpha',0.08,'FontSize',7,...
            'LineWidth',0.5,'TickLength',[0.02 0.02]);
            
        yl_3b = ylim(ax3b);
        y_star_3b = yl_3b(2) - 0.02 * range(yl_3b);
        for w = 1:3
            sym_ch3b = p2stars(p_be_wins_Ch(w));
            if p_be_wins_Ch(w)<.05
                text(ax3b, w, y_star_3b, sym_ch3b,...
                     'HorizontalAlignment','center','VerticalAlignment','top',...
                     'FontSize',8,'Color','k');
            end
        end
        
        ylabel(ax3b,'\DeltazMI','Interpreter','tex');
        title(ax3b,'Temporal Dynamics of \beta-PAC','FontSize',8,'FontWeight','normal','Interpreter','tex');
        
        legend([b3(1), b3(2)], {'No Chew', 'Chew'}, 'Location','northwest','Box','off','FontSize',6);
        
        exportgraphics(fig3b,fullfile(DIR_SUPP,'FigSN3b_BetaPAC_TemporalDynamics.png'),'Resolution',DPI);
        close(fig3b);
        fprintf('  -> FigSN3b_BetaPAC_TemporalDynamics.png\n\n');
    end
end

%% ============================================================================
%%SN4 ── Alpha-ERD curso temporal + topograma (S-nuevo 4)
% ============================================================================
fprintf('>>> [SN4]  Alpha-ERD curso temporal + topograma...\n');
have_tfraw4 = isfile(F_TF_CCH) && isfile(F_TF_CNC);
have_mets4  = isfile(F_TF_METRICS);

% --- topoplotIndie.m está en la carpeta del pipeline (junto a S0_paths) ---
if isempty(which('topoplotIndie'))
    addpath(S.dir_pip);
    if ~isempty(which('topoplotIndie'))
        fprintf('  [OK] topoplotIndie añadido al path: %s\n', S.dir_pip);
    else
        fprintf('  [AVISO] topoplotIndie no encontrado en %s → SN4b (topo) se omitirá.\n', S.dir_pip);
    end
end
% --------------------------------------------------------

if ~have_tfraw4 && ~have_mets4
    fprintf('  [SKIP] Sin datos TF disponibles.\n\n');
elseif ~have_tfraw4
    % ── Fallback: TF_band_metrics
    fprintf('  [FALLBACK] TF_band_metrics — 3 ventanas temporales.\n');
    TF4 = load(F_TF_METRICS,'TF_metrics');
    % TF_metrics win order: 1=Early(0-300ms), 2=Late(300-900ms), 3=Mid(200-700ms)
    % Reorder to chronological [Early, Mid, Late] = columns [1, 3, 2]
    alpha_ch4 = squeeze(TF4.TF_metrics(:,1,2,[1,3,2]));
    alpha_nc4 = squeeze(TF4.TF_metrics(:,2,2,[1,3,2]));
    nSt4 = size(alpha_ch4,1);
    mu_ch4f = nanmean(alpha_ch4,1);  se_ch4f = nanstd(alpha_ch4,0,1)/sqrt(nSt4);
    mu_nc4f = nanmean(alpha_nc4,1);  se_nc4f = nanstd(alpha_nc4,0,1)/sqrt(nSt4);
    
    fig4a = figure('Name','Alpha_ERD_Fallback',...
                   'Units','centimeters','Position',[4.5 4.5 fig_width fig_height],'Visible','off');
    ax4a  = axes(fig4a); hold(ax4a,'on');
    
    fill(ax4a,[1:3,3:-1:1],[mu_nc4f-se_nc4f,fliplr(mu_nc4f+se_nc4f)],...
         c_ctrl_fill,'FaceAlpha',.40,'EdgeColor','none','HandleVisibility','off');
    fill(ax4a,[1:3,3:-1:1],[mu_ch4f-se_ch4f,fliplr(mu_ch4f+se_ch4f)],...
         c_case_fill,'FaceAlpha',.40,'EdgeColor','none','HandleVisibility','off');
    
    plot(ax4a,1:3,mu_nc4f,'o-','Color',c_ctrl,'LineWidth',1.0,...
         'MarkerFaceColor',c_ctrl,'MarkerSize',3.5,'DisplayName','No Chew');
    plot(ax4a,1:3,mu_ch4f,'o-','Color',c_case,'LineWidth',1.2,...
         'MarkerFaceColor',c_case,'MarkerSize',3.5,'DisplayName','Chew');
     
    set(ax4a,'XTick',1:3,'XTickLabel',{'Early','Mid','Late'},...
        'TickDir','out','Box','off','XGrid','on','YGrid','on','GridAlpha',0.08,...
        'LineWidth',0.5,'TickLength',[0.02 0.02]);
        
    yl_4af = ylim(ax4a);
    y_star_4af = yl_4af(2) - 0.02 * range(yl_4af);
    for w = 1:3
        [~, p_w] = signrank(alpha_ch4(:,w), alpha_nc4(:,w));
        if p_w < 0.05
            text(ax4a, w, y_star_4af, p2stars(p_w), 'HorizontalAlignment','center',...
                 'VerticalAlignment','top','FontSize',8,'Color','k');
        end
    end
    
    xlabel(ax4a,'Time Window');
    ylabel(ax4a,'\alpha Power (dB)','Interpreter','tex');
    title(ax4a,'\alpha-ERD by Time Window (F1, Fz, F2, FC1, FCz, FC2, AFz)','FontSize',8,'FontWeight','normal','Interpreter','tex');
    legend(ax4a,'show','Location','northeast','Box','off','FontSize',6);
    
    exportgraphics(fig4a,fullfile(DIR_SUPP,'FigSN4a_Alpha_ERD_TimeCourse.png'),'Resolution',DPI);
    close(fig4a);
    fprintf('  -> FigSN4a_Alpha_ERD_TimeCourse.png (fallback)\n\n');
else
    % ── Versión Completa: tfraw
    S_ch4 = load(F_TF_CCH);   S_nc4 = load(F_TF_CNC);
    tf_ch4 = get6D(S_ch4);    tf_nc4 = get6D(S_nc4);
    pwr_ch4 = squeeze(tf_ch4(2,:,:,:,:,1));  
    pwr_nc4 = squeeze(tf_nc4(2,:,:,:,:,1));
    [times_tf4, frex_tf4] = get_tf_meta(S_ch4);
    nSt4 = min(size(pwr_ch4,4), size(pwr_nc4,4));
    
    has_eeg4 = ~isempty(which('topoplot')) || ~isempty(which('topoplotIndie'));
    if ~has_eeg4 && exist(DIR_EEGLAB,'dir')
        addpath(DIR_EEGLAB); eeglab nogui;
        has_eeg4 = ~isempty(which('topoplot')) || ~isempty(which('topoplotIndie'));
    end
    
    chanlocs4 = [];  roi4 = [];
    electrodes_list = {'F1','Fz','F2','FC1','FCz','FC2','AFz'}; %Es esta la lista?
    if has_eeg4
        sf4 = dir(fullfile(DIR_EPOCH,'*_Ch.set'));
        if ~isempty(sf4)
            EEGt4 = pop_loadset('filename',sf4(1).name,'filepath',DIR_EPOCH,'loadmode','info');
            n64_4 = min(size(pwr_ch4,1), numel(EEGt4.chanlocs));
            chanlocs4 = EEGt4.chanlocs(1:n64_4);
            for lab4 = electrodes_list
                ix4 = find(strcmpi({chanlocs4.labels},lab4{1}),1);
                if ~isempty(ix4); roi4(end+1)=ix4; end %#ok
            end
        end
    end
    if isempty(roi4); roi4 = 1:min(7,size(pwr_ch4,1)); end
    
    alpha4 = frex_tf4>=bnd{2}(1) & frex_tf4<=bnd{2}(2);   % alpha desde config (8-12)
    if ~any(alpha4)
        [~,ia4]=min(abs(frex_tf4-bnd{2}(1))); [~,ib4]=min(abs(frex_tf4-bnd{2}(2)));
        alpha4=false(size(frex_tf4)); alpha4(ia4:ib4)=true;
    end
    
    tc_ch4 = squeeze(mean(mean(pwr_ch4(roi4,alpha4,:,1:nSt4),1),2));
    tc_nc4 = squeeze(mean(mean(pwr_nc4(roi4,alpha4,:,1:nSt4),1),2));
    if size(tc_ch4,1)~=numel(times_tf4); tc_ch4=tc_ch4'; tc_nc4=tc_nc4'; end
    mu_ch4 = nanmean(tc_ch4,2);  se_ch4 = nanstd(tc_ch4,0,2)/sqrt(nSt4);
    mu_nc4 = nanmean(tc_nc4,2);  se_nc4 = nanstd(tc_nc4,0,2)/sqrt(nSt4);
    
    % Evaluar significancia punto a punto (Time Course)
    p_tc4 = ones(size(times_tf4));
    for ti = 1:numel(times_tf4)
        [~, p_tc4(ti)] = signrank(tc_ch4(ti,:)', tc_nc4(ti,:)');
    end
    sig_tc4 = p_tc4 < 0.05;
    
    % ── SN4a: Curso temporal de alta resolución
    fig4a = figure('Name','Alpha_ERD_TimeCourse',...
                   'Units','centimeters','Position',[4.5 4.5 fig_width fig_height],'Visible','off');
    ax4a  = axes(fig4a); hold(ax4a,'on');
    
    fill(ax4a,[times_tf4(:);flipud(times_tf4(:))],[mu_nc4+se_nc4;flipud(mu_nc4-se_nc4)],...
         c_ctrl_fill,'FaceAlpha',.40,'EdgeColor','none','HandleVisibility','off');
    fill(ax4a,[times_tf4(:);flipud(times_tf4(:))],[mu_ch4+se_ch4;flipud(mu_ch4-se_ch4)],...
         c_case_fill,'FaceAlpha',.40,'EdgeColor','none','HandleVisibility','off');
    
    plot(ax4a,times_tf4,mu_nc4,'-','Color',c_ctrl,'LineWidth',1.0,'DisplayName','No Chew');
    plot(ax4a,times_tf4,mu_ch4,'-','Color',c_case,'LineWidth',1.2,'DisplayName','Chew');
    
    xline(ax4a,0,'k:','LineWidth',0.8,'HandleVisibility','off');
    yline(ax4a,0,'k-','LineWidth',0.4,'Color',[.4 .4 .4],'HandleVisibility','off');
    
    yla4 = ylim(ax4a);
    patch(ax4a,[200 700 700 200],[yla4(1) yla4(1) yla4(2) yla4(2)],...
          [.90 .90 .90],'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
     
    if any(sig_tc4)
        y_sig_line = yla4(2) - 0.05 * range(yla4); 
        sig_points_y = repmat(y_sig_line, size(times_tf4));
        sig_points_y(~sig_tc4) = NaN;
        plot(ax4a, times_tf4, sig_points_y, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8, 'DisplayName', 'Significant difference');
    end
     
    set(ax4a,'XLim',[-200 1200],'TickDir','out','Box','off',...
        'XGrid','on','YGrid','on','GridAlpha',0.08,'Layer','top',...
        'LineWidth',0.5,'TickLength',[0.02 0.02]);
    xlabel(ax4a,'Time (ms)');
    ylabel(ax4a,'\alpha Power (dB)','Interpreter','tex');
    title(ax4a,'\alpha-ERD Time Course (F1, Fz, F2, FC1, FCz, FC2, AFz)','FontSize',8,'FontWeight','normal','Interpreter','tex');
    legend(ax4a,'show','Location','southwest','Box','off','FontSize',6);
    
    exportgraphics(fig4a,fullfile(DIR_SUPP,'FigSN4a_Alpha_ERD_TimeCourse.png'),'Resolution',DPI);
    close(fig4a);
    fprintf('  -> FigSN4a_Alpha_ERD_TimeCourse.png\n');
    
    % ── SN4b: Topograma ERD Ventana Activa (Edición topoplotIndie + Fondo Blanco + Escala Simétrica)
    if has_eeg4 && ~isempty(chanlocs4)
        act4 = times_tf4>=200 & times_tf4<=700;
        n64_4 = numel(chanlocs4);
        erd_topo4 = zeros(n64_4,1);
        for s = 1:nSt4
            pc4 = squeeze(pwr_ch4(1:n64_4,alpha4,:,s));
            pn4 = squeeze(pwr_nc4(1:n64_4,alpha4,:,s));
            erd_topo4 = erd_topo4 + nanmean(pc4(:,act4),2) - nanmean(pn4(:,act4),2);
        end
        erd_topo4 = erd_topo4/nSt4;
        
        fig4b = figure('Name','Alpha_ERD_Topography',...
                       'Units','centimeters','Position',[5 5 fig_width fig_height],'Visible','off');
        set(fig4b, 'Color', 'w'); 
        ax4b  = axes(fig4b);
        set(ax4b, 'Color', 'w'); 
        
        try
            % 1. Encontrar el límite absoluto para forzar simetría
            clim4b = max(abs(erd_topo4))*1.05;
            
            % 2. Llamar a topoplotIndie usando SOLO los argumentos que soporta
            topoplotIndie(erd_topo4, chanlocs4, 'electrodes', 'on');
            
            % 3. Forzar el límite de color del eje de forma manual
            if exist('clim', 'file')
                clim(ax4b, [-clim4b clim4b]); % MATLAB R2022a o superior
            else
                caxis(ax4b, [-clim4b clim4b]); % Versiones anteriores a R2022a
            end
            
            cb4b = colorbar(ax4b, 'Location', 'eastoutside');
            set(cb4b, 'FontSize', 6, 'TickLength', 0.015, 'LineWidth', 0.4);
            ylabel(cb4b, '\alpha ERD  Ch-Nc (dB)', 'FontSize', 7, 'Interpreter', 'tex');
            title(ax4b, '\alpha-ERD Topography (200-700 ms)', 'FontSize', 8, 'FontWeight', 'normal', 'Interpreter', 'tex');
        catch ME4b
            text(.5,.5,sprintf('Error topoplotIndie:\n%s', ME4b.message),'Units','normalized','HorizontalAlignment','center','FontSize',6, 'Interpreter', 'none'); 
            axis off;
            fprintf('  [ERROR SN4b] Ocurrió un error al ejecutar topoplotIndie:\n  %s\n', ME4b.message);
        end
        exportgraphics(fig4b,fullfile(DIR_SUPP,'FigSN4b_Alpha_ERD_Topography.png'),'Resolution',DPI);
        close(fig4b);
        fprintf('  -> FigSN4b_Alpha_ERD_Topography.png\n');
    else
        fprintf('  [SKIP SN4b] EEGLAB/chanlocs no disponibles o no se encontró el archivo topoplotIndie.m.\n');
    end
    fprintf('\n');
end
%% ============================================================================
%%SN5 ── Diagrama de mediación (S-nuevo 5)
% ============================================================================
fprintf('>>> [SN5]  Mediation path diagram...\n');

% Variables estadísticas
% Valores del Reporte_Supplementary (01-Jun, bandas nuevas). ab = a·b (el +6.00
% impreso en el reporte es un bug; el indirecto real es a·b).
a_path  = -0.0856;
b_path  = +27.98;
ab_path = a_path*b_path;        % ≈ -2.40
ci_ab   = [-10.31, 4.57];
rho_bpac_ies = -0.4465;

% Lienzo CUADRADO estándar para collage
fig5 = figure('Name','Mediation_Path',...
              'Units','centimeters','Position',[5.5 5.5 fig_width fig_height],'Visible','off');
ax5  = axes(fig5,'Visible','off');
axis(ax5,[0 1 0 1]); hold(ax5,'on');

% Geometría de Triángulo Isósceles exacto
nodes5(1).xy = [0.20 0.38]; nodes5(1).lbl = {'\beta-PAC','(Early)'};
nodes5(2).xy = [0.50 0.88]; nodes5(2).lbl = {'\alpha-ERD','(Mid)'};
nodes5(3).xy = [0.80 0.38]; nodes5(3).lbl = {'IES'};

% Dimensiones relativas de las cajas
nw5 = 0.25;  nh5 = 0.16;
nclr5 = {c_case_fill; c_ctrl_fill; [0.90 0.85 0.75]};

% 1. DIBUJAR CAJAS
for n5 = 1:3
    rectangle('Position',[nodes5(n5).xy(1)-nw5/2, nodes5(n5).xy(2)-nh5/2, nw5, nh5],...
              'Curvature',0.15,'FaceColor',nclr5{n5},'EdgeColor',[.4 .4 .4],'LineWidth',1);
    text(ax5,nodes5(n5).xy(1),nodes5(n5).xy(2),nodes5(n5).lbl,...
         'HorizontalAlignment','center','FontSize',7,'Interpreter','tex');
end

% 2. DIBUJAR FLECHAS (Matemáticamente centradas y simétricas)
% La longitud fija (0.20) garantiza que sean idénticas y no toquen las cajas
arrow_len = 0.20; 
draw_centered_arrow(ax5, nodes5(1).xy, nodes5(2).xy, arrow_len); % a path
draw_centered_arrow(ax5, nodes5(2).xy, nodes5(3).xy, arrow_len); % b path
draw_centered_arrow(ax5, nodes5(1).xy, nodes5(3).xy, arrow_len); % rho path

% 3. ETIQUETAS DE RUTAS (Totalmente aisladas de las flechas)
% Path a (Arriba a la izquierda, por fuera del triángulo)
text(ax5, 0.20, 0.68, sprintf('a = %.3f\n(n.s.)',a_path),...
     'HorizontalAlignment','center','FontSize',6.5,'Color',[.4 .4 .4]);

% Path b (Arriba a la derecha, por fuera del triángulo)
text(ax5, 0.80, 0.68, sprintf('b = %+.2f ms/u\np < 0.05',b_path),...
     'HorizontalAlignment','center','FontSize',6.5,'Color',[.1 .1 .1]);

% Correlación directa (rho) (Centrado debajo de la flecha base)
text(ax5, 0.50, 0.28, sprintf('\\rho = %+.3f, p_{raw} = 0.014 (n.s. FDR)',rho_bpac_ies),...
     'HorizontalAlignment','center','FontSize',7,'Color',c_case*0.70,'Interpreter','tex');
 
% 4. ESTADÍSTICAS GLOBALES
text(ax5, 0.50, 0.12, sprintf('Indirect effect: ab = %+.2f ms  95%% CI [%.2f, %.2f] (n.s.)',ab_path,ci_ab(1),ci_ab(2)),...
     'HorizontalAlignment','center','FontSize',6.5,'Color',[.4 .4 .4]);

text(ax5, 0.50, 0.04, 'Steiger Z = +2.65,  p = 0.008  \rightarrow  |\rho_{\beta,IES}| > |\rho_{\alpha,IES}|',...
     'HorizontalAlignment','center','FontSize',7.5,'Color',[0.15 0.35 0.60],'Interpreter','tex');

exportgraphics(fig5,fullfile(DIR_SUPP,'FigSN5_Mediation_Path.png'),'Resolution',DPI);
close(fig5);
fprintf('  -> FigSN5_Mediation_Path.png\n\n');
%% ============================================================================
%%SN6 ── Pipeline diagram (S-nuevo 6)  → REPOSITORIO MANUAL
% ============================================================================
fprintf('>>> [SN6]  Pipeline preprocessing diagram\n');
fprintf('  *** REPOSITORIO MANUAL — Guardar flujograma externo como FigSN6_Pipeline.png ***\n');
fprintf('  Destino: %s\n\n', DIR_SUPP);


%% ============================================================================
%%SN7 ── FOOOF independence scatter (S-nuevo 7)
% ============================================================================
fprintf('>>> [SN7]  FOOOF independence: delta-exponente vs delta-theta...\n');
if ~isfile(F_FOOOF)
    fprintf('  [SKIP] %s\n\n', F_FOOOF);
else
    if ~exist('GR','var') || ~isfield(GR,'Cases')
        W = load(F_FOOOF);  GR = W.(firstfield(W));
        if ~isfield(GR,'Cases') && isfield(GR,'casos'); GR.Cases = GR.casos; end
    end
    f7   = GR.Cases.f(:);
    dexp = double(GR.Cases.exp_Ch(:)) - double(GR.Cases.exp_Nc(:));
    th_m7 = f7>=4 & f7<=7;
    dth7  = nanmean(real(GR.Cases.Res_Ch(th_m7,:)),1)' - nanmean(real(GR.Cases.Res_Nc(th_m7,:)),1)';
    
    v7 = ~isnan(dexp) & ~isnan(dth7);
    x7 = dexp(v7);  y7 = dth7(v7);  n7 = sum(v7);
    [rho7,p7] = corr(x7,y7,'Type','Spearman');
    
    % Forzar límites absolutos
    xlim_7 = [-1.2 0.2];
    ylim_7 = [-2.5 1.5];
    
    Pf7 = polyfit(x7,y7,1);
    
    % CORRECCIÓN: Calcular la línea de tendencia y CI a lo largo de todo el eje X
    xx7 = linspace(xlim_7(1), xlim_7(2), 100);
    yy7 = polyval(Pf7,xx7);
    
    nbt7=2000; rb7=nan(nbt7,1); pb7=nan(nbt7,numel(xx7));
    for k=1:nbt7
        ib7=randi(n7,n7,1);
        rb7(k)=corr(x7(ib7),y7(ib7),'Type','Spearman');
        pb7(k,:)=polyval(polyfit(x7(ib7),y7(ib7),1),xx7);
    end
    ylo7 = prctile(pb7,2.5);  yhi7 = prctile(pb7,97.5);
    
    fig7 = figure('Name','FOOOF_Independence',...
                  'Units','centimeters','Position',[6 6 fig_width fig_height],'Visible','off');
    ax7  = axes(fig7); hold(ax7,'on');
    
    % Sombreado del intervalo de confianza que cruza de lado a lado
    fill(ax7,[xx7 fliplr(xx7)],[ylo7 fliplr(yhi7)],[.5 .5 .5],...
         'FaceAlpha',0.12,'EdgeColor','none');
         
    % Nube de puntos de Casos unificada
    scatter(ax7,x7,y7,24,c_case,'filled','MarkerFaceAlpha',0.60,...
            'MarkerEdgeColor','w','LineWidth',0.4);
    
    % Línea de tendencia (punteada si p > 0.05)
    ls7 = '-'; if p7>.05; ls7='--'; end
    plot(ax7,xx7,yy7,ls7,'Color',[.30 .30 .30],'LineWidth',1.0);
    
    yline(ax7,0,'k-','LineWidth',0.4,'Color',[.5 .5 .5],'HandleVisibility','off');
    xline(ax7,0,'k-','LineWidth',0.4,'Color',[.5 .5 .5],'HandleVisibility','off');
    
    % Leyenda de correlación limpia
    ps7 = sprintf('p = %.3f',p7); 
    ann7 = sprintf('\\rho = %+.3f\n%s', rho7, ps7);
    
    text(ax7,.05,.18,ann7,'Units','normalized',...
         'HorizontalAlignment','left','VerticalAlignment','bottom',...
         'FontSize',6.5,'BackgroundColor','none', 'Interpreter', 'tex');
     
    % Texto de conclusión matemática
    if abs(rho7)<.25 || p7>.05
        istr7 = 'Aperiodic \perp Periodic \theta';  iclr7 = [0.15 0.45 0.20];
    else
        istr7 = 'Correlated Subcomponents';         iclr7 = [0.70 0.15 0.15];
    end
    text(ax7,.05,.06,istr7,'Units','normalized','FontSize',7,'FontWeight','normal','Color',iclr7, 'Interpreter', 'tex');
    
    % Setear explícitamente los límites de los ejes para enmarcar el sombreado
    set(ax7,'XLim',xlim_7,'YLim',ylim_7,'TickDir','out','Box','off','FontSize',7,...
        'XGrid','on','YGrid','on','GridAlpha',0.08,'Layer','top',...
        'LineWidth',0.5,'TickLength',[0.02 0.02]);
        
    % Ejes simplificados 
    xlabel(ax7,'\Delta Aperiodic Exponent','Interpreter','tex');
    ylabel(ax7,'\Delta \theta Periodic Power','Interpreter','tex');
    title(ax7,'Aperiodic vs Periodic Component Independence','FontSize',8,'FontWeight','normal');
    
    exportgraphics(fig7,fullfile(DIR_SUPP,'FigSN7_FOOOF_Independence.png'),'Resolution',DPI);
    close(fig7);
    fprintf('  -> FigSN7_FOOOF_Independence.png\n\n');
end
%% ── RESUMEN FINAL DE PROCESAMIENTO ──────────────────────────────────────────
fprintf('%s\n', repmat('=',1,72));
fprintf('  S_Supp_NewFigures.m — PARÁMETROS DE PUBLICACIÓN COMPLETADOS\n');
fprintf('  Ubicación de Salida: %s\n', DIR_SUPP);
fprintf('%s\n\n', repmat('=',1,72));


%% ============================================================================
%  FUNCIONES LOCALES CONTROLADAS
%% ============================================================================
function s = p2stars(p)
    if     p<.001; s='***';
    elseif p<.01;  s='**';
    elseif p<.05;  s='*';
    else;          s='n.s.';
    end
end

function name = firstfield(S)
    fn = fieldnames(S); name = fn{1};
end

function tf6d = get6D(S)
    fn = fieldnames(S);  tf6d = [];
    for k=1:numel(fn)
        if ndims(S.(fn{k}))==6; tf6d=S.(fn{k}); return; end
    end
    if isempty(tf6d); tf6d=S.(fn{1}); end
end

function [times,frex] = get_tf_meta(S)
    if     isfield(S,'eeg_times'); times=S.eeg_times;
    elseif isfield(S,'times');     times=S.times;
    else;  times=[]; end
    if     isfield(S,'frex');        frex=S.frex;
    elseif isfield(S,'f');           frex=S.f;
    elseif isfield(S,'frequencies'); frex=S.frequencies;
    else;  frex=[]; end
    if isempty(times)||isempty(frex)
        fn=fieldnames(S);
        for k=1:numel(fn)
            if ndims(S.(fn{k}))==6
                if isempty(frex);  frex=1:size(S.(fn{k}),3); end
                if isempty(times); times=linspace(-200,1300,size(S.(fn{k}),4)); end
                break;
            end
        end
    end
end

function draw_centered_arrow(ax, p1, p2, arrow_length)
    % Calcula la distancia total y las componentes del vector director
    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);
    d = sqrt(dx^2 + dy^2);
    
    % Vector unitario de la dirección
    ux = dx / d; 
    uy = dy / d;
    
    % Punto medio geométrico exacto entre las cajas
    mx = (p1(1) + p2(1)) / 2;
    my = (p1(2) + p2(2)) / 2;
    
    % Coordenadas de inicio y fin del tallo de la flecha
    sx = mx - (arrow_length / 2) * ux;
    sy = my - (arrow_length / 2) * uy;
    
    ex = mx + (arrow_length / 2) * ux;
    ey = my + (arrow_length / 2) * uy;
    
    % Dibujar la línea principal (tallo)
    plot(ax, [sx ex], [sy ey], 'Color', [.4 .4 .4], 'LineWidth', 1.0);
    
    % --- Construcción de la cabeza de la flecha (Triángulo perfecto) ---
    head_len = 0.025; % Largo de la punta
    head_wid = 0.022; % Ancho de la base de la punta
    
    % Base del triángulo
    bx = ex - head_len * ux;
    by = ey - head_len * uy;
    
    % Vector normal (perpendicular a la flecha)
    nx = -uy;
    ny = ux;
    
    % Esquinas de la base del triángulo
    c1x = bx + (head_wid / 2) * nx;
    c1y = by + (head_wid / 2) * ny;
    
    c2x = bx - (head_wid / 2) * nx;
    c2y = by - (head_wid / 2) * ny;
    
    % Dibujar y rellenar el triángulo
    fill(ax, [ex c1x c2x], [ey c1y c2y], [.4 .4 .4], 'EdgeColor', 'none');
end
