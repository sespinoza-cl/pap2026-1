%% ============================================================================
%  S6_TF_Behavior.m — Potencia TF (theta/alpha/beta) vs IES y ChewFreq
%  ----------------------------------------------------------------------------
%  Extrae potencia espectral por banda × ventana temporal desde las matrices
%  grupales TF (output S3_TF) y la correlaciona con:
%    (a) IES de la condición Chew y NoChew (conducta)
%    (b) Frecuencia masticatoria individual (chew_metrics.mat)
%
%  ROI EEG: F1 (ch12) y FC1 (ch16) — mismos canales que en S5_PAC_Compute.m
%  Métrica: Potencia Total baseline-corregida (dB), NPL pura.
%  ============================================================================
clear; clc; close all;

%% -- 1. RUTAS ----------------------------------------------------------------
% Fuente única de rutas — garantiza coherencia con el resto del pipeline
P         = S0_paths();
dir_tf    = P.dir_tf;            % C:\...\Desktop\Exp2\EEG\TF  (sólo lectura)
dir_out   = P.fig02b;            % D:\Exp2\Version_Mayo\Plots\Paper\Figure02b_TF_Correlaciones
dir_plots = dir_out;             % plots en la misma carpeta figura

f_beh     = P.file_beh;          % D:\Exp2\Version_Mayo\Pipelines\data_beh_tb_45.mat (N=45, auth.)
f_chew    = P.file_chew;         % C:\...\ChewFreq\chew_metrics.mat
f_rep     = fullfile(dir_out, 'Reporte_TF_Behavior.txt');

% S0_paths ya crea las carpetas; verificar por si se llama sin ella
if ~exist(dir_out, 'dir'), mkdir(dir_out); end

%% -- 2. PARÁMETROS -----------------------------------------------------------
% ROI: F1=12, FC1=16 
ROI_CH   = P.roi_ch;

% Bandas y ventanas — definición ÚNICA desde S0_paths (θ4-7 / α8-12 / β13-30)
BANDS    = P.bands;
BANDS_HZ = P.bands_hz;
nB       = numel(BANDS);

% Ventanas temporales (ms)
WINS     = P.wins;
WINS_MS  = P.wins_ms;
nW       = numel(WINS);
CONDS    = {'Ch','Nc'};

ITPC_WIN_MS  = [0 500];

%% -- 3. CARGAR MATRICES GRUPALES (Igual que en S2_TF_plot) ------------------
fprintf('>>> Cargando Matrices Grupales TF (Casos)...\n');
s_nc = load(fullfile(dir_tf, 'Group_Cases_Nc_tfraw.mat'));
s_ch = load(fullfile(dir_tf, 'Group_Cases_Ch_tfraw.mat'));

% Dim1=2 (NPL), Dim6=1 (Power) -> [64ch, frex, times, sujetos]
casos_nc_pwr = squeeze(s_nc.tfraw_pre_g(2,:,:,:,:,1));
casos_ch_pwr = squeeze(s_ch.tfraw_pre_g(2,:,:,:,:,1));

% Dim1=1 (Total Phase), Dim6=2 (ITPC) -> [64ch, frex, times, sujetos]
casos_ch_itpc = squeeze(s_ch.tfraw_pre_g(1,:,:,:,:,2));

times = s_ch.eeg_times;
frex  = s_ch.frex;
nS    = size(casos_ch_pwr, 4);

fprintf('  Datos cargados: %d Casos.\n', nS);

%% -- 4. EXTRACCIÓN DE MÉTRICAS POR SUJETO -----------------------------------
% TF_metrics:   [nS × nConds × nBands × nWins]  — Potencia NPL (dB)
TF_metrics   = nan(nS, 2, nB, nW);
% ITPC_metrics: [nS × 1 (solo Ch) × nBands]      — ITPC Total
ITPC_metrics = nan(nS, 1, nB);

fprintf('>>> Extrayendo promedios por ROI, Banda y Ventana...\n');
for ib = 1:nB
    fidx = dsearchn(frex', BANDS_HZ{ib}');
    
    % Extraer ITPC (Ventana 0-500ms fija)
    tidx_itpc = dsearchn(times', ITPC_WIN_MS');
    for i = 1:nS
        roi_itpc = casos_ch_itpc(ROI_CH, fidx(1):fidx(2), tidx_itpc(1):tidx_itpc(2), i);
        ITPC_metrics(i, 1, ib) = mean(roi_itpc, 'all', 'omitnan');
    end
    
    % Extraer Power NPL (Por ventana)
    for iw = 1:nW
        tidx = dsearchn(times', WINS_MS{iw}');
        for i = 1:nS
            roi_ch = casos_ch_pwr(ROI_CH, fidx(1):fidx(2), tidx(1):tidx(2), i);
            roi_nc = casos_nc_pwr(ROI_CH, fidx(1):fidx(2), tidx(1):tidx(2), i);
            
            TF_metrics(i, 1, ib, iw) = mean(roi_ch, 'all', 'omitnan'); % Cond 1: Ch
            TF_metrics(i, 2, ib, iw) = mean(roi_nc, 'all', 'omitnan'); % Cond 2: Nc
        end
    end
end

%% -- 5. CARGAR CONDUCTA Y CHEWFREQ ------------------------------------------
fprintf('>>> Cargando datos conductuales y masticatorios...\n');
tmp_beh = load(f_beh);

% Extracción directa como en S1_Conducta.mat
IES_Nc = double(tmp_beh.tb_data_45.casos.nochew.ies(:));
IES_Ch = double(tmp_beh.tb_data_45.casos.chew.ies(:));
Delta_IES = IES_Ch - IES_Nc;
casos_ids = string(tmp_beh.tb_data_45.casos.chew.Participantes);

% ChewFreq
tmp_chw = load(f_chew, 'T_freq');
T_freq  = tmp_chw.T_freq;
ChewFreq = nan(nS, 1);

for i = 1:nS
    idx = find(strcmp(T_freq.Sujeto, casos_ids(i)));
    if ~isempty(idx)
        ChewFreq(i) = mean([T_freq.Freq_Left(idx), T_freq.Freq_Right(idx)], 'omitnan');
    end
end

%% -- 5b. RATIO THETA/ALPHA ---------------------------------------------------
% ratio = theta_NPL_dB - alpha_NPL_dB  (mayor ratio = más theta relativo a alpha)
ratio_metrics = nan(nS, 2, nW);
for iw = 1:nW
    ratio_metrics(:, :, iw) = squeeze(TF_metrics(:, :, 1, iw)) - squeeze(TF_metrics(:, :, 2, iw)); 
end

%% -- 6. GUARDAR MÉTRICAS ----------------------------------------------------
T_TF = table(casos_ids(:), 'VariableNames', {'Sujeto'});
for ib = 1:nB
    for iw = 1:nW
        for ic = 1:2
            col_name = sprintf('TF_%s_%s_%s', BANDS{ib}, WINS{iw}, CONDS{ic});
            T_TF.(col_name) = TF_metrics(:, ic, ib, iw);
        end
    end
end
for iw = 1:nW
    T_TF.(sprintf('Ratio_ThAl_%s_Ch', WINS{iw})) = ratio_metrics(:, 1, iw);
    T_TF.(sprintf('Ratio_ThAl_%s_Nc', WINS{iw})) = ratio_metrics(:, 2, iw);
end
T_TF.IES_Ch   = IES_Ch;
T_TF.IES_Nc   = IES_Nc;
T_TF.Delta_IES = Delta_IES;
T_TF.ChewFreq  = ChewFreq;

save(fullfile(dir_out, 'TF_band_metrics.mat'), 'T_TF', 'TF_metrics', 'ITPC_metrics', ...
    'ratio_metrics', 'casos_ids', 'BANDS', 'BANDS_HZ', 'WINS', 'WINS_MS', 'ITPC_WIN_MS', 'ROI_CH', '-v7.3');
fprintf('>>> TF_band_metrics.mat guardado con éxito.\n');

%% -- 7. CORRELACIONES --------------------------------------------------------
% Matrices de resultados [nB × nW] para cada tipo de correlación
corr_types = {
    'TF_Ch vs IES_Ch',  'Ch', IES_Ch;
    'TF_Nc vs IES_Nc',  'Nc', IES_Nc;
    'Delta TF vs Delta IES', 'delta', Delta_IES;
    'TF_Ch vs ChewFreq', 'Ch', ChewFreq
};
nCorr = size(corr_types, 1);
rho_all = nan(nB, nW, nCorr);
p_all   = nan(nB, nW, nCorr);

for icorr = 1:nCorr
    ctype = corr_types{icorr,2}; beh_y = corr_types{icorr,3};
    for ib = 1:nB
        for iw = 1:nW
            if strcmp(ctype, 'Ch')
                x = TF_metrics(:, 1, ib, iw);
            elseif strcmp(ctype, 'Nc')
                x = TF_metrics(:, 2, ib, iw);
            else  % delta
                x = TF_metrics(:, 1, ib, iw) - TF_metrics(:, 2, ib, iw);
            end
            v = ~isnan(x) & ~isnan(beh_y);
            if sum(v) >= 5
                [rho_all(ib,iw,icorr), p_all(ib,iw,icorr)] = corr(x(v), beh_y(v), 'Type','Spearman');
            end
        end
    end
end

%% -- 7b. CORRELACIONES RATIO THETA/ALPHA ------------------------------------
rho_ratio = nan(nW, 3); p_ratio   = nan(nW, 3);
ratio_corr_types = {
    ratio_metrics(:,1,:), IES_Ch,    'Ratio_Ch vs IES_Ch';
    ratio_metrics(:,1,:) - ratio_metrics(:,2,:), Delta_IES, 'Delta_Ratio vs Delta_IES';
    ratio_metrics(:,1,:), ChewFreq,  'Ratio_Ch vs ChewFreq'
};
for ic2 = 1:3
    x_all = squeeze(ratio_corr_types{ic2,1}); beh_y = ratio_corr_types{ic2,2};
    for iw = 1:nW
        x = x_all(:,iw); v = ~isnan(x) & ~isnan(beh_y);
        if sum(v) >= 5
            [rho_ratio(iw,ic2), p_ratio(iw,ic2)] = corr(x(v), beh_y(v), 'Type','Spearman');
        end
    end
end

%% -- 7c. CORRELACIONES ITPC --------------------------------------------------
itpc_corr_types = {'ITPC_Ch vs IES_Ch', IES_Ch; 'ITPC_Ch vs ChewFreq', ChewFreq};
nICorr = size(itpc_corr_types, 1);
rho_itpc = nan(nB, nICorr); p_itpc   = nan(nB, nICorr);

for ic2 = 1:nICorr
    beh_y = itpc_corr_types{ic2,2};
    for ib = 1:nB
        x = ITPC_metrics(:, 1, ib); 
        v = ~isnan(x) & ~isnan(beh_y);
        if sum(v) >= 5
            [rho_itpc(ib,ic2), p_itpc(ib,ic2)] = corr(x(v), beh_y(v), 'Type','Spearman');
        end
    end
end

%% -- 8. REPORTE TXT (Mantenido intacto por tu excelente formato) -------------
fid = fopen(f_rep, 'w');
W   = @(s) fprintf(fid, '%s\n', s);
W(repmat('=',1,80));
W('  REPORTE: POTENCIA TF vs CONDUCTA (IES) y FRECUENCIA MASTICATORIA');
W(sprintf('  %s', datestr(now)));
W(repmat('=',1,80));
W(sprintf('  Casos N=%d  |  ROI: F1 (ch%d), FC1 (ch%d)', nS, ROI_CH(1), ROI_CH(2)));
W(sprintf('  Potencia NPL/Inducida (dB, baseline-corregida) — Wavelet Morlet'));
W('');
for icorr = 1:nCorr
    W(sprintf('  %s', corr_types{icorr,1}));
    W(repmat('-',1,80));
    hdr = sprintf('  %-12s', '');
    for iw = 1:nW; hdr = [hdr sprintf('  %-22s', WINS{iw})]; end
    W(hdr);
    for ib = 1:nB
        row = sprintf('  %-12s', BANDS{ib});
        for iw = 1:nW
            rho = rho_all(ib,iw,icorr); p = p_all(ib,iw,icorr); sym = sig_sym_local(p);
            row = [row sprintf('  rho=%+.2f p=%.3f%-4s ', rho, p, sym)];
        end
        W(row);
    end
    W('');
end
% --- Sección ITPC ---
W(''); W('  ITPC (Total, ventana 0–500 ms) — Correlaciones Spearman'); W(repmat('-',1,80));
hdr_itpc = sprintf('  %-14s', '');
for ic2 = 1:nICorr; hdr_itpc = [hdr_itpc sprintf('  %-26s', itpc_corr_types{ic2,1})]; end; W(hdr_itpc);
for ib = 1:nB
    row = sprintf('  %-14s', BANDS{ib});
    for ic2 = 1:nICorr
        rho = rho_itpc(ib,ic2); p = p_itpc(ib,ic2); sym = sig_sym_local(p);
        row = [row sprintf('  rho=%+.2f p=%.3f%-4s ', rho, p, sym)];
    end
    W(row);
end
% --- Sección Ratio theta/alpha ---
W(''); W('  RATIO THETA/ALPHA (theta_dB - alpha_dB) — Correlaciones Spearman');
W('  Interpretación: ratio bajo durante Ch = theta suprimido relativo a alpha = WM activo'); W(repmat('-',1,80));
ratio_corr_labels = {'Ratio_Ch vs IES_Ch', 'Delta_Ratio vs Delta_IES', 'Ratio_Ch vs ChewFreq'};
hdr_ratio = sprintf('  %-12s', '');
for iw = 1:nW; hdr_ratio = [hdr_ratio sprintf('  %-22s', WINS{iw})]; end; W(hdr_ratio);
for ic2 = 1:3
    row = sprintf('  %-12s', ratio_corr_labels{ic2});
    for iw = 1:nW
        rho = rho_ratio(iw,ic2); p = p_ratio(iw,ic2); sym = sig_sym_local(p);
        row = [row sprintf('  rho=%+.2f p=%.3f%-4s ', rho, p, sym)];
    end
    W(row);
end
fclose(fid);
fprintf('>>> Reporte guardado: %s\n', f_rep);

%% -- 9. FIGURAS --------------------------------------------------------------
c_ch  = [0.30 0.65 0.45]; c_nc  = [0.55 0.55 0.55]; c_chf = [0.40 0.65 0.80]; DPI   = 300;

% -- Fig 1: Grid TF_Ch vs IES_Ch 
f = figure('Color','w','Position',[60 60 1100 800],'Visible','off');
tl = tiledlayout(f, nB, nW, 'TileSpacing','compact','Padding','compact');
for ib = 1:nB
    for iw = 1:nW
        nexttile; hold on; x = TF_metrics(:,1,ib,iw); v = ~isnan(x) & ~isnan(IES_Ch);
        if sum(v) >= 5
            scatter(x(v), IES_Ch(v), 45, c_ch,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',c_ch*0.6);
            lf = polyfit(x(v), IES_Ch(v), 1); xl = linspace(min(x(v)), max(x(v)), 40);
            plot(xl, polyval(lf,xl), '--','Color',[c_ch 0.7],'LineWidth',1.8);
        end
        title(sprintf('%s-%s\n\\rho=%+.2f  p=%.3g%s', BANDS{ib}, WINS{iw}, rho_all(ib,iw,1), p_all(ib,iw,1), sig_sym_local(p_all(ib,iw,1))), 'FontSize',10,'FontWeight','bold');
        if iw==1; ylabel('IES (Ch)','FontSize',9); end
        if ib==nB; xlabel('Potencia TF (dB)','FontSize',9); end
        grid on; box off;
    end
end
title(tl,'Potencia TF Ch vs IES Ch — Spearman (Casos, ROI F1/FC1)','FontSize',13,'FontWeight','bold');
exportgraphics(f, fullfile(dir_plots,'Fig_TF_Ch_vs_IES_Ch.png'),'Resolution',DPI); close(f);

% -- Fig 2: Grid Delta TF vs Delta IES
f = figure('Color','w','Position',[60 60 1100 800],'Visible','off');
tl = tiledlayout(f, nB, nW, 'TileSpacing','compact','Padding','compact');
for ib = 1:nB
    for iw = 1:nW
        nexttile; hold on; dTF = TF_metrics(:,1,ib,iw) - TF_metrics(:,2,ib,iw); v = ~isnan(dTF) & ~isnan(Delta_IES);
        if sum(v) >= 5
            scatter(dTF(v), Delta_IES(v), 45, [0.4 0.4 0.7],'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',[0.3 0.3 0.6]);
            lf = polyfit(dTF(v), Delta_IES(v), 1); xl = linspace(min(dTF(v)), max(dTF(v)), 40);
            plot(xl, polyval(lf,xl), '--','Color',[0.4 0.4 0.7 0.7],'LineWidth',1.8);
        end
        yline(0,'k:','LineWidth',0.8); xline(0,'k:','LineWidth',0.8);
        title(sprintf('%s-%s\n\\rho=%+.2f  p=%.3g%s', BANDS{ib}, WINS{iw}, rho_all(ib,iw,3), p_all(ib,iw,3), sig_sym_local(p_all(ib,iw,3))), 'FontSize',10,'FontWeight','bold');
        if iw==1; ylabel('\Delta IES (Ch-Nc)','FontSize',9); end
        if ib==nB; xlabel('\Delta TF (Ch-Nc) dB','FontSize',9); end
        grid on; box off;
    end
end
title(tl,'\Delta Potencia TF (Ch-Nc) vs \Delta IES — Spearman (Casos)','FontSize',13,'FontWeight','bold');
exportgraphics(f, fullfile(dir_plots,'Fig_TF_Delta_vs_Delta_IES.png'),'Resolution',DPI); close(f);

% -- Fig 3: TF_Ch vs ChewFreq
f = figure('Color','w','Position',[60 60 1100 800],'Visible','off');
tl = tiledlayout(f, nB, nW, 'TileSpacing','compact','Padding','compact');
for ib = 1:nB
    for iw = 1:nW
        nexttile; hold on; x = TF_metrics(:,1,ib,iw); v = ~isnan(x) & ~isnan(ChewFreq);
        if sum(v) >= 5
            scatter(ChewFreq(v), x(v), 45, c_chf,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',c_chf*0.6);
            lf = polyfit(ChewFreq(v), x(v), 1); xl = linspace(min(ChewFreq(v))-0.05, max(ChewFreq(v))+0.05, 40);
            plot(xl, polyval(lf,xl), '--','Color',[c_chf 0.7],'LineWidth',1.8);
        end
        title(sprintf('%s-%s\n\\rho=%+.2f  p=%.3g%s', BANDS{ib}, WINS{iw}, rho_all(ib,iw,4), p_all(ib,iw,4), sig_sym_local(p_all(ib,iw,4))), 'FontSize',10,'FontWeight','bold');
        if iw==1; ylabel('Potencia TF (dB)','FontSize',9); end
        if ib==nB; xlabel('Freq masticatoria (Hz)','FontSize',9); end
        grid on; box off;
    end
end
title(tl,'Potencia TF Ch vs Frecuencia Masticatoria — Spearman (Casos)','FontSize',13,'FontWeight','bold');
exportgraphics(f, fullfile(dir_plots,'Fig_TF_vs_ChewFreq.png'),'Resolution',DPI); close(f);

% -- Fig 4: Heatmaps independientes rho (Estética Nature)

% Redefinimos los nombres de las bandas para usar símbolos griegos (TeX)
% Asegúrate de que nB y nW coincidan con la longitud de estos arreglos
BANDS_tex = {'\theta', '\alpha', '\beta'}; 
% Si tus ventanas son distintas, ajusta este arreglo:
WINS_tex = {'Early', 'Active', 'Late'}; 

% Matriz de pares: {rho, p-value, Título en Inglés, Nombre del Archivo}
pairs_hmap = {
    rho_all(:,:,1), p_all(:,:,1), 'TF Chew vs. IES Chew', 'Fig_TF_Behavior_RhoMap_IES.png';
    rho_all(:,:,4), p_all(:,:,4), 'TF Chew vs. Masticatory Freq.', 'Fig_TF_Behavior_RhoMap_ChewFreq.png'
};

for k = 1:size(pairs_hmap, 1)
    % 1. Crear figura independiente y ajustada
    f = figure('Color','w','Position',[100 100 450 400],'Visible','off');
    ax = axes('Parent', f);
    
    rho_m = pairs_hmap{k,1}; 
    p_m = pairs_hmap{k,2};
    
    % 2. Dibujar Heatmap
    imagesc(ax, rho_m'); 
    colormap(ax, redblue_local()); 
    clim(ax, [-0.7 0.7]);
    
    % 3. Estética Nature en los Ejes
    set(ax, 'TickDir', 'out', 'Box', 'off', 'FontSize', 12, 'FontName', 'Arial', ...
        'LineWidth', 1.2, 'XTick', 1:nB, 'XTickLabel', BANDS_tex, ...
        'YTick', 1:nW, 'YTickLabel', WINS_tex, 'YDir', 'normal', ...
        'TickLabelInterpreter', 'tex'); % Clave para renderizar \theta, \alpha, etc.
        
    % 4. Configuración del Colorbar
    cb = colorbar(ax); 
    cb.Label.String = 'Spearman \rho';
    cb.Label.FontSize = 13;
    cb.Label.FontName = 'Arial';
    cb.LineWidth = 1.2;
    cb.TickDirection = 'out';
    cb.Box = 'off';
    
    % 5. Textos de correlación y significancia
    for ib = 1:nB
        for iw = 1:nW
            sym = sig_sym_local(p_m(ib,iw)); 
            col = 'k'; 
            % Contraste adaptativo del texto según la intensidad del color
            if abs(rho_m(ib,iw)) > 0.4, col = 'w'; end
            
            text(ax, ib, iw, sprintf('%.2f%s', rho_m(ib,iw), sym), ...
                'HorizontalAlignment','center','FontSize',11, ...
                'FontWeight','bold','Color',col,'FontName','Arial');
        end
    end
    
    % 6. Título en inglés
    title(ax, pairs_hmap{k,3}, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % 7. Exportar cada gráfico independientemente
    exportgraphics(f, fullfile(dir_plots, pairs_hmap{k,4}), 'Resolution', DPI); 
    close(f);
end

% -- Fig 5: ITPC_Ch vs IES_Ch 
c_itpc = [0.80 0.45 0.20];
f = figure('Color','w','Position',[60 60 900 320],'Visible','off');
tl = tiledlayout(f, 1, nB, 'TileSpacing','compact','Padding','compact');
for ib = 1:nB
    nexttile; hold on; x = ITPC_metrics(:,1,ib); v = ~isnan(x) & ~isnan(IES_Ch);
    if sum(v) >= 5
        scatter(x(v), IES_Ch(v), 45, c_itpc,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',c_itpc*0.7);
        lf = polyfit(x(v), IES_Ch(v), 1); xl = linspace(min(x(v)), max(x(v)), 40);
        plot(xl, polyval(lf,xl), '--','Color',[c_itpc 0.7],'LineWidth',1.8);
    end
    title(sprintf('%s ITPC (0-500ms)\n\\rho=%+.2f  p=%.3g%s', BANDS{ib}, rho_itpc(ib,1), p_itpc(ib,1), sig_sym_local(p_itpc(ib,1))), 'FontSize',10,'FontWeight','bold');
    if ib==1; ylabel('IES (Ch)','FontSize',9); end
    xlabel('ITPC','FontSize',9); grid on; box off;
end
title(tl,'ITPC Ch (0-500ms) vs IES Ch — Spearman (Casos, ROI F1/FC1)','FontSize',13,'FontWeight','bold');
exportgraphics(f, fullfile(dir_plots,'Fig_ITPC_vs_IES_Ch.png'),'Resolution',DPI); close(f);

% -- Fig 6: ITPC_Ch vs ChewFreq 
f = figure('Color','w','Position',[60 60 900 320],'Visible','off');
tl = tiledlayout(f, 1, nB, 'TileSpacing','compact','Padding','compact');
for ib = 1:nB
    nexttile; hold on; x = ITPC_metrics(:,1,ib); v = ~isnan(x) & ~isnan(ChewFreq);
    if sum(v) >= 5
        scatter(ChewFreq(v), x(v), 45, c_chf,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',c_chf*0.6);
        lf = polyfit(ChewFreq(v), x(v), 1); xl = linspace(min(ChewFreq(v))-0.05, max(ChewFreq(v))+0.05, 40);
        plot(xl, polyval(lf,xl), '--','Color',[c_chf 0.7],'LineWidth',1.8);
    end
    title(sprintf('%s ITPC (0-500ms)\n\\rho=%+.2f  p=%.3g%s', BANDS{ib}, rho_itpc(ib,2), p_itpc(ib,2), sig_sym_local(p_itpc(ib,2))), 'FontSize',10,'FontWeight','bold');
    if ib==1; ylabel('ITPC','FontSize',9); end
    xlabel('Freq masticatoria (Hz)','FontSize',9); grid on; box off;
end
title(tl,'ITPC Ch (0-500ms) vs Frecuencia Masticatoria — Spearman (Casos)','FontSize',13,'FontWeight','bold');
exportgraphics(f, fullfile(dir_plots,'Fig_ITPC_vs_ChewFreq.png'),'Resolution',DPI); close(f);

% -- Fig 7: Ratio theta/alpha Ch vs IES_Ch
c_ratio = [0.25 0.60 0.35]; 
f = figure('Color','w','Position',[60 60 900 320],'Visible','off');
tl = tiledlayout(f, 1, nW, 'TileSpacing','compact','Padding','compact');
for iw = 1:nW
    nexttile; hold on; x = ratio_metrics(:,1,iw); v = ~isnan(x) & ~isnan(IES_Ch);
    if sum(v) >= 5
        scatter(x(v), IES_Ch(v), 45, c_ratio,'filled', 'MarkerFaceAlpha',0.75,'MarkerEdgeColor',c_ratio*0.7);
        lf = polyfit(x(v), IES_Ch(v), 1); xl = linspace(min(x(v)), max(x(v)), 40);
        plot(xl, polyval(lf,xl),'--','Color',[c_ratio 0.7],'LineWidth',1.8);
    end
    title(sprintf('Ratio %s\n\\rho=%+.2f  p=%.3g%s', WINS{iw}, rho_ratio(iw,1), p_ratio(iw,1), sig_sym_local(p_ratio(iw,1))), 'FontSize',10,'FontWeight','bold');
    if iw==1; ylabel('IES Ch (ms)','FontSize',9); end
    xlabel('\theta/\alpha ratio (dB)','FontSize',9); grid on; box off;
end
title(tl,'Ratio \theta/\alpha (NPL) Ch vs IES Ch — Spearman (Casos, ROI F1/FC1)','FontSize',13,'FontWeight','bold');
exportgraphics(f, fullfile(dir_plots,'Fig_Ratio_ThAl_vs_IES_Ch.png'),'Resolution',DPI); close(f);

fprintf('\n>>> Figuras guardadas en %s\n', dir_plots);
fprintf('>>> S6_TF_Behavior DONE.\n');

%% ============================================================================
%  LOCAL FUNCTIONS
%% ============================================================================
function s = sig_sym_local(p)
    if isnan(p) || p >= 0.05, s = '';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    else,             s = '*';
    end
end
function cmap = redblue_local()
    n = 64; r = [linspace(0.2,1,n)', linspace(0.2,1,n)', ones(n,1)];
    b = [ones(n,1), linspace(1,0.2,n)', linspace(1,0.2,n)']; cmap = [r; b];
end
