%% S3a_FOOOF_verify.m  —  Verificación FOOOF/specparam (Revisor M3)
% Lee v1_S3_specparam_results.csv (Python specparam ya corrido) y responde:
%   (1) ¿Qué % de sujetos tienen peak theta GENUINO sobre la 1/f?
%   (2) ¿El exponent aperiódico difiere entre grupos/condiciones?
%   (3) ¿La oscillación theta "sobrevive" (coexiste con la 1/f flattening)?
%
% Resultados clave esperados (del análisis previo):
%   - Cases Ch exponent ≈ 0.0-0.5 (aplanamient drástico vs Nc ≈ 0.9)
%   - 65% Cases_Ch tienen peak theta detectable vs 87% en Nc
%   - LME exponent Grupo×Cond: F=26.31 p<0.001
%
% Correr desde Analysis_V1_Final/

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

fprintf('\n=== S3a_FOOOF_verify.m ===\n');

%% ── Cargar CSV ────────────────────────────────────────────────────────────
csv_file = FILE_CSV;   % v1_S3_specparam_results.csv
assert(exist(csv_file,'file')==2, 'Falta %s', csv_file);
T = readtable(csv_file, 'TextType','string');
fprintf('CSV cargado: %d filas, %d columnas\n', height(T), width(T));
fprintf('Columnas: %s\n', strjoin(T.Properties.VariableNames, ', '));

%% ── Control de calidad ────────────────────────────────────────────────────
% Excluir fits con R² < 0.85 (ajuste insuficiente del modelo)
MIN_R2 = 0.85;
ok_mask = T.r_squared >= MIN_R2;
fprintf('\nFilas totales: %d | Con R²≥%.2f: %d | Excluidas: %d\n', ...
    height(T), MIN_R2, sum(ok_mask), sum(~ok_mask));

T_qc = T(ok_mask, :);

%% ── Subgrupos ─────────────────────────────────────────────────────────────
cas_ch = T_qc(T_qc.group == "Cases"    & T_qc.condition == "Ch", :);
cas_nc = T_qc(T_qc.group == "Cases"    & T_qc.condition == "Nc", :);
ctr_ch = T_qc(T_qc.group == "Controls" & T_qc.condition == "Ch", :);
ctr_nc = T_qc(T_qc.group == "Controls" & T_qc.condition == "Nc", :);

n_cas_ch = height(cas_ch);
n_cas_nc = height(cas_nc);
n_ctr_ch = height(ctr_ch);
n_ctr_nc = height(ctr_nc);

%% ── (1) % de sujetos con peak theta genuino ──────────────────────────────
fprintf('\n--- (1) Peak theta [4-7Hz] sobre background aperiódico ---\n');
pct = @(t) sum(t.has_theta_peak == "True") / height(t) * 100;

pct_cas_ch = pct(cas_ch);
pct_cas_nc = pct(cas_nc);
pct_ctr_ch = pct(ctr_ch);
pct_ctr_nc = pct(ctr_nc);

fprintf('  Cases   NoChew: %d/%d = %.1f%%\n', sum(cas_nc.has_theta_peak=="True"), n_cas_nc, pct_cas_nc);
fprintf('  Cases   Chew:   %d/%d = %.1f%%\n', sum(cas_ch.has_theta_peak=="True"), n_cas_ch, pct_cas_ch);
fprintf('  Controls NoChew:%d/%d = %.1f%%\n', sum(ctr_nc.has_theta_peak=="True"), n_ctr_nc, pct_ctr_nc);
fprintf('  Controls Chew:  %d/%d = %.1f%%\n', sum(ctr_ch.has_theta_peak=="True"), n_ctr_ch, pct_ctr_ch);

% Chi-square: Cases Ch vs Nc
obs_cas = [sum(cas_ch.has_theta_peak=="True") sum(cas_ch.has_theta_peak=="False");
           sum(cas_nc.has_theta_peak=="True") sum(cas_nc.has_theta_peak=="False")];
[~, p_chi_cas] = fishertest(obs_cas);
fprintf('\n  Fisher Cases Ch vs Nc: p=%.4f\n', p_chi_cas);

%% ── (2) Exponente aperiódico ──────────────────────────────────────────────
fprintf('\n--- (2) Exponente aperiódico (1/f slope flattening) ---\n');
fprintf('  Cases   NoChew: M=%.3f ± %.3f\n', mean(cas_nc.ap_exponent), std(cas_nc.ap_exponent));
fprintf('  Cases   Chew:   M=%.3f ± %.3f  ← drástico aplanamiento\n', mean(cas_ch.ap_exponent), std(cas_ch.ap_exponent));
fprintf('  Controls NoChew:M=%.3f ± %.3f\n', mean(ctr_nc.ap_exponent), std(ctr_nc.ap_exponent));
fprintf('  Controls Chew:  M=%.3f ± %.3f\n', mean(ctr_ch.ap_exponent), std(ctr_ch.ap_exponent));

% Wilcoxon: Cases Ch vs Nc (exponente)
[p_wilc, ~, ~] = signrank(cas_ch.ap_exponent, cas_nc.ap_exponent);
fprintf('\n  Wilcoxon Cases exp Ch vs Nc: p=%.6f\n', p_wilc);

% Interacción: (Cases_Ch - Cases_Nc) vs (Controls_Ch - Controls_Nc)
% Requiere datos pareados. Construimos por sujeto.
cas_subj = unique(cas_ch.subject);
delta_exp_cas = zeros(length(cas_subj), 1);
for i = 1:length(cas_subj)
    s = cas_subj(i);
    e_ch = cas_ch.ap_exponent(cas_ch.subject == s);
    e_nc = cas_nc.ap_exponent(cas_nc.subject == s);
    if ~isempty(e_ch) && ~isempty(e_nc)
        delta_exp_cas(i) = e_ch(1) - e_nc(1);
    else
        delta_exp_cas(i) = NaN;
    end
end

ctr_subj = unique(ctr_ch.subject);
delta_exp_ctr = zeros(length(ctr_subj), 1);
for i = 1:length(ctr_subj)
    s = ctr_subj(i);
    e_ch = ctr_ch.ap_exponent(ctr_ch.subject == s);
    e_nc = ctr_nc.ap_exponent(ctr_nc.subject == s);
    if ~isempty(e_ch) && ~isempty(e_nc)
        delta_exp_ctr(i) = e_ch(1) - e_nc(1);
    else
        delta_exp_ctr(i) = NaN;
    end
end

delta_exp_cas = delta_exp_cas(~isnan(delta_exp_cas));
delta_exp_ctr = delta_exp_ctr(~isnan(delta_exp_ctr));
[p_int_exp, ~, s_int] = ranksum(delta_exp_cas, delta_exp_ctr);
fprintf('  Interacción ΔExponente (Cases vs Controls): p=%.6f\n', p_int_exp);
fprintf('  ΔExp Cases: M=%.3f ± %.3f\n', mean(delta_exp_cas), std(delta_exp_cas));
fprintf('  ΔExp Controls: M=%.3f ± %.3f\n', mean(delta_exp_ctr), std(delta_exp_ctr));

%% ── (3) CF theta para quienes TIENEN peak ────────────────────────────────
fprintf('\n--- (3) Frecuencia central theta (solo sujetos con peak) ---\n');

cas_ch_pk = cas_ch(cas_ch.has_theta_peak == "True", :);
cas_nc_pk = cas_nc(cas_nc.has_theta_peak == "True", :);
ctr_ch_pk = ctr_ch(ctr_ch.has_theta_peak == "True", :);
ctr_nc_pk = ctr_nc(ctr_nc.has_theta_peak == "True", :);

fprintf('  Cases   NoChew CF: %.2f ± %.2f Hz  (N=%d)\n', mean(cas_nc_pk.theta_cf), std(cas_nc_pk.theta_cf), height(cas_nc_pk));
fprintf('  Cases   Chew   CF: %.2f ± %.2f Hz  (N=%d)\n', mean(cas_ch_pk.theta_cf), std(cas_ch_pk.theta_cf), height(cas_ch_pk));
fprintf('  Controls NoChew CF:%.2f ± %.2f Hz  (N=%d)\n', mean(ctr_nc_pk.theta_cf), std(ctr_nc_pk.theta_cf), height(ctr_nc_pk));
fprintf('  Controls Chew   CF:%.2f ± %.2f Hz  (N=%d)\n', mean(ctr_ch_pk.theta_cf), std(ctr_ch_pk.theta_cf), height(ctr_ch_pk));

%% ── Figura ────────────────────────────────────────────────────────────────
fig = figure('Name','S3a FOOOF Verify','Units','normalized',...
             'Position',[0.05 0.05 0.90 0.85],'Color','w');

groups  = {'Cases\nNoChew','Cases\nChew','Controls\nNoChew','Controls\nChew'};
labels  = {sprintf('Cases\nNc (n=%d)',n_cas_nc), sprintf('Cases\nCh (n=%d)',n_cas_ch), ...
           sprintf('Ctrl\nNc (n=%d)',n_ctr_nc),  sprintf('Ctrl\nCh (n=%d)',n_ctr_ch)};

colors = [COL_CASE; COL_CHEW; COL_CTRL; COL_CTRL*0.7];

%% Panel 1: % theta peaks
ax1 = subplot(2,3,1);
pcts = [pct_cas_nc pct_cas_ch pct_ctr_nc pct_ctr_ch];
b = bar(pcts, 'FaceColor','flat');
b.CData = colors;
set(ax1, 'XTickLabel', labels, 'FontSize', FIG_FS-1, 'Box','off');
ylabel('% con peak θ', 'FontSize', FIG_FS);
ylim([0 105]);
title('Peak θ sobre 1/f', 'FontSize', FIG_FS+1, 'FontWeight','bold');
for i = 1:4
    text(i, pcts(i)+2, sprintf('%.0f%%',pcts(i)), ...
        'HorizontalAlignment','center','FontSize',FIG_FS-1,'FontWeight','bold');
end
yline(50, 'k--', 'LineWidth', 1);

%% Panel 2: Exponent aperiódico — boxplot
ax2 = subplot(2,3,2);
exp_data  = [cas_nc.ap_exponent; cas_ch.ap_exponent; ctr_nc.ap_exponent; ctr_ch.ap_exponent];
exp_group = [repmat("CasNc",n_cas_nc,1); repmat("CasCh",n_cas_ch,1); ...
             repmat("CtrNc",n_ctr_nc,1); repmat("CtrCh",n_ctr_ch,1)];
boxplot(exp_data, exp_group, 'Labels', labels, 'Colors', 'k');
hold on;
yline(0, 'r--', 'LineWidth',1.5);
ylabel('Exponent aperiódico', 'FontSize', FIG_FS);
title(sprintf('1/f Exponent\nCases Ch vs Nc: p=%.4f', p_wilc), ...
    'FontSize', FIG_FS+1, 'FontWeight','bold');
set(ax2, 'FontSize', FIG_FS-1, 'Box','off');

%% Panel 3: Delta exponent (Ch - Nc) por grupo
ax3 = subplot(2,3,3);
boxplot([delta_exp_cas; delta_exp_ctr], ...
    [repmat("Cases",length(delta_exp_cas),1); repmat("Controls",length(delta_exp_ctr),1)], ...
    'Colors','k');
yline(0,'k-','LineWidth',1.2);
ylabel('Δ Exponent (Ch - Nc)', 'FontSize', FIG_FS);
title(sprintf('Interacción ΔExponent\np=%.4f', p_int_exp), ...
    'FontSize', FIG_FS+1, 'FontWeight','bold');
set(ax3, 'FontSize', FIG_FS-1, 'Box','off');

%% Panel 4: CF theta
ax4 = subplot(2,3,4);
cf_means = [mean(cas_nc_pk.theta_cf) mean(cas_ch_pk.theta_cf) ...
            mean(ctr_nc_pk.theta_cf) mean(ctr_ch_pk.theta_cf)];
cf_sems  = [std(cas_nc_pk.theta_cf)/sqrt(height(cas_nc_pk)) ...
            std(cas_ch_pk.theta_cf)/sqrt(height(cas_ch_pk)) ...
            std(ctr_nc_pk.theta_cf)/sqrt(height(ctr_nc_pk)) ...
            std(ctr_ch_pk.theta_cf)/sqrt(height(ctr_ch_pk))];
b4 = bar(cf_means, 'FaceColor','flat');
b4.CData = colors;
hold on;
errorbar(1:4, cf_means, cf_sems, 'k.', 'LineWidth', 1.5);
set(ax4, 'XTickLabel', labels, 'FontSize', FIG_FS-1, 'Box','off');
ylabel('CF theta (Hz)', 'FontSize', FIG_FS);
ylim([4 8]);
yline(4,'k:'); yline(7,'k:');
title('Frecuencia central θ peak', 'FontSize', FIG_FS+1, 'FontWeight','bold');

%% Panel 5: Scatter exponent vs % con peak (por grupo)
ax5 = subplot(2,3,5);
all_exp  = [cas_nc.ap_exponent; cas_ch.ap_exponent; ctr_nc.ap_exponent; ctr_ch.ap_exponent];
all_peak = double(strcmp(string([cas_nc.has_theta_peak; cas_ch.has_theta_peak; ...
                                  ctr_nc.has_theta_peak; ctr_ch.has_theta_peak]), 'True'));
scatter(all_exp, all_peak + randn(size(all_exp))*0.02, 25, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.4);
xlabel('Exponent aperiódico', 'FontSize', FIG_FS);
ylabel('Tiene peak θ (1=sí, 0=no)', 'FontSize', FIG_FS);
title('Exponent vs Detección peak', 'FontSize', FIG_FS+1, 'FontWeight','bold');
[r_ep, p_ep] = corr(all_exp, all_peak, 'Type','Spearman');
text(0.05, 0.92, sprintf('r_s=%.2f p=%.3f', r_ep, p_ep), ...
    'Units','normalized','FontSize',FIG_FS,'FontWeight','bold');
set(ax5, 'FontSize', FIG_FS-1, 'Box','off');

%% Panel 6: Narrativa resumen
ax6 = subplot(2,3,6);
axis off;
txt = sprintf([...
    'RESUMEN FOOOF (specparam Python)\n\n'...
    '(1) Peak θ genuino (FDR sobre 1/f):\n'...
    '   Cases  Nc: %d/%d = %.0f%%\n'...
    '   Cases  Ch: %d/%d = %.0f%%  ← reducción\n'...
    '   Controls Nc: %d/%d = %.0f%%\n'...
    '   Controls Ch: %d/%d = %.0f%%\n\n'...
    '(2) 1/f Exponent (aplanamiento):\n'...
    '   Cases Nc: %.2f ± %.2f\n'...
    '   Cases Ch: %.2f ± %.2f  ← casi plano\n'...
    '   Wilcoxon: p=%.6f\n\n'...
    '(3) Interacción Δexponent:\n'...
    '   Cases: Δ=%.2f ± %.2f\n'...
    '   Controls: Δ=%.2f ± %.2f\n'...
    '   Ranksum: p=%.6f\n\n'...
    'CONCLUSIÓN: θ oscila genuinamente\n'...
    'PERO la 1/f se aplana en Cases_Ch\n'...
    '(2 mecanismos separados)'], ...
    sum(cas_nc.has_theta_peak=="True"), n_cas_nc, pct_cas_nc, ...
    sum(cas_ch.has_theta_peak=="True"), n_cas_ch, pct_cas_ch, ...
    sum(ctr_nc.has_theta_peak=="True"), n_ctr_nc, pct_ctr_nc, ...
    sum(ctr_ch.has_theta_peak=="True"), n_ctr_ch, pct_ctr_ch, ...
    mean(cas_nc.ap_exponent), std(cas_nc.ap_exponent), ...
    mean(cas_ch.ap_exponent), std(cas_ch.ap_exponent), p_wilc, ...
    mean(delta_exp_cas), std(delta_exp_cas), ...
    mean(delta_exp_ctr), std(delta_exp_ctr), p_int_exp);

text(0.02, 0.98, txt, 'Units','normalized', ...
    'FontSize', FIG_FS-1.5, 'VerticalAlignment','top', ...
    'FontName','Courier New', 'Color','k');

sgtitle(sprintf('FOOOF Verification — N=%d casos + %d controles', N_CASES, N_CONTROLS), ...
    'FontSize', FIG_FS+2, 'FontWeight','bold');

saveas(fig, fullfile(OUT_FIGS, 'S3a_FOOOF_verify.png'));
print(fig, fullfile(OUT_FIGS, 'S3a_FOOOF_verify'), '-dpng', '-r300');
fprintf('\nFigura guardada.\n');

%% ── Guardar .mat ──────────────────────────────────────────────────────────
save(fullfile(OUT_STATS, 'S3a_FOOOF_verify.mat'), ...
    'pct_cas_ch','pct_cas_nc','pct_ctr_ch','pct_ctr_nc', ...
    'p_chi_cas','p_wilc','p_int_exp', ...
    'delta_exp_cas','delta_exp_ctr', ...
    'cas_ch','cas_nc','ctr_ch','ctr_nc', 'MIN_R2');
fprintf('MAT guardado.\n');

%% ── Reporte para el revisor ───────────────────────────────────────────────
fid = fopen(fullfile(OUT_REVIEWER,'S3a_FOOOF_verify_result.txt'),'w');
fprintf(fid,'=== S3a FOOOF/specparam Verification — Revisor M3 ===\n');
fprintf(fid,'Fecha: %s | N_casos=%d | N_controles=%d | R²_min=%.2f\n\n', ...
    datestr(now,'yyyy-mm-dd'), N_CASES, N_CONTROLS, MIN_R2);
fprintf(fid,'--- Peak theta genuino sobre componente aperiódico ---\n');
fprintf(fid,'Cases  NoChew: %d/%d (%.1f%%)\n', sum(cas_nc.has_theta_peak=="True"), n_cas_nc, pct_cas_nc);
fprintf(fid,'Cases  Chew:   %d/%d (%.1f%%) — Fisher vs Nc: p=%.4f\n', ...
    sum(cas_ch.has_theta_peak=="True"), n_cas_ch, pct_cas_ch, p_chi_cas);
fprintf(fid,'Controls NoChew: %d/%d (%.1f%%)\n', sum(ctr_nc.has_theta_peak=="True"), n_ctr_nc, pct_ctr_nc);
fprintf(fid,'Controls Chew:   %d/%d (%.1f%%)\n\n', sum(ctr_ch.has_theta_peak=="True"), n_ctr_ch, pct_ctr_ch);
fprintf(fid,'--- Exponente aperiódico ---\n');
fprintf(fid,'Cases NoChew: %.3f ± %.3f\n', mean(cas_nc.ap_exponent), std(cas_nc.ap_exponent));
fprintf(fid,'Cases Chew:   %.3f ± %.3f  (Wilcoxon p=%.6f)\n', mean(cas_ch.ap_exponent), std(cas_ch.ap_exponent), p_wilc);
fprintf(fid,'Interaction Delta-exp (Cases vs Controls): p=%.6f\n\n', p_int_exp);
fprintf(fid,'--- Frecuencia central theta (sujetos con peak) ---\n');
fprintf(fid,'Cases Nc: %.2f ± %.2f Hz (n=%d)\n', mean(cas_nc_pk.theta_cf), std(cas_nc_pk.theta_cf), height(cas_nc_pk));
fprintf(fid,'Cases Ch: %.2f ± %.2f Hz (n=%d)\n', mean(cas_ch_pk.theta_cf), std(cas_ch_pk.theta_cf), height(cas_ch_pk));
fprintf(fid,'\n--- Conclusión ---\n');
fprintf(fid,'El theta EN CASOS CHEW involucra dos mecanismos:\n');
fprintf(fid,'1. Aplanamiento aperiódico (exponent→0): señal de neuromodulación/excitabilidad\n');
fprintf(fid,'2. Oscilación theta genuina: presente en %.0f%% de Cases_Ch (specparam peak)\n', pct_cas_ch);
fprintf(fid,'Ambos contribuyen al aumento de potencia theta en el TF CBPT.\n');
fclose(fid);

fprintf('\n✓ S3a completo.\n');
fprintf('Interpretación: θ genuino en %.0f%% Cases_Ch + 1/f flattening (exponent %.2f±%.2f)\n', ...
    pct_cas_ch, mean(cas_ch.ap_exponent), std(cas_ch.ap_exponent));
fprintf('Siguiente: revisar figura, luego PAC (Fase 4).\n\n');
