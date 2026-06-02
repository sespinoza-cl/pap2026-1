%% ============================================================================
%  05_FOOOF_LME.m — Modelo lineal de efectos mixtos del exponente aperiódico
%  ----------------------------------------------------------------------------
%  Reemplaza el LME ad-hoc que antes se corría en el Live Editor (solo quedaba
%  el .txt Stats_FOOOF_Exp2.txt, sobre N=46 → INCORRECTO).
%
%  Modelo:  Exponent ~ 1 + Group*Condition + (1 | ID)
%  Muestra: 30 casos + 15 controles = 45 sujetos  →  90 observaciones (DF=86).
%
%  Lee el exponente desde FOOOF_Workspace.mat (mismo que usa 04_FOOOF_Figuras.m):
%    GR.Cases.exp_Ch/exp_Nc   (vector nCas×1)
%    GR.Controls.exp_Ch/exp_Nc(vector nCtr×1)
%
%  El efecto within-group (Casos Ch vs Nc, signrank) se reporta en
%  04_FOOOF_Figuras.m; aquí va la INTERACCIÓN Grupo×Condición.
%  ============================================================================
clear; clc; close all;

P = S0_paths();
assert(exist(P.file_fooof_ws,'file')==2, 'No se encontró FOOOF_Workspace.mat:\n  %s', P.file_fooof_ws);
W  = load(P.file_fooof_ws);
GR = W.GR;

exp_cas_ch = GR.Cases.exp_Ch(:);     exp_cas_nc = GR.Cases.exp_Nc(:);
exp_ctr_ch = GR.Controls.exp_Ch(:);  exp_ctr_nc = GR.Controls.exp_Nc(:);

nCas = numel(exp_cas_ch);
nCtr = numel(exp_ctr_ch);
nTot = nCas + nCtr;

fprintf('>>> FOOOF LME — N casos=%d | N controles=%d | total=%d | obs=%d\n', ...
    nCas, nCtr, nTot, 2*nTot);
if nCas ~= 30 || nCtr ~= 15
    warning(['Muestra esperada = 30 casos + 15 controles = 45 (90 obs). ' ...
             'El FOOOF_Workspace tiene %d + %d. Regenerar el workspace con la ' ...
             'muestra final (excluir E3S3, E3S5) antes de reportar el LME.'], nCas, nCtr);
end

% ── Tabla larga: cada sujeto aporta 2 filas (Nc, Ch) ───────────────────────
Exponent  = [exp_cas_nc; exp_cas_ch; exp_ctr_nc; exp_ctr_ch];
Group     = categorical([ repmat({'Cases'},   2*nCas, 1); ...
                          repmat({'Controls'},2*nCtr, 1) ]);
Condition = categorical([ repmat({'Nc'}, nCas,1); repmat({'Ch'}, nCas,1); ...
                          repmat({'Nc'}, nCtr,1); repmat({'Ch'}, nCtr,1) ]);
id_cas = (1:nCas)';   id_ctr = (nCas + (1:nCtr))';
ID     = categorical([ id_cas; id_cas; id_ctr; id_ctr ]);

T = table(Exponent, Group, Condition, ID);

% ── LME ────────────────────────────────────────────────────────────────────
lme = fitlme(T, 'Exponent ~ 1 + Group*Condition + (1 | ID)');
disp(lme);

% ── Reporte ─────────────────────────────────────────────────────────────────
f_rep = fullfile(P.fig03, 'Reporte_FOOOF_LME.txt');
fid = fopen(f_rep, 'w');
fprintf(fid, '================================================================\n');
fprintf(fid, '  FOOOF — LME exponente aperiódico (05_FOOOF_LME.m)\n');
fprintf(fid, '  Generado: %s\n', datestr(now));
fprintf(fid, '  Modelo: Exponent ~ 1 + Group*Condition + (1|ID)\n');
fprintf(fid, '  N = %d casos + %d controles = %d  (obs=%d)\n', nCas, nCtr, nTot, 2*nTot);
fprintf(fid, '================================================================\n\n');
co = lme.Coefficients;
fprintf(fid, '%-32s %10s %10s %8s %5s %12s\n', 'Term','Estimate','SE','t','DF','p');
for i = 1:height(co)
    fprintf(fid, '%-32s %10.4f %10.4f %8.3f %5d %12.3e\n', ...
        char(co.Name{i}), co.Estimate(i), co.SE(i), co.tStat(i), co.DF(i), co.pValue(i));
end
fprintf(fid, '\nNota: el término de interacción Group:Condition cuantifica la diferencia\n');
fprintf(fid, 'del efecto de la masticación (Ch vs Nc) entre Casos y Controles.\n');
fclose(fid);
fprintf('>>> Reporte LME: %s\n', f_rep);
