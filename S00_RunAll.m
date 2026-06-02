%% ============================================================================
%  00_RUN_ALL.m — Ejecuta TODO el pipeline del Paper 2 en el orden correcto
%  ----------------------------------------------------------------------------
%  Uso:  abrir en MATLAB y pulsar Run (o  >> run('S00_RunAll.m')  ).
%
%  Notas:
%   • Cada script Sxx hace 'clear' al inicio, por eso aquí NO se usan variables
%     de loop entre llamadas (secuencia plana de run()). El addpath y el parpool
%     sobreviven a los 'clear'.
%   • El pool de workers se crea UNA vez y lo reutilizan 03 (cluster-perm) y 06 (PAC).
%   • Tiempos por script: cada Sxx imprime su propio tiempo interno.
%   • Orden con dependencias: 02 y 06 generan los workspaces que consume 07.
%  ============================================================================
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
cd(here); addpath(here);              % asegura que S0_paths.m y los scripts se encuentren

% Pool de workers — Ryzen 5900X (12 núcleos físicos)
if isempty(gcp('nocreate'))
    parpool('local', min(12, feature('numcores')));
end

fprintf('\n############ PIPELINE PAPER 2 — INICIO  %s ############\n', datestr(now));

run('S01_Conducta.m');          fprintf('\n>>>>> [1/9] S01_Conducta OK <<<<<\n');
run('S02_TF_Correlaciones.m');  fprintf('\n>>>>> [2/9] S02_TF_Correlaciones OK <<<<<\n');
run('S03_TF_ClusterPerm.m');    fprintf('\n>>>>> [3/9] S03_TF_ClusterPerm OK <<<<<\n');
run('S04_FOOOF_Figuras.m');     fprintf('\n>>>>> [4/9] S04_FOOOF_Figuras OK <<<<<\n');
run('S05_FOOOF_LME.m');         fprintf('\n>>>>> [5/9] S05_FOOOF_LME OK <<<<<\n');
run('S06_PAC.m');               fprintf('\n>>>>> [6/9] S06_PAC OK <<<<<\n');
% S06b_PAC_EEG (PAC EEG-EEG θ-γ) REMOVIDO del paper — análisis γ artefactual / bandas inconsistentes.
% S07 detecta la ausencia de PAC_EEG_Workspace.mat y salta sus secciones D/E automáticamente.
run('S07_Supplementary.m');     fprintf('\n>>>>> [7/9] S07_Supplementary OK <<<<<\n');
run('S08_ChewFreq.m');          fprintf('\n>>>>> [8/9] S08_ChewFreq OK <<<<<\n');
run('S09_Figuras_Final.m');     fprintf('\n>>>>> [9/9] S09_Figuras_Final OK <<<<<\n');

fprintf('\n############ PIPELINE PAPER 2 — FIN  %s ############\n', datestr(now));
fprintf('Salidas en: %s\n', fullfile(here,'outputs'));
