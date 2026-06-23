%% run_pending_matlab.m — Correr todos los análisis MATLAB pendientes P2V1 R1
% Orden: (1) Gap1 Δχ×EMG, (2) M7 neural baseline, (3) trial counts
%
% Prerrequisitos:
%   - S6b_A3_muscle.mat ya existe (emg_rms computado)
%   - v1_S2_TF_data.mat ya existe
%   - v1_S3_specparam_results.csv ya existe
%   - EEGLAB en el path (para Strial_counts.m Intento 2)
%
% Uso:
%   cd('C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1\Analysis_V1_Final')
%   run S0_config
%   run code/run_pending_matlab

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

fprintf('\n======================================================\n');
fprintf(' run_pending_matlab.m — P2V1 R1 análisis pendientes\n');
fprintf(' %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('======================================================\n');

%% 1. Gap #1: Δχ × EMG RMS (no requiere EEGLAB)
fprintf('\n[1/3] SGap1_chi_emg ...\n');
run(fullfile(ROOT_FINAL, 'code', 'SGap1_chi_emg.m'));

%% 2. M7: Equivalencia neural baseline (no requiere EEGLAB)
fprintf('\n[2/3] SM7_baseline_neural ...\n');
run(fullfile(ROOT_FINAL, 'code', 'SM7_baseline_neural.m'));

%% 3. Trial counts (requiere EEGLAB para Intento 2)
fprintf('\n[3/3] Strial_counts ...\n');
run(fullfile(ROOT_FINAL, 'code', 'Strial_counts.m'));

fprintf('\n======================================================\n');
fprintf(' COMPLETADO. Resultados en:\n');
fprintf('   %s\n', OUT_REVIEWER);
fprintf('   SGap1_chi_emg_result.txt\n');
fprintf('   SM7_baseline_neural_result.txt\n');
fprintf('   Strial_counts_result.txt\n');
fprintf('======================================================\n');
