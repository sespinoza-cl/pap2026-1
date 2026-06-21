% run_S4b.m — runner: añade EEGLAB y ejecuta S4b_PAC_ROI.
% Variables esperadas desde -batch: SMOKE (true/false). Default full.
addpath('D:\EEGLAB');
evalc('eeglab nogui');            % inicializa paths sin GUI, silencioso
close all force
THIS = fileparts(mfilename('fullpath'));
addpath(THIS); addpath(fileparts(THIS));
if ~exist('SMOKE','var'), SMOKE = false; end
run(fullfile(THIS,'S4b_PAC_ROI.m'));
