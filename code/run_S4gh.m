% run_S4gh.m — comodulograma descriptivo (S4g) + stats 2D (S4h)
addpath('D:\EEGLAB'); evalc('eeglab nogui'); close all force
THIS = fileparts(mfilename('fullpath'));
addpath(THIS); addpath(fileparts(THIS));
fprintf('\n########## S4g (descriptivo) ##########\n');
run(fullfile(THIS,'S4g_comodulogram.m'));
fprintf('\n########## S4h (stats 2D) ##########\n');
run(fullfile(THIS,'S4h_comodulogram_stats.m'));
fprintf('\n########## F6 COMPLETO ##########\n');
