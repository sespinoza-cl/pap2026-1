%% run_figures.m — Reproduce all paper figures from pre-computed workspaces
%
%   USAGE
%     1. Open MATLAB (R2020b or later recommended).
%     2. Set the Current Folder to the repo root (where this file lives).
%     3. Run:  >> run_figures
%
%   OUTPUT
%     All figures are saved as 300-dpi PNG files in outputs/.
%     See the folder map below for each figure's location.
%
%   REQUIREMENTS
%     - MATLAB R2020b+  (tested on R2024b)
%     - Statistics and Machine Learning Toolbox
%     - Signal Processing Toolbox
%     - No EEGLAB needed  (pre-computed workspaces are included)
%
%   WHAT RUNS (repo mode — no raw EEG required)
%     Script                  Output folder               Paper figure
%     S01_Conducta            outputs/Figure01_Behavior   Fig. 1 A-B
%     S04_FOOOF_Figuras       outputs/Figure03_FOOOF      Fig. 1 C-D + Supp S3
%     S05_FOOOF_LME           outputs/Figure03_FOOOF      Reporte_FOOOF_LME.txt
%     S07_Supplementary       outputs/Supplementary       Supp Figs S2-S4
%     plots/P_PAC_Panels      outputs/figures             Fig. 3 A-B-D
%     plots/P_Rayleigh_Panels outputs/figures             Supp Fig. S4 (polar)
%     plots/P_Supp_Panels     outputs/figures             Supp Fig. S3
%
%   WHAT REQUIRES RAW EEG (full pipeline — see README.md)
%     S02_TF_Correlaciones  — 6-GB TF matrices
%     S03_TF_ClusterPerm    — epoched EEG (.set files)
%     S06_PAC               — continuous cleaned EEG (.set files)
%
% ─────────────────────────────────────────────────────────────────────────

clear; clc; close all;
fprintf('=== pap2026-1: Reproducing paper figures (%s) ===\n\n', datestr(now));

% Add repo root and plots subfolder to path
repo = fileparts(mfilename('fullpath'));
addpath(repo, fullfile(repo,'plots'));

% ── 1. Behavioural results (Fig. 1 A-B) ──────────────────────────────────
fprintf('[1/6] S01_Conducta — behavioural figures...\n');
run(fullfile(repo, 'S01_Conducta.m'));
close all;

% ── 2. FOOOF spectral decomposition figures (Fig. 1 C-D, Supp S3) ────────
fprintf('[2/6] S04_FOOOF_Figuras — FOOOF spectra...\n');
run(fullfile(repo, 'S04_FOOOF_Figuras.m'));
close all;

% ── 3. FOOOF LME statistics (Reporte_FOOOF_LME.txt) ──────────────────────
fprintf('[3/6] S05_FOOOF_LME — mixed-effects model...\n');
run(fullfile(repo, 'S05_FOOOF_LME.m'));
close all;

% ── 4. PAC panels (Fig. 3 A-B-D) ─────────────────────────────────────────
fprintf('[4/6] P_PAC_Panels — PAC heatmap, scatter, violin...\n');
run(fullfile(repo, 'plots', 'P_PAC_Panels.m'));
close all;

% ── 5. Rayleigh preferred-phase polar plots (Supp S4) ────────────────────
fprintf('[5/6] P_Rayleigh_Panels — masticatory phase preference...\n');
run(fullfile(repo, 'plots', 'P_Rayleigh_Panels.m'));
close all;

% ── 6. Supplementary panels + mediation + Steiger (Supp S2-S3) ───────────
fprintf('[6/6] S07_Supplementary — supplementary figures...\n');
run(fullfile(repo, 'S07_Supplementary.m'));
close all;

% ── Done ──────────────────────────────────────────────────────────────────
fprintf('\n=== All figures saved to outputs/ ===\n');
fprintf('  Figure 1 (Behavior + FOOOF) : outputs/Figure01_Behavior/\n');
fprintf('                                outputs/Figure03_FOOOF/\n');
fprintf('  Figure 3 (PAC)              : outputs/figures/  (P_Fig3*.png)\n');
fprintf('  Supplementary               : outputs/Supplementary/\n');
fprintf('  Rayleigh polar plots        : outputs/figures/  (P_Rayleigh*.png)\n\n');
fprintf('Statistical report (LME)      : outputs/Figure03_FOOOF/Reporte_FOOOF_LME.txt\n\n');
fprintf('Note: Figures 1C-D and 3A require compositing in Inkscape\n');
fprintf('(individual panels provided; see README.md).\n');
