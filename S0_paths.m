function P = S0_paths()
%% S0_paths — Central configuration for the pap2026-1 repository
%
%   P = S0_paths()
%
%   Centralises all paths and methodological parameters.
%   Called at the top of every script.
%
%   ── MODES ─────────────────────────────────────────────────────────────────
%   REPO MODE (default, DATA_ROOT = ''):
%     Pre-computed workspaces in data/workspaces/ and outputs/.
%     run_figures.m uses this mode — no raw EEG needed.
%
%   FULL PIPELINE MODE (DATA_ROOT set below):
%     S02, S03, S06 read raw EEG/TF data from DATA_ROOT.
%     Requires EEGLAB and the EEG dataset (available on request).

% ╔══════════════════════════════════════════════════════════════════════╗
% ║  Set DATA_ROOT to your local raw-EEG folder to run the full pipeline ║
% ║  Leave empty ('') to reproduce figures from pre-computed workspaces   ║
% ╚══════════════════════════════════════════════════════════════════════╝
DATA_ROOT     = '';          % e.g. 'C:\Data\Exp2'  or  '/data/Exp2'
P.eeglab_path = '';          % e.g. 'C:\EEGLAB'

% ── Repo root and output root ──────────────────────────────────────────────
P.dir_pip  = fileparts(mfilename('fullpath'));   % repo root
P.dir_out  = fullfile(P.dir_pip, 'outputs');
P.dir_data = fullfile(P.dir_pip, 'data');
P.dir_ws   = fullfile(P.dir_pip, 'data', 'workspaces');

% ── Output subdirectories ─────────────────────────────────────────────────
P.fig01     = fullfile(P.dir_out, 'Figure01_Behavior');
P.fig02     = fullfile(P.dir_out, 'Figure02_TF');
P.fig02b    = fullfile(P.dir_out, 'Figure02b_TF_Correlaciones');
P.fig03     = fullfile(P.dir_out, 'Figure03_FOOOF');
P.fig04     = fullfile(P.dir_out, 'Figure04_PAC');
P.fig05     = fullfile(P.dir_out, 'Figure05_ChewFreq');
P.dir_supp  = fullfile(P.dir_out, 'Supplementary');
P.dir_paper = P.dir_out;

% ── Light data (always in repo) ───────────────────────────────────────────
P.file_beh       = fullfile(P.dir_data, 'data_beh_tb_45.mat');
P.file_incluidos = fullfile(P.dir_data, 'incluidos45.mat');
P.dir_lists      = P.dir_data;

% ── Pre-computed workspaces ───────────────────────────────────────────────
%   FOOOF decomposition (used by S04_FOOOF_Figuras, S05_FOOOF_LME)
P.file_fooof_ws  = fullfile(P.dir_ws, 'FOOOF_Workspace.mat');
P.file_fooof_txt = fullfile(P.fig03,  'Reporte_FOOOF_LME.txt');

%   PAC and TF workspaces live in their natural outputs/ subdirectories
%   (outputs/Figure04_PAC/ and outputs/Figure02b_TF_Correlaciones/).
%   S07_Supplementary and the plot scripts read them from P.fig04 / P.fig02b.

% ── Raw EEG paths (only used by S02, S03, S06 — full pipeline mode) ───────
if ~isempty(DATA_ROOT)
    P.dir_tf    = fullfile(DATA_ROOT, 'EEG', 'TF');
    P.dir_epoch = fullfile(DATA_ROOT, 'EEG', 'Epochs');
    P.dir_cont  = fullfile(DATA_ROOT, 'EEG', 'Prepro', 'S5_Final');
    P.file_chew = fullfile(DATA_ROOT, 'EEG', 'ChewFreq', 'chew_metrics.mat');
    P.dir_pac   = fullfile(DATA_ROOT, 'EEG', 'PAC', 'Continuous');
else
    P.dir_tf    = '';
    P.dir_epoch = '';
    P.dir_cont  = '';
    P.file_chew = '';
    P.dir_pac   = '';
end

% ╔══════════════════════════════════════════════════════════════════════╗
% ║  METHODOLOGICAL PARAMETERS — single definition for the entire paper  ║
% ╚══════════════════════════════════════════════════════════════════════╝

% ROI: F1 (ch12), FC1 (ch16)
P.roi_ch     = [12, 16];
P.roi_labels = {'F1', 'FC1'};

% Canonical frequency bands (literature convention)
P.bands    = {'Theta', 'Alpha', 'Beta'};
P.bands_hz = {[4 7], [8 12], [13 30]};   % θ 4-7 | α 8-12 | β 13-30 Hz

% Post-stimulus windows (ms). Column order: Early=1, Late=2, Mid=3
P.wins    = {'Early', 'Late', 'Mid'};
P.wins_ms = {[0 300], [300 900], [200 700]};

% Excluded participants (insufficient EMG/EEG quality)
P.excluded = {'E3S3', 'E3S5'};

% Consistent colour palette across all figures
P.clr_cas = [0.80 0.22 0.22];   % Cases
P.clr_ctr = [0.22 0.45 0.72];   % Controls
P.clr_ch  = [0.20 0.60 0.45];   % Chew condition
P.clr_nc  = [0.60 0.60 0.60];   % No-Chew condition

% ── Create output directories if they do not exist ────────────────────────
dirs = {P.dir_out, P.fig01, P.fig02, P.fig02b, P.fig03, ...
        P.fig04, P.fig05, P.dir_supp};
for k = 1:numel(dirs)
    if ~exist(dirs{k}, 'dir'); mkdir(dirs{k}); end
end

end
