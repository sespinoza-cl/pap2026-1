%% S0_config.m — Configuración maestra Analysis_V1_Final
% Pipeline P2V1 Revisión Revisor 1 | N=31 casos + 15 controles
% Fuente única de verdad: todos los scripts en code/ hacen run(S0_config.m)
% Actualizado: 2026-06-17

%% ── Rutas absolutas ──────────────────────────────────────────────────────
ROOT_FINAL = fileparts(mfilename('fullpath'));   % .../Analysis_V1_Final/
ROOT_P2V1  = fileparts(ROOT_FINAL);             % .../P2V1/

% --- Datos de entrada (grandes, no copiados — se referencian) ---
DATA_PAC       = fullfile(ROOT_P2V1, 'Data_PAC');              % clean_emg.set
DIR_EPOCHS     = fullfile(DATA_PAC, 'Epochs');                  % *_ep.set
FILE_CHEW      = fullfile(ROOT_FINAL, 'data', 'computed', 'chew_metrics.mat');
DATA_TFRAW     = fullfile(ROOT_P2V1, 'Analysis_v1', 'outputs', 'stats'); % Group_*_tfraw.mat
FILE_FULLSPACE = fullfile(ROOT_P2V1, 'Analysis_paper', 'outputs', 'stats', 'S2_TF_fullspace_stats.mat');

% --- Datos computados base (copiados, pequeños) ---
DATA_COMPUTED  = fullfile(ROOT_FINAL, 'data', 'computed');

% --- Archivos v1_ (N=31) como fuente de verdad ---
FILE_BEH    = fullfile(DATA_COMPUTED, 'v1_S1_behavior_stats.mat');
FILE_TF     = fullfile(DATA_COMPUTED, 'v1_S2_TF_data.mat');
FILE_FOOOF  = fullfile(DATA_COMPUTED, 'v1_S3_specparam_stats.mat');
FILE_LME    = fullfile(DATA_COMPUTED, 'v1_S3_LME_stats.mat');
FILE_PAC    = fullfile(DATA_COMPUTED, 'v1_S4_PAC_stats.mat');
FILE_CSV    = fullfile(DATA_COMPUTED, 'v1_S3_specparam_results.csv');

% --- Outputs de este pipeline ---
OUT_STATS    = fullfile(ROOT_FINAL, 'outputs', 'stats');
OUT_FIGS     = fullfile(ROOT_FINAL, 'outputs', 'figures');
OUT_REVIEWER = fullfile(ROOT_FINAL, 'outputs', 'reviewer');

for d_ = {OUT_STATS, OUT_FIGS, OUT_REVIEWER}
    if ~exist(d_{1},'dir'), mkdir(d_{1}); end
end
clear d_

addpath(fullfile(ROOT_FINAL, 'code'));
addpath(ROOT_FINAL);

%% ── Sujetos (N=31 + 15) ─────────────────────────────────────────────────
% N=31 casos: E3S1-S33 excluyendo E3S3 y E3S5.
% E3S12 INCLUIDO (criterio original del manuscrito V1)
CASES = {'E3S1','E3S2','E3S4','E3S6','E3S7','E3S8','E3S9','E3S10',...
         'E3S11','E3S12','E3S13','E3S14','E3S15','E3S16','E3S17','E3S18',...
         'E3S19','E3S20','E3S21','E3S22','E3S23','E3S24','E3S25','E3S26',...
         'E3S27','E3S28','E3S29','E3S30','E3S31','E3S32','E3S33'};
CONTROLS = {'E3C1','E3C2','E3C3','E3C4','E3C5','E3C6','E3C7','E3C8',...
            'E3C9','E3C10','E3C11','E3C12','E3C13','E3C14','E3C15'};
N_CASES    = numel(CASES);    % 31
N_CONTROLS = numel(CONTROLS); % 15

%% ── Canales ──────────────────────────────────────────────────────────────
EEG_N      = 64;
EMG_CHANS  = [65 66];   % bipolar masétero L e I en los clean_emg.set
EMG_CHANS_RAW = [65 67]; % canales EMG en los raw .set (superior bilateral)

% ROI a priori (solo referencia histórica — NO usar en análisis nuevos)
ROI_APRIORI = {'F1','F2','FC1','FC2'};

% ROI CANÓNICO — fuente única de verdad para todos los análisis (S4b, S4g, etc.)
% 18 electrodos: interacción Cases×Cond FDR<0.05 sobre theta_topo (Morlet NPL dB,
% WIN_LATE) en v1_S2_TF_data.mat. Derivado/reproducido por code/S2d_theta_ROI.m
% → outputs/stats/ROI_canonical.mat (verificado 2026-06-19: reproduce este set exacto).
% Robustez: método de épocas (eegfilt+Hilbert²) da 13/18 (subconjunto frontal). (fix I4)
ROI_CBPT = {'Fp1','AF7','AF3','F1','F3','F5','F7','FC1','Fpz',...
             'Fp2','AF8','AF4','AFz','Fz','F2','F4','FC2','FCz'};

%% ── Suffijos de archivos ─────────────────────────────────────────────────
EMG_SUFFIX    = '_Ch_clean_emg.set';
EMG_SUFFIX_NC = '_Nc_clean_emg.set';
EP_SUFFIX     = '_Ch_ep.set';
EP_SUFFIX_NC  = '_Nc_ep.set';

%% ── Bandas de frecuencia (Hz) ────────────────────────────────────────────
% Theta 4-7 Hz con banda de GUARDA intencional 7-8 Hz (theta/alfa).
% Verificado empíricamente (S2d_theta_ROI.m, 2026-06-19): ensanchar a 4-8 Hz
% DILUYE la interacción Casos×Cond en sitios midline FM-theta (pierde Fz/FCz/FC1/FC2),
% porque 7-8 Hz aporta actividad de borde-alfa que no discrimina grupos (interacción
% alfa ≈ 0 electrodos). 4-7 preserva el ROI midline cognitivo → blinda anti-artefacto.
% Pico FOOOF empírico ≈6.5 Hz (central en 4-7).
BAND_THETA = [4   7];
BAND_ALPHA = [8  13];
BAND_BETA  = [13 30];
BAND_EMG   = [0.5  2];    % rango de fase masticatoria
BAND_MUSC  = [30  40];    % proxy de contaminación muscular

%% ── Ventanas temporales (ms) ─────────────────────────────────────────────
WIN_EARLY    = [100   900];
WIN_LATE     = [900  1300];
WIN_BASE     = [-500 -100];
WIN_EPOCH    = [-1000 1500];
% Ventana de análisis estadístico (trim para CBPT y gráficas)
% Trial: estímulo 200ms + ISI 2000ms. TF computado en [-2000,2496]ms para
% margen de borde wavelet. Se restringe a [-200,1500]ms para el análisis.
WIN_ANALYSIS = [-200 1500];

%% ── Parámetros análisis ──────────────────────────────────────────────────
FS         = 500;       % frecuencia de muestreo (post-downsample)
N_PERM     = 5000;
ALPHA_VOXEL = 0.01;    % threshold voxel para CBPT
ALPHA_CLUST = 0.05;    % threshold cluster
MI_BINS    = 18;
N_SURR     = 500;
RNG_SEED   = 42;

%% ── TF ───────────────────────────────────────────────────────────────────
FREX       = linspace(1, 40, 400);
N_FREX     = numel(FREX);
N_CYCLES   = logspace(log10(3), log10(10), N_FREX);
BASE_WIN   = [-500 -100];

%% ── Colores (Paleta Okabe-Ito — proyecto P2V1) ───────────────────────────
COL_CASE  = [0.000, 0.620, 0.451];   % #009E73  bluish-green  (Cases)
COL_CTRL  = [0.000, 0.447, 0.698];   % #0072B2  blue          (Controls)
COL_NC    = [0.533, 0.533, 0.533];   % #888888  gris neutro   (NoChew)
COL_THETA = [0.835, 0.369, 0.000];   % #D55E00  vermillion    (banda θ)
COL_ALPHA = [0.337, 0.706, 0.914];   % #56B4E9  sky blue      (banda α)
COL_BETA  = [0.800, 0.475, 0.655];   % #CC79A7  reddish-purple(banda β)
COL_CHEW  = COL_CASE;                % alias
COL_SIG   = [0.9  0.2  0.2];        % rojo para marcadores de sig

%% ── Figura ───────────────────────────────────────────────────────────────
FIG_DPI  = 300;
FIG_FONT = 'Arial';
FIG_FS   = 10;
FIG_LW   = 1.5;

fprintf('[S0_config] P2V1 Final | N=%d casos + %d controles | ROI_CBPT=%d elecs (%s...) | %s\n', ...
    N_CASES, N_CONTROLS, numel(ROI_CBPT), strjoin(ROI_CBPT(1:3),'+'), datestr(now,'yyyy-mm-dd HH:MM'));
