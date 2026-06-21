# REPRODUCE — Paper2 V1 Final (mapa de reproducibilidad)
Cómo regenerar TODOS los resultados/figuras desde los datos. Config única: `S0_config.m` (MATLAB) /
`code/_figstyle.py` (estilo de figuras). Diseño: 46 adultos sanos (31 casos mastican + 15 controles).
Marco final: **θ-PAC cortical genuino + aperiódico=excitabilidad; dominancia β = leakage EMG**.

## Requisitos
- **MATLAB R2025b** + Statistics/Signal Processing Toolboxes; **EEGLAB** en `D:\EEGLAB` (solo para .set/topoplots).
- **Python 3.10** + `numpy scipy h5py matplotlib mne fooof statsmodels pandas pillow` (vía `uv`/pip).
- Datos: `.mat` canónicos en `outputs/stats/` (incluidos, ~6.5 MB). Crudos EEG/EMG `.set` en `Data_PAC/`
  y `D:\Exp1\Exp1` (Paper 1) — NO en el repo (pesados/privados; pedir al autor).

## Orden de ejecución
### A) Pipeline canónico (MATLAB) — produce los .mat de `outputs/stats/`
| Paso | Script | Produce |
|---|---|---|
| 0 | `S0_config.m` | parámetros (bandas, ROI-18, ventanas, seed=42) |
| 1 | `code/S1_behavior_analysis.m` | `v1_S1_behavior_stats.mat` (RT/IES/d, LME) |
| 1b| `code/S0c_chew_engagement.m` | `chew_engagement_check.mat` (manip-check masetero) |
| 2 | `code/S2a/S2b/S2c/S2d_*.m` | CBPT + TF group + ROI-18 (`S2*`, `ROI_canonical.mat`) |
| 3 | `code/S3a/S3c_FOOOF_*.m` (+ `fooof_fit.py`) | `FOOOF_Workspace_V1.mat` (χ, pico θ) |
| 4 | `code/S4b_PAC_ROI.m` (`compute_pac_*.m`, `aaft.m`, `emg_bilateral.m`) | `v1_S4b_PAC_ROI.mat` (zMI doble null) |
| 4g/h | `code/S4g/S4h_*.m` (`cluster_stat_2d.m`) | comodulogramas `S4g/S4h_*.mat` |
| 5 | `code/S6_artifact_controls.m`, `S6b_A3_muscle.m` | `S6_*`, `S6b_*` (A1–A6) |
| 6 | `code/S7_dose_response.m`, `S8_controls_PAC.m` | `S7_*`, `S8_*` (dosis, control negativo) |

### B) Figuras canónicas (Python + MATLAB) — leen los .mat, escriben `outputs/figures/`
| Script | Figuras |
|---|---|
| `code/figA_py.py` | conducta (RT/IES), manip-check, PAC zMI (doble null, 2×2, bandas), roses θ/β |
| `code/figA2_py.py` | TF maps (casos/ctrl/interacción), FOOOF (PSD/exponente/pico), comodulograma |
| `code/figA_topo.m` | topo θ interacción + por grupo (EEGLAB) |
| `code/figB_py.py` | suplementarios: spatial, masetero×Δθ, dosis, MI-ventana, comodulograma 2D |
| estilo | `code/_figstyle.py` (box+scatter, symlog, Okabe-Ito, cuadradas PNG@300) |

### C) Análisis de revisión (Python) — NUEVOS, responden al panel/vacíos
| Script | Resultado | .mat/.csv |
|---|---|---|
| `code/rev_dissociation.py` / `rev_dissociation2.py` | Δneural×conducta entre-sujetos = n.s. (todas las métricas) | — |
| `code/rev_singletrial.py` | **puente trial-level**: θ-late×RT (LMM −41.8 ms, p=0.0007) | `rev_singletrial.csv` |
| `code/rev_singletrial_fig.py` | Fig.6 (quintiles + especificidad ventana) | — |
| `code/rev_beta_cmc.m` (MATLAB) | CMC β EEG-EMG (magnitud al piso) | `rev_beta_cmc.mat` |
| `code/rev_cmc_topography.m` (MATLAB) | **topografía CMC β = temporal/leakage** (no central) | `rev_cmc_topography.mat` |
| `code/rev_paper1_aperiodic.py` | **cross-study**: aperiódico aplana offline (Paper 1, Δ−0.09, p=0.005) | `rev_paper1_aperiodic.npz` |
| `code/fig_suppl_offline_aperiodic.py` | Suppl: aperiódico offline vs online | — |

## Resultados clave → dónde verificarlos
Números canónicos consolidados: `outputs/reviewer/RESULTS_CANONICAL.md`, `CLAIMS_SUMMARY.md`,
`RESULTS_TO_PASTE.txt`. Análisis de revisión: `REVISION_ANALYSES.md`, `REVISION_SUMMARY.md`,
`NARRATIVE_online_vs_offline.md`, `LIT_SUPPORT.md`. Figuras: `outputs/figures/FIGURE_INDEX.md`.

## Notas
- `parfor` en PAC (S4b): usar `eeglab nogui`, NO `addpath(genpath(eeglab))` dentro de parfor.
- Forzar UTF-8 en Python (`PYTHONUTF8=1`) por nombres con diacríticos.
- Seeds fijas (42) para jitter/surrogados → reproducible.
