# pap2026-1 — Code and Data for:

**"Rhythmic mastication entrains frontocentral β-band amplitude via corticomuscular phase–amplitude coupling during working memory"**

Espinoza S., Moraga-Espinoza D., Morreal-Ortega L., El-Deredy W. (2025). *Scientific Reports.*  
DOI: 
Repository: [github.com/sespinoza-cl/pap2026-1](https://github.com/sespinoza-cl/pap2026-1)

---

## Quick start — reproduce all paper figures

> **Requirements:** MATLAB R2020b or later · Statistics and Machine Learning Toolbox · Signal Processing Toolbox · No EEGLAB needed.

```matlab
% 1. Clone or download this repository
% 2. Open MATLAB, navigate to the repo root
% 3. Run:
run_figures
```

That's it. All figures are saved to `outputs/` as 300-dpi PNG files. No raw EEG data required — all necessary pre-computed workspaces are included in this repository (~1 MB total).

---

## Repository structure

```
pap2026-1/
├── run_figures.m              ← ENTRY POINT — run this
├── S0_paths.m                 ← central config (paths + parameters)
│
├── code — analysis pipeline (S01–S09)
│   ├── S01_Conducta.m         Behavioural stats (IES, effect sizes)
│   ├── S02_TF_Correlaciones.m TF–behaviour correlations       [*]
│   ├── S03_TF_ClusterPerm.m   Cluster-permutation (5 000 perms)[*]
│   ├── S04_FOOOF_Figuras.m    FOOOF spectral decomposition
│   ├── S05_FOOOF_LME.m        LME: Group × Condition (N = 45)
│   ├── S06_PAC.m              Phase–amplitude coupling (PAC)   [*]
│   ├── S07_Supplementary.m    Steiger Z, mediation, ITPC, Rayleigh
│   ├── S08_ChewFreq.m         Chewing-frequency analysis       [*]
│   └── S09_Figuras_Final.m    Figure compositor
│
├── plots/                     Figure-specific plot scripts (read workspaces)
│   ├── P_PAC_Panels.m         Fig. 3 — PAC heatmap, scatter, violin
│   ├── P_Rayleigh_Panels.m    Supp. Fig. S4 — polar phase histograms
│   └── P_Supp_Panels.m        Supp. Fig. S3 — FOOOF, β-specificity, α-ERD
│
├── data/
│   ├── data_beh_tb_45.mat     Behavioural data (N = 45, IES, RT, ACC)
│   ├── incluidos45.mat        Participant list (30 Cases + 15 Controls)
│   └── workspaces/
│       └── FOOOF_Workspace.mat  Spectral decomposition (used by S04, S05)
│
└── outputs/                   Created at runtime; pre-computed .mat files
    ├── Figure02_TF/
    │   ├── chk_all_pwr_clean.mat    TF power matrix (N = 45)
    │   └── chk_all_itpc_clean.mat   ITPC matrix (N = 45)
    ├── Figure02b_TF_Correlaciones/
    │   └── TF_band_metrics.mat      Band-power × IES correlations
    └── Figure04_PAC/
        ├── PAC_4Groups_Workspace.mat  Main PAC results (used by S07, plots)
        └── PAC_EEG_Workspace.mat      EEG-phase PAC workspace
```

`[*]` Scripts marked require raw EEG data (see [Full pipeline](#full-pipeline) below).

---

## What each figure contains

| Paper figure | Script(s) | Output location |
|---|---|---|
| **Fig. 1A–B** — Behaviour (IES, RT) | `S01_Conducta` | `outputs/Figure01_Behavior/` |
| **Fig. 1C–D** — FOOOF PSD + exponent | `S04_FOOOF_Figuras` | `outputs/Figure03_FOOOF/` |
| **Fig. 2** — TF cluster-permutation | `S03_TF_ClusterPerm` [*] | `outputs/Figure02_TF/` |
| **Fig. 3A** — PAC heatmap (ΔzMI) | `P_PAC_Panels` | `outputs/figures/P_Fig3A_PAC_Heatmap.png` |
| **Fig. 3B** — β-PAC × IES scatter | `P_PAC_Panels` | `outputs/figures/P_Fig3B_BetaEarly_vs_IES.png` |
| **Fig. 3D** — Violin 2×2 Group×Cond | `P_PAC_Panels` | `outputs/figures/P_Fig3D_Violin_BetaMid.png` |
| **Supp. S2** — Supplementary stats | `S07_Supplementary` | `outputs/Supplementary/` |
| **Supp. S3** — FOOOF residual + specificity + α-ERD | `P_Supp_Panels` | `outputs/figures/FigSN*.png` |
| **Supp. S4** — Rayleigh polar plots | `P_Rayleigh_Panels` | `outputs/figures/P_Rayleigh_*.png` |

All figures requiring cross-panel compositing (Figs. 1, 3) are provided as individual panels ready for assembly in Inkscape or similar.

---

## Statistical outputs

The MATLAB console and `outputs/Figure03_FOOOF/Reporte_FOOOF_LME.txt` report:

- **Behavioural:** IES Wilcoxon signed-rank, Cohen's *d*, Mann–Whitney *U* (baseline equivalence)
- **FOOOF LME:** `Exponent ~ Group × Condition + (1|ID)`, *t*(86), *p*-value, *N* = 45
- **PAC:** Δ*z*MI Wilcoxon per band × window, Rayleigh *R* and *p*
- **Brain–behaviour:** Spearman ρ (β-PAC × IES; α-power × IES), Steiger *Z*, bootstrapped mediation (2 000 iterations)

---

## Full pipeline

Scripts S02, S03, S06, and S08 require the raw EEG dataset (~12 GB total). To run them:

1. Request the preprocessed EEG files from the corresponding author: `sebastian.espinoza@uv.cl`
2. In `S0_paths.m`, set `DATA_ROOT` to your local copy of the dataset.
3. Install [EEGLAB](https://sccn.ucsd.edu/eeglab/) and set `P.eeglab_path` accordingly.
4. Run scripts in order: S01 → S02 → S03 → S04 → S05 → S06 → S07 → S08 → S09  
   (or run `S00_RunAll.m` which executes the full sequence).

> **Note:** S06_PAC uses `parfor` over 200 surrogates; a parallel pool accelerates it ~10× (tested on Ryzen 5900X, 12 cores). The Parallel Computing Toolbox is optional — S06 degrades gracefully to serial execution.

---

## Experimental design

### Participants
- **N = 45**: 30 Cases (chewing intervention) + 15 Controls (no chewing)
- Healthy right-handed adults without temporomandibular disorders (DC/TMD screened)
- Two cases excluded for EEG/EMG quality: E3S3, E3S5

### Task
- Visuospatial 2-back working memory task (200-ms stimuli, 2000-ms ISI, 25% targets)
- Performance metric: Inverse Efficiency Score (IES = median RT / proportion correct)

### Conditions
- **Cases:** Block 1 = No-Chew (baseline) → Block 2 = Chew (specially formulated gum)
- **Controls:** Block 1 + Block 2 = No-Chew (practice-effect control)

### EEG/EMG recording
- BioSemi ActiveTwo: 64 scalp channels + 8 ExG channels (bilateral masseter EMG)
- Sampling rate: 1024 Hz → downsampled to 256 Hz
- Reference: common average (post-ICA)

### Preprocessing pipeline
1. 0.5 Hz HPF (zero-phase FIR)
2. ZapLine-Plus (50 Hz line noise)
3. iCanClean spatial regression (masseter EMG reference, 20 Hz HPF)
4. ASR (burst criterion *k* = 15, 3-min resting baseline)
5. AMICA ICA (trained on temporary 1.5 Hz HPF copy)
6. ICLabel classification + EMG-envelope correlation rejection
7. Spherical interpolation → common average re-reference

### Key parameters
| Parameter | Value |
|---|---|
| ROI | F1, FC1 (channels 12 & 16) |
| Frequency bands | θ: 4–7 Hz · α: 8–12 Hz · β: 13–30 Hz |
| Post-stimulus windows | Early: 0–300 ms · Mid: 200–700 ms · Late: 300–900 ms |
| PAC method | Tort MI (18 bins, 200 surrogates → *z*MI) |
| Cluster permutations | 5 000 · threshold *p* < 0.05 · cluster mass statistic |
| FOOOF mode | Fixed aperiodic (3–40 Hz, no knee) |

---

## Citation

```bibtex
@article{espinoza2025,
  title   = {Rhythmic mastication entrains frontocentral β-band amplitude via
             corticomuscular phase–amplitude coupling during working memory},
  author  = {Espinoza, Sebastian and Moraga-Espinoza, Daniel and
             Morreal-Ortega, Luis and El-Deredy, Wael},
  journal = {Scientific Reports},
  year    = {2025},
  doi     = {10.1038/s41598-025-27606-5}
}
```

---

## License

Code: MIT License.  
Data: CC BY 4.0 — please cite the paper above if you use these materials.

---

## Contact

Sebastian Espinoza · `sebastian.espinoza@uv.cl`  
Universidad de Valparaíso, Chile
