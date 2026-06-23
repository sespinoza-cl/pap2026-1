# pap2026-1 — Code and Data for:

**"Rhythmic mastication couples to frontocentral theta and flattens the cortical aperiodic spectrum during working memory"**

Espinoza S., Moraga-Espinoza D., Morreal-Ortega L., El-Deredy W. (2026). *Manuscript under revision (R1).*
Repository: [github.com/sespinoza-cl/pap2026-1](https://github.com/sespinoza-cl/pap2026-1)

> **This branch (`v1-final`) supersedes the earlier β-CMC framing.** Re-analysis with a per-electrode
> corticomuscular-coherence control showed the broadband β dominance is **residual masseter EMG (diffuse,
> temporal topography), not corticomuscular coherence**. The cortically interpretable coupling is
> **frontal-medial theta** (frontal, AAFT-robust, condition-specific, genuine FOOOF peak), accompanied by an
> **aperiodic (1/f) flattening** indexing increased cortical excitability.

---

## Quick start

> **Requirements:** MATLAB R2025b + Statistics & Signal Processing Toolboxes + EEGLAB.
> Python is only required for FOOOF/specparam fitting (`fooof_fit.py`, called automatically by `S3b_FOOOF_fromPython.m`).

All pre-computed result workspaces are in `outputs/stats/`. To regenerate figures from them:

```matlab
% From Analysis_V1_Final/ in MATLAB:
run S0_config
run code/S2c_TF_GroupFigure   % TF + topo
run code/S3c_FOOOF_paper_figs  % FOOOF panels
run code/S4b_PAC_figs          % PAC panels
run code/rev_cmc_topography    % CMC topography
run code/figA_topo             % topoplots
```

Full step-by-step map (every result → script → data): see **[`REPRODUCE.md`](REPRODUCE.md)**.

---

## What's here
```
pap2026-1/
├── REPRODUCE.md              ← full reproducibility map (start here)
├── S0_config.m               ← single config (bands θ4–7/α8–13/β13–30, ROI-18, windows, seed=42)
├── code/                     ← analysis + figure scripts (MATLAB S*.m + fooof_fit.py)
│   ├── S1..S8 (canonical pipeline) + fooof_fit.py (FOOOF bridge for S3b)
│   └── rev_*  (revision analyses: single-trial bridge, CMC topography, cross-study aperiodic)
├── outputs/
│   ├── stats/                ← pre-computed .mat (reproduce figures without raw EEG)
│   ├── figures/              ← 300-dpi PNG panels (square, Okabe-Ito)
│   └── reviewer/             ← canonical results, claims, revision analyses, panel review (docs)
└── (raw EEG/EMG .set NOT included — see Data availability)
```

## Design
46 healthy adults: 31 Cases (chewed in block 2) + 15 Controls (no-chew practice control). Visuospatial 2-back.
BioSemi 64+8 (bilateral masseter EMG), 1024→256 Hz. ROI-18 frontal-medial. Bands θ 4–7 / α 8–13 / β 13–30 Hz.

## Data availability
Pre-computed workspaces are included. Raw/preprocessed EEG/EMG (`.set`, large) available from the corresponding
author: `sebastian.espinoza@uv.cl`. The analysis is driven by a single config (`S0_config.m`) with fixed seeds,
so every reported number regenerates deterministically.

## Citation
```bibtex
@article{espinoza2026,
  title  = {Rhythmic mastication couples to frontocentral theta and flattens the cortical
            aperiodic spectrum during working memory},
  author = {Espinoza, Sebastian and Moraga-Espinoza, Daniel and Morreal-Ortega, Luis and El-Deredy, Wael},
  year   = {2026}, note = {Manuscript under revision}
}
```

## License
Code: MIT. Data: CC BY 4.0. Contact: Sebastian Espinoza · `sebastian.espinoza@uv.cl` · Universidad de Valparaíso.
