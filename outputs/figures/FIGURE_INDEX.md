# FIGURE_INDEX — Paper2 V1 paneles (marco DUAL θ+β-CMC) · 2026-06-20
Estándar: cuadradas 3.5×3.5 in · **PNG @300 dpi** · Arial 12/11/9 · Okabe-Ito · etiquetas SIN TMD.
QA programático: 24/24 cuadradas (aspect 0.90–1.00 roses/topos incluidas), 300 dpi, RGB(A). Scripts:
`code/figA_py.py` · `code/figA2_py.py` · `code/figA_topo.m` · `code/figB_py.py` · estilo `code/_figstyle.py`.

## PRINCIPALES (16)
| Panel | Fuente .mat | Claim |
|---|---|---|
| behavior_RT_by_block_group | v1_S1_behavior_stats | RT 591→519 casos, Wilcoxon p=0.0001 d=−0.58 |
| behavior_IES_by_block_group | v1_S1_behavior_stats | IES 994→812, p=0.0008 d=−0.49 |
| manipcheck_masseter_SNR_ChNc_group | chew_engagement_check | masetero +6.1 dB 31/31 (p=1.2e-6) vs ctrl n.s. |
| tf_cases_chew_minus_nochew | S2c_TF_GroupFigure (map_cas,mask_cas) | cluster CBPT casos p=0.005 |
| tf_controls_chew_minus_nochew | S2c_TF_GroupFigure (map_ctr) | controles n.s. p=1.0 |
| tf_interaction | S2c_TF_GroupFigure (map_int,mask_int) | interacción θ-α p=0.010 |
| topo_theta_interaction | v1_S2_TF_data (theta_topo_*) | frontal-medial, Fpz/AFz FDR |
| fooof_PSD_aperiodic_fit_groups | FOOOF_Workspace_V1 (GR) | PSD + ajuste aperiódico |
| fooof_exponent_cases_vs_controls | FOOOF_Workspace_V1 (GR.exp) | χ 1.11→0.74 casos p<0.0001; ctrl n.s. |
| fooof_theta_peak_overlay | FOOOF_Workspace_V1 (GR.Res) | pico θ genuino 97%, CF≈6.6 Hz |
| pac_zMI_theta_doublenull | v1_S4b_PAC_ROI (zC_ch,zA_ch) | θ-PAC 22/31 circ + 21/31 AAFT |
| pac_zMI_theta_chew_vs_nochew_negctrl | v1_S4b + S8_controls_PAC | 2×2: 22/31 ≫ 7/31 ≈ 1/15 |
| pac_rayleigh_theta_rose | v1_S4b_PAC_ROI (pref_ch,rayl) | θ Z=13.47, R=0.66 |
| **pac_zMI_band_comparison_doublenull** ⬆ | v1_S4b_PAC_ROI (zC_ch) | broadband, β el más fuerte (29/31), Friedman p<0.0001 |
| **pac_rayleigh_beta_rose** ⬆ | v1_S4b_PAC_ROI (pref_ch[β]) | β Z=13.16 (stream CMC) |
| comodulogram_descriptive_MI | S4g_comodulogram | descriptivo, broadband fa 4–30 Hz |

⬆ = promovido a principal (evidencia de banda; el β se reinterpretó luego como leakage EMG).

### Figura 6 — puente conductual↔EEG single-trial (fig:bridge en 6_figs.tex)
| Panel | Fuente | Claim |
|---|---|---|
| bridge_thetaLate_RT_quintiles | rev_singletrial.csv | RT↓ por quintil de θ-late within-subj; LMM slope −41.8 ms p=0.0007, 23/31 |
| bridge_window_specificity | rev_singletrial.csv | link específico de ventana tardía (pre/early n.s.); t-test Fisher-z p=0.0016 |
Script: `code/rev_singletrial.py` (+ `rev_singletrial_fig.py`). Va en RESULTADOS (data del estudio actual).

## SUPLEMENTARIOS (8 generados)
| Panel | Fuente | Claim |
|---|---|---|
| topoS_theta_cases | v1_S2_TF_data | θ Ch−Nc casos |
| topoS_theta_controls | v1_S2_TF_data | θ Ch−Nc controles (≈0) |
| artifactS_spatial_frontal_vs_temporal | S6_artifact_controls | A1 frontal sig vs temporal 0/6 |
| artifactS_masseterRMS_x_dtheta | S6b_A3_muscle | A3 EMG×Δθ ρ=−0.22 n.s. (Spearman) |
| doseS_engagement_x_dtheta | S7_dose_response | dosis-respuesta n.s. |
| pacS_MI_by_window | v1_S4b_PAC_ROI (MIw_ch) | θ MI Base/Early/Late, Late max (Friedman p=0.0007) |
| comodulogramS_2D_cluster | S4h_comodulogram_stats | cluster 2D sig (broadband) |
| comodulogramS_2D_zmap | S4h_comodulogram_stats | mapa z 2D |
| fooofS_aperiodic_offline_paper1 | rev_paper1_aperiodic.npz (D:\Exp1) | **re-análisis cross-study**: aperiódico offline (Paper 1) χ 1.45→1.35, Δ=−0.09, 21/30, p=0.005 |
| fooofS_aperiodic_online_paper2 | FOOOF_Workspace_V1 (GR) | aperiódico online (Paper 2) χ 1.11→0.74, Δ=−0.38, p<0.0001 (par para la fig cross-study) |

**Cross-study (Suppl Fig `fig:aperiodic_crossstudy` en 8_suppl.tex):** los dos paneles `fooofS_aperiodic_*`
forman la figura de convergencia aperiódica offline↔online. Etiquetada como re-análisis original (no resultado
de Espinoza 2025). Script: `code/fig_suppl_offline_aperiodic.py`. Narrativa: `outputs/reviewer/NARRATIVE_online_vs_offline.md`.

## DIFERIDOS (documentado por qué; baja prioridad / requieren fuente extra)
- `artifactS_sharpness_x_zMI` (B4): el vector sharpness×zMI no está en los .mat guardados (sólo rho_A4/p_A4 escalares). Requiere recomputar A4.
- `artifactS_fooof_dissociation` (B5): redundante con fooof_exponent + theta_peak; sin fuente limpia de banda 30–40 Hz por sujeto.
- `topoS_alpha_*`, `topoS_theta_alpha_*` (B9/B10): topos α/θ-α NO precomputados por canal (sólo θ). Requieren recomputo per-canal multibanda desde epochs (.set) — costo alto, valor suplementario.

## NOTA
Nada depositado en `figures_ok/` (carpeta de validación manual del usuario). Para asomar a Fpz/AFz se usó
emarker2 sobre los electrodos FDR. Paleta y tamaños idénticos en todas para ensamblaje en LaTeX por el usuario.
