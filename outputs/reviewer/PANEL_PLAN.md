# PANEL_PLAN — Paper2 V1 figuras (Fases A+B) · marco DUAL θ-cognitivo + β-CMC
Fecha: 2026-06-20 · Para aprobación ANTES de generar. Cada panel = ARCHIVO INDIVIDUAL, cuadrado, nombre explicativo.

## ESTÁNDAR ÚNICO (retrofit obligatorio — los scripts viejos NO lo cumplen)
- **Cuadradas 1:1**, 3.5×3.5 in (88.9 mm). Export **SOLO PNG @300dpi** (confirmado usuario 2026-06-20: nada de PDF/SVG). (Los scripts viejos usan 4×5, 2.8×4, 3.2×3.2 → corregir.)
- **Tipografía homogénea Arial:** título 12 pt · label de eje 11 pt · ticks 9 pt · leyenda 9 pt. Idéntico en TODAS.
- **Paleta Okabe-Ito (S0_config):** CASE #009E73 · CTRL #0072B2 · NC #888 · θ #D55E00 · α #56B4E9 · β #CC79A7 · SIG rojo.
- **Etiquetas SIN TMD:** "Chewing group (Cases, n=31)" · "Control / no-chew (n=15)". Nunca clínico/paciente/TMD.
- Salida → `outputs/figures/` (NO `figures_ok/`). Validar cada una con `figure-qa`.
- Estilo a imitar: `outputs/figures_ok/`, `..\Paper_plots\Final_Figures\`.

## ⚠ CAMBIO POR MARCO DUAL (proponer): el set legacy es θ-céntrico
El marco dual eleva β-CMC a co-titular. Por eso **promuevo a PRINCIPAL** lo que era suplementario y **agrego β**:
- `pac_zMI_band_comparison_doublenull` (θ/α/β, β el más fuerte) → de Suppl a **MAIN** (es la evidencia del stream B).
- `pac_rayleigh_beta_rose` (β Z=13.16) → **MAIN nuevo**, junto al theta rose, para el acoplamiento β.
- Comodulograma descriptivo: rotular ejes mostrando que el cluster es broadband fa=4–30 Hz (sostiene β + θ).

---

## FASE A — PANELES PRINCIPALES

| # | Archivo (.png/.svg/.pdf) | Qué muestra | Fuente .mat → variables | Script base | Motor |
|---|---|---|---|---|---|
| A1 | `behavior_RT_by_block_group` | RT B1/B2 pareado × grupo; Casos 591→519 (Wilcoxon p=0.0001, d=−0.58) | data/computed/v1_S1_behavior_stats.mat (+data_beh_tb.mat) | S1_behavior_figs.py | py |
| A2 | `behavior_IES_by_block_group` | IES B1/B2 × grupo; Casos 994→812 (p=0.0008, d=−0.49) | idem | S1_behavior_figs.py | py |
| A3 | `manipcheck_masseter_SNR_ChNc_group` | EMG masetero Ch−Nc; Casos +6.1 dB (31/31, p=1.2e-6) vs Ctrl ≈0 | chew_engagement_check.mat → snrCh,snrNc,casd,ctrd,p_cas,p_ctr,p_bt | nuevo | py |
| A4 | `tf_cases_chew_minus_nochew` | TF Chew−NoChew casos + contornos CBPT | S2c_TF_GroupFigure.mat → map_cas,mask_cas,times_anal,FREX_TF | S2c_export_panels.m | MATLAB |
| A5 | `tf_controls_chew_minus_nochew` | TF controles (n.s.) | → map_ctr,mask_ctr | S2c_export_panels.m | MATLAB |
| A6 | `tf_interaction` | TF interacción Grupo×Cond + contornos | → map_int,mask_int | S2c_export_panels.m | MATLAB |
| A7 | `topo_theta_interaction` | Topo θ interacción ROI-18; Fpz/AFz FDR marcados | S2c → topo_int,t_topo_int,sig_int_fdr **o** ROI_canonical.mat → theta_topo_int_M,sig_int_M,ROI_IDX | S2c/S2d | MATLAB+EEGLAB |
| A8 | `fooof_PSD_aperiodic_fit_groups` | PSD + ajuste aperiódico, casos vs controles | FOOOF_Workspace_V1.mat → GR(struct),ROI_FOOOF,FIT_RANGE | S3c_FOOOF_paper_figs.m | MATLAB |
| A9 | `fooof_exponent_cases_vs_controls` | Exponente χ 1.110→0.735 casos (p<0.0001); ctrl n.s. | → GR exponents | S3c_FOOOF_paper_figs.m | MATLAB |
| A10 | `fooof_theta_peak_overlay` | Pico θ genuino (97%, CF≈6.58 Hz) sobre aperiódico | → GR periodic | S3c_FOOOF_paper_figs.m | MATLAB |
| A11 | `pac_zMI_theta_doublenull` | zMI θ vs circ (22/31) + AAFT (21/31) | v1_S4b_PAC_ROI.mat → zC_ch,zA_ch[θ],BAND_NAMES,BANDS_HZ | S4c_PAC_figures.py | py |
| A12 | `pac_zMI_theta_chew_vs_nochew_negctrl` | 2×2 mecanismo: Casos-Ch 13.3 (22/31) ≫ Casos-Nc 1.08 (7/31) ≈ Ctrl 0.52 (1/15) | v1_S4b → zC_ch,zC_nc[θ] + S8_controls_PAC.mat → zC,k,nv | S4c + nuevo | py |
| A13 | `pac_rayleigh_theta_rose` | Rosa fase preferida θ; Z=13.47 | v1_S4b → pref_ch,rayl_R,rayl_Z,rayl_p[θ] | S4d_Rayleigh_rose.py | py |
| **A14** | `pac_zMI_band_comparison_doublenull` ⬆ | θ/α/β zMI doble null; **β el más fuerte** (88–120); Friedman p<0.0001 | v1_S4b → zC_ch,zA_ch[θ,α,β] | S4c (extender) | py |
| **A15** | `pac_rayleigh_beta_rose` ⬆ nuevo | Rosa fase β; Z=13.16 (acoplamiento β = stream CMC) | v1_S4b → pref_ch,rayl_Z[β] | S4d (extender) | py |
| A16 | `comodulogram_descriptive_MI` | Comodulograma MI (rotular **"descriptivo"**, broadband fa=4–30) | S4g_comodulogram.mat → MI_mean,F_AMP,F_PHASE | S4g_comodulogram.m | MATLAB |

## FASE B — PANELES SUPLEMENTARIOS (todos individuales, mismo estándar)

| # | Archivo | Qué muestra | Fuente .mat → variables | Motor |
|---|---|---|---|---|
| B1 | `pacS_band_specificity_zMI` | (si A14 sube a main, B1 queda como versión extendida/por-sujeto) | v1_S4b | py |
| B2 | `artifactS_spatial_frontal_vs_temporal` | A1: frontal sig vs temporal 0/6 | S6_artifact_controls.mat → fm_idx,tp_idx,n_fm_sig,n_tp_sig,p_el,p_fdr,labels | py |
| B3 | `artifactS_masseterRMS_x_dtheta` | A3 EMG masetero real × Δθ ρ=−0.22 (n.s.) | S6b_A3_muscle.mat → emg_rms,dtheta,r_rms,p_rms | py |
| B4 | `artifactS_sharpness_x_zMI` | A4 sharpness × zMI ρ=0.17 (n.s.) | S6_artifact_controls.mat → rho_A4,p_A4 | py |
| B5 | `artifactS_fooof_dissociation` | A2 disociación χ + pico θ vs banda alta | FOOOF_Workspace_V1.mat + S6 | MATLAB |
| B6 | `doseS_engagement_x_effects` | S7 intensidad × efectos (n.s.) | S7_dose_response.mat → eng,dth,drt,dies,miL,zmi,R,pairs | py |
| B7 | `pacS_MI_by_window` | MI θ Base/Early/Late (Late max; Friedman p=0.0007) | v1_S4b → MIw_ch (ventanas) | py |
| B8a-c | `topoS_theta_{cases,controls,interaction}` | Topo θ × 3 contrastes (paneles separados) | S2c topo_* | MATLAB+EEGLAB |
| B9a-c | `topoS_alpha_{...}` | Topo α × 3 — **⚠ requiere fuente topo por-banda** (S2c parece θ-céntrico; verificar/recomputar) | (pendiente) | MATLAB+EEGLAB |
| B10a-c | `topoS_theta_alpha_{...}` | Topo θ−α × 3 | (pendiente) | MATLAB+EEGLAB |
| B11 | `comodulogramS_2D_cluster` | Cluster 2D significativo (broadband) | S4h_comodulogram_stats.mat → sig_mask,clust_label,p_clusters,F_AMP,F_PHASE | MATLAB |
| B12 | `comodulogramS_2D_zmap` | Mapa z 2D | S4h → z_mean,z_maps_out,t_map,p_map | MATLAB |

## RIESGOS / PENDIENTES TÉCNICOS (resolver antes o durante generación)
1. **Topoplots (A7, B8–B10)** requieren EEGLAB + chanlocs en MATLAB; usar `eeglab nogui`.
2. **FOOOF GR struct** está anidado en v7.3 → inspeccionar `GR` antes de A8–A10.
3. **topoS por banda (B9/B10)**: S2c_TF_GroupFigure.mat parece θ-céntrico (topo_cas/ctr/int). Verificar si hay
   topo α/β; si no, recomputar desde v1_S2_TF_data.mat (S2c método) — coste medio.
4. **Behaviour**: S1_behavior_figs.py separa Cases/Controls en panels distintos; A1/A2 los fusionan en un
   panel agrupado pareado por grupo → reescritura ligera.

## ORDEN DE EJECUCIÓN PROPUESTO
py primero (rápido, sin EEGLAB): A1–A3, A11–A15, B2–B4, B6–B7 → luego MATLAB (TF/topo/FOOOF/comodulo):
A4–A10, A16, B5, B8–B12. Validar cada una con `figure-qa`. Nada a `figures_ok/`.
