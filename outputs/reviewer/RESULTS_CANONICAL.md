# RESULTS_CANONICAL — Paper2 V1 (corrida a prueba de balas)
Fecha: 2026-06-19 · Fuente única: S0_config.m · θ=4–7 (guarda 7–8) α=8–13 β=13–30 · ROI-18 · N=31+15 · seed=42

> Reemplaza los números STALE de PAC_results_summary.txt y FOOOF_results_summary.txt.
> Cada número recomputado esta corrida desde datos canónicos. Sidecars de provenance junto a cada .mat.

## DISEÑO (sin TMD — ambos grupos sanos)
Manipulación de masticación entre sujetos: "Casos" (n=31) masticaron en el 2º bloque; "Controles" (n=15)
repitieron la tarea sin masticar (control de práctica). Ambos: B1=NoChew, orden fijo. La interacción
Grupo×Bloque aísla el efecto de masticación controlando práctica (I6).
**Manipulation check** (masetero 1–2.5 Hz, S0c_chew_engagement.m → chew_engagement_check.mat):
Casos Ch−Nc=+6.10 dB (31/31, Wilcoxon p=1.2×10⁻⁶) · Controles −0.72 dB (5/15, n.s.) · entre grupos p=7×10⁻⁸.
NO usar lenguaje clínico/TMD.

---

## ROI canónico (fix I4)
- **ROI-18** = interacción Cases×Cond FDR<0.05 sobre `theta_topo` (Morlet NPL dB, WIN_LATE) en
  v1_S2_TF_data.mat. Reproducido exactamente por `code/S2d_theta_ROI.m` → `ROI_canonical.mat/.txt`.
  Electrodos: Fp1 AF7 AF3 F1 F3 F5 F7 FC1 Fpz Fp2 AF8 AF4 AFz Fz F2 F4 FC2 FCz.
- Robustez: método de épocas (eegfilt+Hilbert²) → 13 electrodos, **todos subconjunto** de los 18.
- Casos Ch>Nc (theta, FDR): **Fpz, AFz** (responde R1.M6).

## Decisión de banda (2026-06-19)
- θ=4–7 con **7–8 Hz como guarda declarada**. Se evaluó θ=4–8 (convención WM) y S2d mostró que
  diluye la interacción midline (pierde Fz/FCz → ROI frontopolar de 11). 4–7 preserva el ROI
  cognitivo midline → mejor defensa anti-artefacto. (Argumento de especificidad, R1.M9.)

---

## Conducta (S1) — v1_S1_behavior_stats.mat (N=31 casos, 15 ctrl) · ap_=N=30 casi idéntico
| | NoChew B1 | Chew B2 | Δ casos | test (casos, Wilcoxon) |
|---|---|---|---|---|
| RT casos | 591±112 ms | 519±123 ms | **−72.4** | **p=0.0001** (d=−0.58) |
| RT controles | 584±98 | 571±176 | −13.5 | — |
| IES casos | 994±496 | 812±437 | **−182** | **p=0.0008** (d=−0.49) |
| IES controles | 939±251 | 984±679 | +45 | — |
| Acc casos | 0.736 | 0.787 | +0.05 | — |
- LME RT: Bloque p=6.3×10⁻⁵ · Grupo p=0.86 (sin dif. basal) · **Grupo×Bloque p=0.054 (tendencia)**.
- LME IES: Grupo×Bloque p=0.106 (n.s.).
- **DOS COSAS DISTINTAS:** (1) la mejora conductual con masticación dentro de casos ES significativa
  (Wilcoxon RT/IES, d≈−0.5); (2) lo n.s. es la correlación PAC×conducta (A6). No confundir.
- **STALE corregido:** el resumen viejo usaba t-pareado (p_rt_delta=0.14, p_ies_delta=0.66) — test
  equivocado para RT/IES sesgados; el Wilcoxon (canónico) da p=0.0001 / 0.0008.
- **Caveats (R5):** el beneficio diferencial casos>controles es solo tendencia (Grupo×Bloque p=0.054);
  confundido con orden fijo NoChew→Chew (práctica, I6).

## PAC principal — `S4b_PAC_ROI.m` → `v1_S4b_PAC_ROI.mat`
EMG = promedio bilateral 65+66 (fallback canal vivo). Fase = f_chew_ind±0.5 Hz. Doble null:
circular-shift (min 5 s) + AAFT (N_SURR=500 c/u). zMI vs null de SU propia condición.
f_chew recomputado: M=1.547 ± 0.152 Hz (30/31; 1 sujeto vía fallback EMG-PSD).

### P1 — ¿PAC theta > null? (Cases-Chew, N=31)
| Null | M z | Mdn z | k>1.96 | binomial | Wilcoxon (1-sided) |
|---|---|---|---|---|---|
| circular-shift | 13.32 | 9.16 | **22/31** | ≈0 | 7.82×10⁻⁶ |
| AAFT (no-sinusoidalidad) | 12.85 | 8.20 | 21/31 | ≈0 | 1.03×10⁻⁵ |
| **supera AMBOS** | — | — | **21/31** | ≈0 | — |
**Veredicto:** PAC theta genuino; sobrevive AAFT (rompe forma de onda) → no es artefacto de
transiente no-sinusoidal. **22/31 reemplaza el STALE 15/31.**

### P2 — ¿Específico de banda? (zMI por banda, doble null)
| Banda | M z circ | k circ | p circ | M z AAFT | k AAFT |
|---|---|---|---|---|---|
| theta | 13.32 | 22/31 | 7.8×10⁻⁶ | 12.85 | 21/31 |
| alpha | 22.55 | 24/31 | 3.3×10⁻⁶ | 22.27 | 24/31 |
| beta | 88.10 | 29/31 | 1.0×10⁻⁶ | 120.44 | 29/31 |
Friedman zMI×banda: circ p<0.0001 | aaft p<0.0001.
**Veredicto (honesto, I2/C2):** el PAC es **broadband**, beta el más fuerte. La especificidad theta
**NO** se sostiene en la magnitud del PAC; se sostiene en (a) la topografía de interacción theta
(frontal-medial, ROI-18), (b) la especificidad de condición (P4), (c) el pico theta FOOOF genuino.

### P3 — ¿Específico de ventana? (MI Tort por ventana, theta)
| Banda | MI Base | MI Early | MI Late | B→E | B→L | E→L |
|---|---|---|---|---|---|---|
| theta | 5.3e-4 | 3.4e-4 | 1.03e-3 | 0.008 | 0.086 | 0.001 |
| alpha | 1.07e-3 | 5.1e-4 | 1.46e-3 | 0.004 | 0.024 | 0.000 |
| beta | 1.34e-3 | 1.00e-3 | 1.29e-3 | 0.001 | 0.491 | 0.003 |
Friedman theta×ventana p=0.0007.
**Veredicto:** theta MI máximo en Late (Late>Early p=0.001); Late vs Base no sig (p=0.086).

### P4 — ¿Difiere Ch vs Nc? (null POR condición — fix I3)
| Medida | Ch | Nc | Wilcoxon |
|---|---|---|---|
| zMI θ circ | M=13.32 | M=1.08 (7/31 sig) | **p=0.0001** |
| zMI θ AAFT | M=12.85 | M=1.13 | p=0.0001 |
| MI θ abs | — | — | p=0.0001 |
**Veredicto:** el PAC theta está impulsado por la masticación; casi ausente en NoChew. Pilar mecanístico.

### Rayleigh fase preferida (Cases-Chew)
θ: R=0.659 Z=13.47 p<0.0001 · α: R=0.605 Z=11.33 p<0.0001 · β: R=0.652 Z=13.16 p<0.0001.
(θ reproduce el histórico R=0.653/Z=13.24.)

### PAC × conducta (Spearman, Cases-Chew, N=31)
zMI θ circ × RT ρ=+0.01 p=0.97 · MI θ Late × RT ρ=−0.22 p=0.23 · × IES ρ=−0.08 p=0.69 (todos n.s.).
**El ρ=−0.377 (MI×medRT) histórico NO se reproduce canónicamente → CONTRADICHO.**

---

## Comodulograma F6 — `S4g`/`S4h`
- S4g descriptivo (MI Tort) + S4h stats 2D (N_SURR=500, min-shift 12.5 s, sign-flip CBPT 5000).
- Cluster 2D significativo p<0.0001 pero **broadband**: fa=[4–30 Hz] (θ+α+β), fp=[0.4–1.8 Hz].
- Reportado como **descriptivo** (C2/A4); confirma que el acoplamiento no es θ-exclusivo.

## Controles de artefacto A1–A6 → ver ARTIFACT_CONTROLS.md
Veredicto: el confound motor NO sostiene el efecto. A1 frontal-medial (Fpz pFDR=0.013, AFz 0.038;
temporal 0/6) · A2 χ 1.11→0.74 + pico θ 97% · A3 **EMG masetero real × Δθ ρ=−0.22 n.s.** ·
A4 AAFT 21/31 + sharpness ρ=0.17 n.s. · A5 Ch≫Nc p=0.0001 · A6 conducta n.s.

## [PENDIENTE] S2a re-verificar · FOOOF re-verificar (ya en 4–7) · F8 reconciliación · carta R1
