# TITLE + RESULTS RESTRUCTURE — Paper2 V1 (marco DUAL θ+β-CMC) · 2026-06-20

## TÍTULO — opciones (todas: SANO, sin "entrains via θ-PAC específico", honestas)
Título actual (θ-céntrico, a reemplazar):
> "Rhythmic mastication entrains frontocentral theta via phase–amplitude coupling during active working memory"

**Opción 1 (ELEGIDA — dual, sobria):**
> "Rhythmic mastication couples to frontocentral theta and broadband cortical dynamics during working memory"
- Capta lo dual (θ frontal-medial + acoplamiento broadband/CMC) sin afirmar θ-específico ni β=CMC.

**Opción 2 (mecanística, aperiódico al frente):**
> "Chewing couples cortical activity to a peripheral motor rhythm and flattens the aperiodic spectrum during working memory"
- Resalta el hallazgo aperiódico (el más robusto) + el acoplamiento a ritmo periférico; menos foco WM/θ.

**Opción 3 (CMC-forward):**
> "Rhythmic mastication drives broadband corticomuscular coupling and frontal-midline theta during visuospatial working memory"
- Más β-CMC; riesgo: "corticomuscular coupling" debe leerse como análogo a CMC (guardarraíl).

Aplicado a `main.tex`: Opción 1.

## RESULTADOS — outline reordenado (cada subsección anclada a su figura)
Orden por cadena de claims: conducta → manip-check → TF → FOOOF → PAC(θ+banda) → control negativo.

1. **Behaviour — chewing speeds performance in chewers** → `behavior_RT_by_block_group`,
   `behavior_IES_by_block_group`. RT 591→519 Wilcoxon p=0.0001 d=−0.58; IES 994→812 p=0.0008.
   Baseline equiv (p=0.98). Diferencial Grupo×Bloque = tendencia (p=0.054). (Sin correlación PAC×conducta: n.s.)
2. **Manipulation check** → `manipcheck_masseter_SNR_ChNc_group`. Masetero +6.1 dB 31/31 (p=1.2e-6) vs ctrl ≈0.
3. **Chewing increases frontal-medial theta (TF + topo)** → `tf_cases_…`, `tf_controls_…`, `tf_interaction`,
   `topo_theta_interaction`. Casos cluster p=0.005; ctrl p=1.0; interacción p=0.010; Fpz/AFz FDR; temporal 0/6.
   (En eje tiempo θ/α/β todas suben → broadband; especificidad θ = topográfica.)
4. **Genuine theta peak over a flatter 1/f** → `fooof_PSD_aperiodic_fit_groups`,
   `fooof_exponent_cases_vs_controls`, `fooof_theta_peak_overlay`. Pico θ 97% CF≈6.6; χ 1.110→0.735 p<0.0001.
5. **Cortical amplitude couples to the masticatory phase (θ + broadband)** → `pac_zMI_theta_doublenull`,
   `pac_rayleigh_theta_rose`, `pac_zMI_band_comparison_doublenull` (β strongest), `pac_rayleigh_beta_rose`,
   `comodulogram_descriptive_MI`. θ 22/31 circ + 21/31 AAFT; broadband β>α>θ Friedman p<0.0001; Rayleigh θ Z=13.47.
   Marco dual: stream β = acoplamiento corticomuscular (análogo CMC); stream θ = WM (especificidad por topo+cond+FOOOF).
6. **The coupling follows real chewing (condition + negative control)** → `pac_zMI_theta_chew_vs_nochew_negctrl`,
   (suppl `topoS_theta_controls`). 2×2: Ch 13.32 (22/31) ≫ Nc 1.08 (7/31) ≈ Controls 0.52 (1/15, p=0.094).
7. **Artifact controls** (→ suppl): `artifactS_spatial_frontal_vs_temporal`, `artifactS_masseterRMS_x_dtheta`
   (ρ=−0.22 n.s.), comodulogram 2D. Confound motor NO sostiene el efecto.

## NÚMEROS A PURGAR del 3_results.tex actual (STALE/contradichos)
- RT 617.9→543.7 / p=0.0321 → 591→519 Wilcoxon p=0.0001.
- θ×IES ρ=−0.3654 p=0.0216 → **ELIMINAR** (canónico n.s.).
- "dose–response" / "θ-band enhancement" / "recovery" → fuera.
- PAC 1.0 Hz × 6.5 Hz (post-hoc) → individual f_chew±0.5 Hz, θ 4–7.
- MI median 1.40e-4 + 6/31 sig + Rayleigh Z=5.097 R=0.405 → 22/31 doble null + Rayleigh θ Z=13.47 R=0.659.
- Heading concatenado "Frontocentral θ reflects Changes…" → separar (Minor 3).
- θ peak × chewing freq ρ=0.324 p=0.038 → conservar (descriptivo, sin "enhancement").

## INTRO — fix crítico (falsación de artefacto)
El intro actual dice: "si EMG conduce, PAC pico en alta frecuencia + lateral/temporal; PAC θ-específico
frontocentral = origen cortical". **Contradicho:** nuestro PAC es broadband β-dominante. Reescribir la lógica
de falsación: el origen no-muscular se sostiene por (a) **topografía** frontal-medial vs temporal 0/6,
(b) supervivencia a **AAFT** (rompe forma de onda), (c) **especificidad de condición** + control negativo,
NO por especificidad de frecuencia del PAC. Y θ 4–8 → 4–7.
