# REVISION_ANALYSES — Respuesta empírica a los 2 bloqueantes del panel · 2026-06-20
Scripts: `code/rev_dissociation.py`, `code/rev_dissociation2.py` (#1) · `code/rev_beta_cmc.m` → `outputs/stats/rev_beta_cmc.mat` (#2).

================================================================================
## #1 — DISOCIACIÓN θ/aperiódico ↔ conducta (panel DA-1, D2, E1)
================================================================================
**Pregunta:** ¿algún cambio neural predice la mejora conductual a nivel individual (intra-casos, n=31)?
Incluye el "enfoque elegante v1" (frecuencia-pico por banda + amplitud media/mediana del ROI).

**Cambio-cambio (Δneural × Δconducta), Spearman:**
- Δθ-late × ΔRT ρ=+0.089 p=0.635 · Δθ-late × ΔIES ρ=+0.234 p=0.205 — n.s.
- Δexponente × ΔRT ρ=+0.172 p=0.355 · × ΔIES ρ=+0.136 p=0.466 — n.s.
- zMI θ (PAC) × ΔRT/ΔIES — n.s. (reproduce canónico)

**Enfoque v1 (pico-frecuencia + mean/median amplitud ROI por banda × conducta), 36 tests:**
- **0/36 nominalmente significativos** (esperados por azar ~1.8; Bonferroni α=0.0014: 0 sobreviven).
- Más cercanos: α-peakHz(chew)×ΔIES p=0.052; β-Δpeak×ΔIES p=0.077 — ambos n.s.

**Único hit en todo el barrido (~44 tests):** θ-late-power(topo, chew) × ΔIES ρ=+0.419 p=0.019 — NO sobrevive
corrección, NO se replica con mean/median del espectro, y su signo es **competitivo** (más θ → menos mejora),
no facilitador. Además Δexponente × Δθ-late ρ=+0.338 p=0.063 → θ sube CON el aplanamiento broadband.

**VEREDICTO #1:** los efectos neural y conductual están **desacoplados a nivel individual**, de forma robusta
a la métrica (potencia, exponente, PAC, pico-frecuencia, amplitud media/mediana). **No hay ancla cognitiva
rescatable.** → La revisión debe **reposicionar** la contribución como neurofisiológica/metodológica y declarar
el vínculo cognitivo como abierto (no insinuarlo). El alza de θ es compatible con excitabilidad/arousal broadband,
no con una firma θ-cognitiva selectiva.

================================================================================
## #2 — COHERENCIA CORTICOMUSCULAR β REAL (EEG-EMG) (panel D1, DA-2)
================================================================================
**Análisis:** magnitude-squared coherence (mscohere) entre EEG ROI frontocentral y EMG masetero (ch 65/66),
Cases-Chew, n=31. CMC clásica = coherencia fase-fase, canónica en β.

**Resultados:**
- **Magnitud CMC muy baja:** media β=0.001, mediana 0.001 (CMC clásica típica 0.05–0.30). θ=0.003 ≥ α=0.002 ≥
  β=0.001 → **NO es β-dominante**. El "10/31 sig>CL95" es engañoso: con datos largos el umbral analítico CL95
  cae al piso de ruido (0.001) y sobre-detecta; en magnitud absoluta NO hay CMC clásica fuerte.
- **PERO correlación banda-específica:** **CMCβ × zMI_β(PAC) ρ=+0.545, p=0.002** (significativo);
  CMCβ × zMI_θ(PAC) ρ=+0.232, p=0.209 (n.s.). Es decir, la variación inter-individual del β-PAC SÍ rastrea
  la (poca) coherencia corticomuscular β, y de forma específica de banda.

**CAVEAT METODOLÓGICO CLAVE:** el EEG de estos `.set` fue limpiado agresivamente (ASR + ICA + iCanClean) para
**remover componentes EMG-correlacionados** → eso ATENÚA por diseño la CMC genuina. Además el ROI frontal-medial
no es el sitio óptimo de CMC masetero (más central/lateral). Por tanto la magnitud baja es un **límite inferior
conservador**, no prueba de ausencia. Una CMC justa requiere EEG mínimamente limpiado + ROI central.

**VEREDICTO #2 (provisional):** NO se puede afirmar CMC clásica fuerte (magnitud al piso, no β-dominante).
El β-PAC covaría con CMC β (ρ=0.55, banda-específico).

### #2b — TOPOGRAFÍA de la CMC β (control definitivo: genuina vs leakage) — `rev_cmc_topography.m`
Prueba clave: la CMC corticomuscular GENUINA pica en sitios CENTRALES/motores; el leakage EMG puro pica en
TEMPORALES (cerca del masetero). Resultado (n=31, Cases-Chew):
- central (10 el) media=0.0011 · temporal (10 el) media=0.0015 · **Wilcoxon central vs temporal p=0.050
  (temporal ≥ central)** · **canal pico = T8 (temporal, 0.0018)**.
- Topoplot: distribución DIFUSA al piso de ruido, **sin foco central/motor**, con leve peso temporal/posterior.
  → `outputs/figures/artifactS_betaCMC_topography.png`.

**VEREDICTO #2 FINAL (importante, reorienta el marco):** NO hay evidencia topográfica de CMC genuina
(central). La coherencia β es difusa y, si acaso, **temporal-leaning** → patrón de **leakage EMG residual**,
no de drive corticomuscular. Por tanto la correlación β-PAC↔β-CMC (ρ=0.55) refleja más plausiblemente
**leakage compartido**, no CMC compartida. → **El "stream β-CMC" del marco dual NO está soportado.** La
supervivencia a AAFT del β-PAC sólo descarta artefacto de forma de onda, NO la conducción de volumen del EMG.

**IMPLICACIÓN DE MARCO (decisión del usuario):** el dato empuja DE VUELTA a una historia **θ-cortical**
(frontal-medial, temporal 0/6, AAFT-robusta, condición-específica, pico FOOOF = genuina), con la **dominancia
β/broadband atribuida explícita y honestamente a leakage EMG residual** (topografía difusa/temporal, sin pico
central), NO a CMC. Esto es más defendible y converge con la crítica DA-2 del panel. Cambios de texto que
implicaría: abstract y Discusión (subsección "corticomuscular beta coupling") pasan de "análogo a CMC" a
"la dominancia β es consistente con leakage EMG residual; el acoplamiento corticalmente interpretable es θ".

================================================================================
## IMPACTO EN LA DECISIÓN DEL PANEL
================================================================================
- **DA-1 (CRITICAL):** confirmado empíricamente → **reposicionar** (no des-clamar a medias). El paper es, con
  honestidad, una contribución **metodológica + neurofisiológica**: masticación recupera EEG válido y produce
  una huella cortical broadband (aperiódico↓, θ frontal-medial↑) + acoplamiento a la fase masticatoria, **sin
  consecuencia cognitiva individual demostrable**. Esa es la historia defendible.
- **D1/DA-2:** parcialmente resuelto → moderar el rótulo β-CMC pero **aportar la correlación β-PAC↔β-CMC (ρ=0.55)**
  como evidencia de un componente corticomuscular real; declarar el caveat de limpieza.
- **Siguiente paso opcional (para fortalecer #2):** recomputar CMC sobre EEG mínimamente limpiado (pre-ICA/
  pre-iCanClean) y ROI central → magnitud justa. Si ahí la CMC β sube y sigue correlacionando, el stream β-CMC
  pasa de "componente parcial" a "hallazgo".
