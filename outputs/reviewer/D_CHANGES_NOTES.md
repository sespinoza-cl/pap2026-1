# D_CHANGES_NOTES — Manuscrito Paper2 V1: cambios aplicados (Fase D) · 2026-06-20

Marco DUAL θ-cognitivo + β-CMC (β-PAC **análogo** a CMC, no "=CMC"). Números CANÓNICOS
(RESULTS_CANONICAL.md). Diseño SANO sin TMD. Edición IN PLACE en `Paper2_v1_overleaf/Text_parts/`.

## TRACK CHANGES (vía elegida: latexdiff)
Snapshot de la versión enviada en `Paper2_v1_overleaf/Text_parts_submitted/` (copia intacta).
Para generar el PDF con marcas:
```
latexdiff Text_parts_submitted/3_results.tex Text_parts/3_results.tex > diff_results.tex
```
o, a nivel de documento, compilar `main.tex` con includes apuntando a cada versión y correr
`latexdiff` sobre los .tex aplanados. Se mantiene además la versión "clean" (la actual, sin marcas).
(Alternativa `changes` package descartada para no ensuciar 5 archivos con \added/\deleted a mano.)

## ARCHIVOS EDITADOS
- **main.tex** — título → "Rhythmic mastication couples to frontocentral theta and broadband cortical
  dynamics during working memory" (dual, sin "entrains via θ-PAC").
- **1_Abstract.tex** — reescrito: dual; quita "coupling strength scaling with chewing frequency" y
  "enhanced θ"; añade aperiódico, β-CMC, control negativo; Wilcoxon.
- **2_intro.tex** — θ 4–8→4–7 (×2); **lógica de falsación reescrita**: el control de artefacto ya NO
  descansa en frecuencia (PAC broadband β-dominante) sino en topografía + AAFT + condición.
- **3_results.tex** — REESCRITO completo y reordenado (manip-check → conducta → TF/topo → FOOOF →
  PAC θ+banda → control negativo). Números canónicos. ELIMINADO el ρ=−0.365 θ×IES (contradicho).
  Quitado "dose-response"/"enhancement"/"recovery". Heading concatenado corregido (Minor 3).
- **4_discu.tex** — apertura reescrita (sin scaling, sin co-ocurrencia conductual sobre-vendida);
  añadida subsección DUAL "corticomuscular beta coupling and executive theta" (β análogo CMC + aperiódico
  E/I + respiración como precedente); muscle-control reescrito (topo+AAFT+condición); Limitations
  ampliadas (tendencia p=0.054, broadband, PAC×conducta n.s., sin crossover, CMC por analogía);
  "two independent experiments" → convergencia cross-study con Espinoza2025.
- **5_methods.tex** — PAC primaria: f_chew individual ±0.5 Hz (no MI-optimizado), θ 4–7, bandas control;
  añadido doble null (circular+AAFT, 500 surr) + subsección **Manipulation check** + subsección
  **Condition specificity and negative control**; EMG-validation reescrita (sin MI 1.66e-4 stale,
  sin lógica de frecuencia; topo+AAFT+condición).
- **6_figs.tex** — REESCRITO: 5 figuras con los paneles nuevos y captions canónicos/duales/sin-TMD.
  Labels: fig:beh, fig:tf, fig:topo, fig:fooof, fig:cfc (todas referenciadas en results).
- **7_others.tex** — Data/Code availability reforzado (config único, seeds, provenance).

## BIBLIOGRAFÍA
Añadidas a `biblio.bib` (claves usadas): MarisOostenveld2007, Gao2017, Gerber2016, Conway1995,
Kristeva2007, vanWijk2012 (resueltas por DOI vía opencite). Voytek2015→**VoytekKnight2015** (ya existía).
Aru2015→**ARU2015**, Zerbi2019→**zerbi2019** (corrección de mayúsculas). Verificado: **0 citas sin resolver,
0 \ref de figura sin label** en todo el manuscrito.
Refs nuevas adicionales disponibles en `outputs/reviewer/biblio_new.bib` (40, para la discusión si se amplía).

## FIGURAS
PNGs nuevos copiados a `Paper2_v1_overleaf/p2_Figs/` (16 paneles principales). El usuario compone/ajusta
en LaTeX según su preferencia (paneles individuales). Suplementarios en `Analysis_V1_Final/outputs/figures/`
(topoS_theta_*, artifactS_*, doseS_*, pacS_MI_by_window, comodulogramS_2D_*). Ver FIGURE_INDEX.md.

## PENDIENTE (no bloqueante, requiere decisión/recurso del usuario)
- Suppl (8_suppl.tex) NO reescrito a fondo: revisar que sus números/figuras suplementarias apunten a
  los .mat canónicos y a los paneles `*S_*`.
- Paneles suplementarios diferidos: topoS α/θ-α (recomputo per-canal), artifactS_sharpness, fooof_dissociation
  (ver FIGURE_INDEX.md "DIFERIDOS").
- Correr `latexdiff` y compilar en Overleaf para el PDF con marcas + versión clean.
