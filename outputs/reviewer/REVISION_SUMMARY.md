# REVISION SUMMARY + RE-REVIEW — Paper2 V1 · 2026-06-20
Verificación (modo re-review) del manuscrito revisado contra la decisión del panel (Phase 2 = MAJOR REVISION,
ver PANEL_REVIEW.md). Matriz de trazabilidad: crítica → acción → evidencia/ubicación → ¿verificado? → residual.
Read-only sobre el manuscrito; ancla a archivos reales.

## R&R TRACEABILITY MATRIX

| # (sev) | Crítica del panel | Acción tomada | Evidencia / ubicación | Verificado | Residual |
|---|---|---|---|---|---|
| **DA-1 (CRITICAL)** | Lectura cognitiva vs motor/arousal sin disociar; ¿algún correlato neuro-conductual? | (a) Barrido disociativo (~44 corr) confirmó nulo entre-sujetos; (b) **análisis single-trial**: θ-late frontal-medial ↔ RT (LMM −41.8 ms, p=0.0007; 23/31; Fisher-z p=0.0016) | Results §"single-trial"; Fig.6 `fig:bridge`; `rev_singletrial.py`; REVISION_ANALYSES.md #1 | **SÍ (sustancial)** | Link es post-respuesta (estado, no prospectivo); dirección no resuelta — declarado |
| **DA-2 / D1 (MAJOR)** | β-PAC = ¿CMC genuina o artefacto reetiquetado? "medir o moderar" | **Medido**: topografía CMC β → difusa, pico temporal T8, central≤temporal (p=0.05) = leakage. **Pivote de marco**: β→leakage, no CMC | Discu §"broadband beta…residual EMG"; Methods §CMC control; Suppl `fig:betacmc`; `rev_cmc_topography.m`; REVISION_ANALYSES #2b | **SÍ (resuelto)** | EEG de cuero cabelludo no separa 100% volumen-conducción — declarado |
| **DA-3 (MAJOR)** | Cherry-picking de θ (banda más débil en magnitud) | Texto declara explícitamente θ = más débil en magnitud; relevancia por topografía+condición+FOOOF+**link trial-level**; β reatribuido a leakage | Results §PAC; Discu; Limitations | **SÍ** | — |
| **D2 / E1 (MAJOR)** | θ-cognitivo sin ancla funcional; gap título/abstract↔evidencia | Ancla funcional aportada (single-trial θ↔RT); título pivotado (θ + aperiódico, sin "broadband/β"); abstract realineado | main.tex título; 1_Abstract; Results | **SÍ (sustancial)** | El ancla es trial-level/estado, no causal — declarado |
| **P1 (MAJOR)** | "¿y qué?" metodológico vs sustantivo | Contribución clarificada: método (EEG en masticación) + θ-cortical genuino con correlato conductual trial-level + aperiódico (cascada cross-study con Paper 1) | Abstract; Discu; Suppl `fig:aperiodic_crossstudy` | **SÍ** | Beneficio diferencial sigue siendo tendencia |
| **M1 (MAJOR)** | Potencia de n=15 para la interacción | **CERRADO**: IC del contraste diferencial + sensibilidad. RT ΔΔ=−58.9 ms CI[−133.5,+15.8] d=−0.60; IES d=−0.50; MDES@80%=d 0.90 (diseño solo detecta interacciones grandes) → la tendencia es límite de potencia, no ausencia de efecto | Results §conducta; Limitations; Methods §stats | **SÍ (resuelto)** | Crossover/mayor n lo zanjaría (declarado) |
| **M2 (MAJOR)** | Confound orden fijo NoChew→Chew | Declarado como limitación; controlado por grupo no-masticador; prominencia subida | Limitations; Results | **SÍ (declarado)** | No eliminable sin crossover — declarado |
| **M3 (MINOR)** | MI: longitud de datos × MI; solape banda tras f_chew±0.5 | Methods describe f_chew±0.5 individual (no MI-optimizado); doble null | Methods §PAC | PARCIAL | Reportar n_trials/longitud por sujeto si lo pide el editor |
| **M4 (MINOR)** | FOOOF R²/error por condición | Offline R²=0.99 reportado; online ajuste estándar | REVISION_ANALYSES; rev_paper1_aperiodic | PARCIAL | Añadir R² online por condición al Suppl |
| **M5 / DA-3 (MINOR)** | Comparaciones múltiples / tabla de contrastes | θ 4–7 a priori declarado; nulos reportados; corrección discutida | Methods §stats | PARCIAL | Tabla maestra a priori vs exploratorio (opcional) |
| **D3 (MINOR)** | Refs CMC/PAC-periferia | biblio +6 (Conway, Kristeva, vanWijk, Gerber, Gao, Maris); biblio_new.bib 40 | biblio.bib | SÍ | — |
| **D4 (MINOR)** | Robustez θ 4–8 | Justificado en Methods (4–7 preserva interacción midline) | Methods §Minor1 (carta) | PARCIAL | Suppl con 4–8 (opcional) |
| **DA-4 (MINOR)** | Pico genuino ≠ relevancia funcional | No se usa el 97% como evidencia de función; la función viene del link trial-level | Results/Discu | SÍ | — |

## NUEVA DECISIÓN EDITORIAL (re-review)
**De MAJOR REVISION → MINOR REVISION.** Los tres bloqueantes MAJOR/CRITICAL quedaron sustancialmente resueltos
con datos, no con retórica:
- **DA-1** (el crítico): ahora hay un puente conductual↔EEG real (single-trial θ-late↔RT) + explicación honesta
  del nulo entre-sujetos (fiabilidad/potencia). El IRON RULE (DA CRITICAL ⇒ no Accept) se levanta porque la
  disociación dejó de estar "sin resolver": hay vínculo a nivel de ensayo, con caveats declarados.
- **DA-2/D1**: el β-CMC se midió y se reinterpretó como leakage — la crítica se volvió un hallazgo del propio
  paper (control topográfico). Ya no hay overclaim.
- La honestidad del conjunto (nulos declarados, pivote de marco basado en datos, caveats explícitos) es ahora
  una fortaleza saliente.

## RESIDUALES PARA CERRAR (prioridad)
1. ~~M1 IC + sensibilidad~~ **CERRADO 2026-06-20** (ver fila M1).
2. **M4 — fácil (opcional):** añadir R²/error del ajuste FOOOF online por condición al Suppl.
3. **Menores opcionales:** tabla maestra de contrastes (M5); Suppl θ 4–8 (D4); n_trials por sujeto (M3).

Ningún residual es bloqueante; todos son menores u opcionales.

## VEREDICTO (actualizado, M1 cerrado)
Los tres MAJOR/CRITICAL y el último MAJOR parcial (M1) están resueltos con datos. Posición: **MINOR REVISION**,
sin residuales bloqueantes. La narrativa final es coherente y honesta:
método (EEG en masticación) → **θ-cortical genuino con correlato conductual trial-level** → **β = leakage**
(control topográfico) → **aperiódico = excitabilidad** (cascada cross-study online↔offline).
Provenance del cierre M1: `code/` cálculo CI/d/MDES sobre `v1_S1_behavior_stats.mat` (TTestIndPower).
