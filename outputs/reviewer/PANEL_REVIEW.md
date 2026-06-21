# PANEL REVIEW — Paper2 V1 (re-submission tras R1) · 2026-06-20
Skill: academic-paper-reviewer (modo full) · ejecución inline · manuscrito revisado en `Paper2_v1_overleaf/`.
Read-only: este informe NO modifica el manuscrito. Todo se ancla a pasajes/números reales.

================================================================================
## FASE 0 — Análisis de campo y configuración del panel
================================================================================
- **Disciplina primaria:** neurociencia cognitiva / electrofisiología humana (EEG-EMG).
- **Secundaria:** cognición corporizada / sensorimotor-cognición; coherencia corticomuscular.
- **Paradigma:** experimental, entre-sujetos (n=31 chew / 15 control), tarea 2-back.
- **Tipo metodológico:** EEG inducido + TF (Morlet/CBPT), specparam/FOOOF, PAC (Tort MI + doble null), circular stats.
- **Tier de revista objetivo (estimado):** tier-2 sólido — *NeuroImage*, *eNeuro*, *Psychophysiology*,
  *Journal of Neuroscience Methods* (por la contribución metodológica de EEG durante masticación). Scientific
  Reports (wlscirep.cls actual) es coherente.
- **Madurez:** re-submission tras una R1 sustantiva; análisis canónico recomputado.

**Panel configurado (5 + EIC):**
| Rol | Identidad simulada | Foco |
|---|---|---|
| EIC | Editor de revista neuro tier-2, sensible a sobre-claims | fit, significancia, originalidad, tono |
| R1 Metodología | Electrofisiólogo cuantitativo (stats EEG, CFC, poder) | diseño, poder, PAC, FOOOF, confounds |
| R2 Dominio | Investigador de oscilaciones/WM + CMC | marco teórico, literatura, β-CMC, aperiódico |
| R3 Perspectiva | Cognición corporizada / MoBI / ritmos corporales | impacto, naturalismo, "¿y qué?" cognitivo |
| DA | Abogado del diablo | reto al argumento central, falacias, contraejemplo |

================================================================================
## FASE 1 — Cinco informes independientes
================================================================================

### [EIC] Editor-in-Chief
**Resumen.** Contribución doble: (1) metodológica — recuperar EEG válido durante masticación activa con
controles de artefacto serios (AAFT, topografía, condición); (2) empírica — masticación rítmica acompaña
mejora conductual intra-grupo, aplanamiento aperiódico, theta frontal-medial y acoplamiento fase-masticatoria.
La re-submission respondió a R1 con honestidad poco común: retractó una correlación previa contradicha,
reportó nulos, y declaró límites. Esto eleva mucho la credibilidad.
**Fortalezas.** Diseño con control de práctica + manipulation check (masetero +6.1 dB, 31/31); doble null +
control negativo (S8, 1/15) → mecanismo limpio; pipeline reproducible (config único, seeds).
**Preocupaciones de nivel editorial.**
- (E1) **Brecha título/abstract ↔ evidencia.** El título ("couples to frontocentral theta and broadband
  cortical dynamics") es defendible, pero el *gancho cognitivo* (memoria de trabajo) descansa sobre un efecto
  conductual DIFERENCIAL que es sólo tendencia (Grupo×Bloque p=0.054) y SIN enlace individual PAC↔conducta
  (n.s.). El lector espera "masticar ayuda a la WM"; lo que se demuestra es "masticar cambia el EEG y mejora
  el desempeño dentro de quienes mastican, sin vínculo neuro-conductual individual". Hay que cerrar esa brecha
  en abstract/intro para no prometer de más.
- (E2) **Significancia.** Para tier-2 alcanza; para tier-1 (Nature Comms/JNeurosci) faltaría el eslabón causal
  o el anclaje conductual. Encaja mejor en revista metodológica/oscilaciones.
**Recomendación EIC:** *Minor-to-Major Revision* (ver síntesis).

### [R1] Revisor de Metodología
**Fortalezas.** Doble null (circular + AAFT) es el estándar correcto para PAC y desactiva el artefacto
no-sinusoidal; el control negativo S8 completa un 2×2 elegante; Wilcoxon en lugar de t para RT/IES sesgados
es correcto; specparam separa pico θ del 1/f.
**Debilidades (con localización y fix):**
- (M1, MAJOR) **Poder del grupo control (n=15) para la interacción.** La interacción Grupo×Bloque (p=0.054) y
  Grupo×Condición (TF p=0.010) descansan sobre 15 controles. Con n=15 el test de interacción está sub-potenciado;
  un p=0.054 "tendencia" puede ser tanto efecto real débil como ruido. *Fix:* reportar potencia observada /
  análisis de sensibilidad (qué tamaño de efecto detectaría el diseño al 80%), e intervalos de confianza del
  contraste diferencial, no sólo p. (Results §conducta; Methods §stats.)
- (M2, MAJOR) **Confound de orden fijo NoChew→Chew.** Ambos grupos hacen B1=NoChew, B2=(Chew|NoChew). La mejora
  intra-casos confunde masticación con práctica/tiempo-en-tarea; el grupo control lo mitiga pero no lo elimina
  (su no-mejora podría deberse a fatiga que cancela la práctica). *Fix:* declarar explícitamente y, si hay datos,
  modelar tendencia temporal dentro de bloque (¿la mejora es escalonada al iniciar el chew o gradual?). (Ya está
  como limitación; subir su prominencia.)
- (M3, MINOR) **MI Tort con N=18 bins y datos continuos:** ok, pero reportar longitud de datos por sujeto y su
  relación con MI (MI sesga con pocos ciclos). El f_chew±0.5 Hz individual es defendible; confirmar que el ancho
  de banda de amplitud (4–7) no solapa con el de fase tras el ±0.5. (Methods §PAC.)
- (M4, MINOR) **FOOOF:** rango de ajuste 3–40 Hz, modo fixed (sin knee). Reportar R² del ajuste y error por
  condición; un aplanamiento del exponente puede confundirse con peor ajuste si el R² cae en chew. (Methods §spectral.)
- (M5, MINOR) **Comparaciones múltiples:** se declara "ROI/banda/ventana a priori → sin corrección". Aceptable,
  pero hay muchos tests (bandas θ/α/β × ventanas × grupos × Rayleigh × correlaciones). Aportar una tabla maestra
  de todos los contrastes con su estatus a priori vs exploratorio.
**Rigor score (0–100):** 72 (sólido; limitado por poder de interacción y orden fijo, ambos declarados).

### [R2] Revisor de Dominio (oscilaciones / WM / CMC)
**Fortalezas.** El giro al marco DUAL es teóricamente más honesto que el θ-PAC monolítico previo. La cobertura
del aperiódico-E/I (Gao, Donoghue, Voytek) y del precedente de ritmos periféricos (Zelano, Tort respiración)
está bien elegida y soporta que "PAC a ritmo periférico" no es artefacto.
**Debilidades:**
- (D1, MAJOR) **El rótulo β-PAC "análogo a CMC" hace trabajo conceptual pesado con respaldo indirecto.** Se
  invoca CMC (coherencia fase-fase) por analogía para un PAC fase-amplitud. Aunque el hedge es correcto, el
  manuscrito no mide CMC clásico (coherencia EEG-EMG en β) que SÍ podría calcularse con estos datos (hay EMG).
  *Fix decisivo:* calcular la coherencia corticomuscular β real (EEG-EMG) y mostrar que coexiste / correlaciona
  con el β-PAC. Eso convertiría una analogía en evidencia y blindaría el stream β. Sin ello, un revisor de CMC
  dirá que la "firma corticomuscular" es especulativa. (Discu §dual reading.)
- (D2, MAJOR) **El stream θ-cognitivo pierde su ancla funcional.** La especificidad θ se sostiene por topografía
  + condición + pico FOOOF, pero NO por función (PAC×conducta n.s.). Entonces ¿en qué sentido es "ejecutivo"?
  La única conexión WM es la co-localización con literatura fm-theta. *Fix:* o bien anclar θ a la conducta con
  un análisis más sensible (p.ej. θ *power* late × IES, no PAC×IES), o moderar "executive/working-memory" a
  "fronto-medial theta, una banda asociada a WM en la literatura", sin atribuir función demostrada aquí.
- (D3, MINOR) **Faltan referencias clave de CMC oromotor y de PAC-a-periferia** que el propio equipo ya tiene
  staged (biblio_new.bib: Glories2021/2023, Ishii2016, vanWijk2012 — varias ya citadas; verificar que Kristeva,
  Conway estén en Discu y no sólo en methods). Considerar Bourguignon/Piitulainen (CMC) si aplica.
- (D4, MINOR) θ=4–7 con guarda 7–8: la decisión está justificada en methods, pero reportar el resultado con
  4–8 en suplementario reforzaría robustez (anticipar al reviewer escéptico de la banda).
**Contribución al dominio score:** 70 (marco coherente; el β-CMC pide medición directa, no analogía).

### [R3] Revisor de Perspectiva (cognición corporizada / MoBI)
**Fortalezas.** El framing naturalista (EEG durante movimiento orofacial real) es valioso y poco común;
el precedente respiratorio sitúa la masticación en una clase amplia de ritmos cuerpo→córtex. Buen potencial
de impacto metodológico (speech, deglución, respiración).
**Debilidades:**
- (P1, MAJOR) **Prueba del "¿y qué?".** Si el efecto conductual diferencial es n.s. y no hay enlace
  PAC↔conducta, la pregunta "¿masticar mejora la cognición?" queda sin responder; lo que queda es "masticar
  deja una huella cortical broadband + θ". Es interesante, pero el manuscrito debe decidir si su contribución
  es (a) metodológica (recuperar EEG en masticación) — y entonces liderar con eso — o (b) sustantiva
  (masticación-cognición) — y entonces aceptar que la evidencia es parcial. Actualmente promete (b) y entrega
  sobre todo (a)+huella neural.
- (P2, MINOR) **Generalización.** Adultos sanos jóvenes, n total 46, control 15. Útil declarar a qué población
  NO se generaliza (mayores, disfunción masticatoria) — ya está, mantener.
- (P3, MINOR) **Implicación práctica/ética:** mínima discusión de aplicabilidad (¿chewing como intervención
  cognitiva?). Dado el nulo conductual diferencial, evitar insinuar aplicaciones; el riesgo de sobre-venta
  mediática es real con "chewing gum makes you smarter".
**Impacto/originalidad score:** 68 (originalidad metodológica alta; contribución sustantiva acotada).

### [DA] Abogado del Diablo
**Contra-argumento más fuerte (≈250 palabras).**
El manuscrito vende una historia sensorimotor→cognición, pero su propia recomputación canónica la desarma en
los dos eslabones que importarían: (1) el beneficio conductual *atribuible a masticar* (vs práctica) NO es
significativo (Grupo×Bloque p=0.054); (2) el acoplamiento neural NO predice la conducta a nivel individual
(todos los ρ n.s.). Lo que queda demostrado es compatible con una lectura mucho más mundana y deflacionaria:
masticar es un acto motor rítmico que (i) introduce un acoplamiento corticomuscular broadband dominado por β
—exactamente lo que la literatura de CMC espera de cualquier salida motora sostenida— y (ii) eleva el arousal/
excitabilidad global (aplanamiento 1/f), lo que infla potencia en varias bandas incluida θ. Bajo esta lectura,
el "θ frontal-medial ejecutivo" no es una firma cognitiva sino el reflejo banda-θ de un cambio de excitabilidad
broadband co-localizado por casualidad con el generador fm-theta; y el "stream β-CMC" no es un hallazgo sino
el nombre noble del artefacto motor que el resto del paper se esfuerza por excluir. La distinción cortical vs
muscular se apoya en AAFT + topografía, que descartan *leakage miogénico directo* pero NO descartan que el
efecto sea *motor-cortical* (drive descendente) más que *cognitivo*. En suma: el paper podría estar describiendo,
con instrumentación excelente, la neurofisiología de mover la mandíbula —no de pensar mejor mientras se mueve.

**Issue list:**
- (DA-1, **CRITICAL**) **Disociación lectura cognitiva vs motora no resuelta.** Nada en los datos separa
  "θ aumenta porque sirve a la WM" de "θ aumenta porque la excitabilidad global subió por el acto motor".
  El nulo PAC↔conducta y el nulo dosis-respuesta favorecen la lectura motora/arousal. *Requisito:* el manuscrito
  debe (a) NO afirmar función cognitiva del θ como conclusión, y (b) ofrecer al menos un test que disocie
  (p.ej. ¿el θ-late predice IES aunque el PAC no? ¿la mejora intra-casos correlaciona con Δθ o Δexponente?).
  Si ninguno separa, reformular la contribución como neurofisiológica/metodológica.
- (DA-2, MAJOR) **β-CMC como reetiquetado del confound.** El paper excluye EMG como artefacto y luego celebra
  el β-PAC como "corticomuscular". O el β-coupling es señal motora genuina (entonces medir CMC clásico y
  abrazarlo, D1) o es residuo motor (entonces no es un "hallazgo"). Elegir y comprometerse.
- (DA-3, MAJOR) **Cherry-picking de ventana/banda para θ.** θ se vuelve "el relevante" sólo tras seleccionar
  ventana late + topografía + FOOOF; en magnitud θ es el MÁS DÉBIL (β≫α≫θ). Riesgo de narrar a posteriori la
  banda que conviene. *Fix:* declarar explícitamente que θ no es el efecto dominante y que su relevancia es
  inferida de la literatura, no de la fuerza del efecto.
- (DA-4, MINOR) **"Genuine theta peak 97%"** es un buen control pero no implica relevancia funcional; no usarlo
  como evidencia de que θ "importa" para la cognición.
**Alternativas ignoradas:** arousal global inespecífico (no ejecutivo); drive motor descendente; movimiento
ocular/postural co-variando con masticación.
**Stakeholders ausentes:** lector clínico (¿implicaciones?), modelador (direccionalidad), escéptico de CFC.

================================================================================
## FASE 2 — Síntesis editorial y decisión
================================================================================

### Matriz de consenso
| Tema | EIC | R1 | R2 | R3 | DA | ¿Consenso? |
|---|---|---|---|---|---|---|
| Rigor metodológico / controles de artefacto | + | + | + | + | + | **Consenso fuerte (fortaleza)** |
| Honestidad de la re-submission (nulos declarados) | + | + | + | + | ~ | Consenso |
| Efecto diferencial conductual sólo tendencia (poder/orden) | E1 | M1,M2 | — | P1 | DA-1 | **Consenso (debilidad MAJOR)** |
| Sin enlace PAC↔conducta ⇒ ancla cognitiva débil | E1 | — | D2 | P1 | DA-1 | **Consenso (debilidad MAJOR/CRITICAL)** |
| β-CMC por analogía, no medido | — | — | D1 | — | DA-2 | Consenso parcial (MAJOR) |
| Selección post-hoc de θ (banda más débil) | — | M5 | D4 | — | DA-3 | Consenso parcial |
| Contribución: ¿metodológica vs sustantiva? | E2 | — | — | P1 | DA-1 | Consenso (posicionamiento) |

### Issue CRÍTICO (DA-1) — IRON RULE: la decisión NO puede ser Accept.
La disociación lectura motora/arousal vs cognitiva no está resuelta y los nulos canónicos favorecen la lectura
deflacionaria. No invalida el paper, pero obliga a (a) moderar las afirmaciones funcionales y (b) intentar al
menos un test disociativo, antes de aceptar.

### DECISIÓN EDITORIAL: **MAJOR REVISION**
Manuscrito metodológicamente sólido y notablemente honesto, pero su narrativa cognitiva excede lo que los datos
canónicos sostienen. Es publicable en revista tier-2 tras (1) realinear claims con evidencia y (2) fortalecer o
reposicionar el eslabón β-CMC y el eslabón cognitivo.

### REVISION ROADMAP (priorizado)
**P1 — Bloqueantes (resolver para aceptar):**
1. **Disociar o desclamar (DA-1, E1, D2, P1).** Añadir test(s) que separen θ-cognitivo de arousal/motor:
   p.ej. Δθ-power-late × ΔIES intra-casos; Δexponente × ΔRT; ¿el θ-late predice conducta aunque el PAC no?
   Si ninguno disocia, reescribir abstract/intro/discusión para que la contribución principal sea
   neurofisiológica/metodológica y el vínculo cognitivo quede como abierto, no insinuado.
2. **β-CMC: medir o moderar (D1, DA-2).** Calcular coherencia corticomuscular β clásica (EEG-EMG) y mostrar
   coexistencia con β-PAC; o, si no, bajar el rótulo a "broadband coupling dominado por β, compatible con drive
   motor" sin invocar "corticomuscular signature" como hallazgo.
3. **Poder e incertidumbre de la interacción (M1, M2).** Reportar IC del contraste diferencial, análisis de
   sensibilidad de potencia (n=15), y prominencia de la limitación de orden fijo.
**P2 — Importantes:**
4. Declarar explícitamente que θ es la banda más débil en magnitud y que su relevancia es inferida (DA-3, M5);
   tabla maestra de contrastes a priori vs exploratorios.
5. FOOOF: reportar R²/error por condición para descartar que el aplanamiento sea peor ajuste (M4).
6. Realinear título/abstract con la contribución real; evitar insinuar aplicaciones (P3).
**P3 — Menores:**
7. Longitud de datos por sujeto × MI; solape de bandas tras f_chew±0.5 (M3).
8. Robustez θ 4–8 en suplementario (D4); verificar refs CMC en Discusión (D3).

### Nota de calibración del panel
Puntajes orientativos (0–100): Rigor 72 · Dominio 70 · Impacto 68 · (Novedad metodológica ~80, Novedad
sustantiva ~55). Coherente con Major Revision en tier-2. La honestidad de la re-submission es un activo real y
debe preservarse: la ruta correcta NO es inflar claims, sino alinear el marketing del paper con su evidencia.
