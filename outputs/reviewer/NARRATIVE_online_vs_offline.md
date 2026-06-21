# NARRATIVE — Por qué theta "brilla" OFFLINE (Paper 1) y no ONLINE (Paper 2)
2026-06-20 · Reconciliación Espinoza et al. 2025 (Sci Rep, pre-task chewing) ↔ Paper2 V1 (online chewing).

## 0. El aparente puzzle (lo que el usuario esperaba)
Paper 1 (offline): masticar 1 min ANTES de la tarea → theta frontocentral ↑ + conectividad frontoparietal
(PLI) ↑ + conducta ↑; theta como "marca" task-específica de WM. Misma tarea (2-back) masticando EN LÍNEA
(Paper 2) → se esperaba lo mismo, pero theta sale **broadband, β-dominante, y desacoplada de la conducta**.

## 1. El hecho que DESACTIVA la contradicción (clave)
**En AMBOS estudios theta sube, y en NINGUNO theta predice la conducta a nivel individual.**
- Paper 1 (PMC12698685): theta power/PLI ↑ post-chewing, PERO "no direct correlation between theta
  power/connectivity and behavioral metrics (RT/accuracy)". Lo único que correlaciona es theta-peak-freq ×
  chewing-freq (tuning de frecuencia), NO theta-power × desempeño.
- Paper 2 (canónico): theta ↑, PERO Δθ/Δexponente/PAC/pico/amplitud × conducta = **0/44 tras corrección**.
→ La "marca de WM" de Paper 1 nunca fue una correlación individual theta↔conducta; fue (a) un AUMENTO
task-específico de theta y (b) un tuning de frecuencia. Paper 2 reproduce el aumento de theta y el tuning
(ρ=0.324 chewing-freq × theta-peak). **No hay, pues, una capacidad predictiva perdida entre offline y online:
nunca existió.** Esto ya tranquiliza: los dos papers son consistentes en el punto que importa.

## 2. Las diferencias REALES offline vs online (tres ejes)
| Eje | OFFLINE (Paper 1) | ONLINE (Paper 2) |
|---|---|---|
| **Ventana** | temprana 150–550 ms (encoding/atención) | tardía 900–1300 ms (mantenimiento) |
| **Contexto espectral** | theta limpio + PLI frontoparietal (banda-específico) | broadband, β el más fuerte, aperiódico aplanado |
| **Acto motor** | TERMINADO (gap >60 s, sin EMG/coupling) | CONCURRENTE (β-PAC, EMG, coupling corticomuscular) |
| **Conducta** | beneficio claro | beneficio diferencial sólo tendencia (p=0.054) |
| **Anestesia** | sin efecto → vía CENTRAL (no sensorial periférica) | (no testeada) |

## 3. El mecanismo unificador: UNA cascada temporal, dos ventanas
**Semilla común (excitabilidad):** masticar eleva la excitabilidad cortical. Evidencia: el aplanamiento del
exponente aperiódico (Paper 2, χ 1.110→0.735) es el índice agudo de ese cambio E/I (Gao 2017; Donoghue 2020;
Miskovic 2024). **MATIZ (auditoría Gemini):** la anestesia de Paper 1 fue **tópica gingival** → bloquea el
tacto superficial pero **NO** la propiocepción profunda (husos musculares maseterinos, mecanorreceptores
periodontales profundos, de integración bilateral en tronco). Por tanto el efecto-nulo de anestesia sólo
descarta una contribución **táctil-superficial**; la aferencia trigeminal **profunda rítmica sigue intacta**.
NO se puede afirmar un driver "puramente central": la vía trigeminal periférica profunda no está excluida.
Redacción correcta: "la facilitación no depende del tacto gingival superficial, consistente con una
contribución central y/o propioceptiva profunda".

**ONLINE (durante la masticación, Paper 2):** el acto motor rítmico CONCURRENTE domina el EEG e impone:
1. **Acoplamiento corticomuscular broadband** (β-PAC el más fuerte; β-PAC covaría con CMC β real ρ=0.55,
   banda-específico) → una firma motora que offline está ausente.
2. **Excitabilidad global** (aperiódico↓) que infla potencia en VARIAS bandas, theta incluida.
→ La theta cognitiva queda **embebida y enmascarada** en un estado broadband motor+arousal, no aislable como
"theta de WM". Además, la demanda motor-cognitiva concurrente **compite por los mismos recursos
fronto-mediales** (interferencia tipo supresión articulatoria; Baddeley 1986; y la CMC se degrada bajo carga
dual, Omejc 2025 — la competencia es recíproca). Esa competencia es la razón plausible de que el beneficio
diferencial sea sólo tendencia (el costo motor concurrente compensa parte del beneficio de arousal).

**OFFLINE (tras la masticación, Paper 1):** el acto motor TERMINÓ; desaparecen EMG, β-PAC y el coupling
corticomuscular broadband. Lo que **persiste** es la cola de excitabilidad/arousal sembrada por masticar
(estado E/I más excitable, tono LC-NE). En esa ventana "limpia", sin motor concurrente, ese estado residual
se expresa como un realce theta **focalizado, temprano y task-específico** con integración frontoparietal
(PLI) → theta "brilla" precisamente porque ya no está enterrada bajo el coupling motor concurrente.

**En una frase:** *online captamos el ACTO (motor broadband + arousal, theta enmascarada y motoramente
competida); offline captamos la HUELLA (sólo la cola de excitabilidad central, que aflora como theta limpia
de WM).* Dos ventanas de la misma cascada masticación→excitabilidad central→theta, separadas por la presencia
o ausencia del motor concurrente.

## 4. Predicciones falsables (fortalecen el relato; varias ya confirmadas)
- (a) ONLINE: el coupling debe ser máximo en β/broadband y la theta-especificidad mínima. **Confirmado**
  (β≫α≫θ; Friedman p<0.0001).
- (b) El aperiódico aplanado debe ser el hilo común; offline debería persistir aunque el coupling motor no.
  **Parcial** (medido online; Paper 1 no reportó aperiódico → análisis pendiente sobre datos de Paper 1).
- (c) OFFLINE la theta debe ser temprana + frontoparietalmente integrada, sin CMC concurrente. **Confirmado**
  (150–550 ms + PLI; gap >60 s sin EMG).
- (d) Remover el acto motor (offline) debe quitar el enmascaramiento y revelar la theta cognitiva. **Confirmado**
  por contraste entre los dos diseños.
- (e) **CONFIRMADO (2026-06-20):** se computó el exponente aperiódico en los datos de Paper 1 (offline,
  n=30, FOOOF, ROI frontocentral, chew vs no-chew durante la tarea). **El aperiódico TAMBIÉN se aplana
  offline:** χ chew=1.352 vs nochew=1.445, Δ=−0.093, 21/30 sujetos, Wilcoxon two-sided p=0.005
  (one-sided chew<nochew p=0.0025), R²=0.99. Online (Paper 2) Δ=−0.375 p<0.0001. **Mismo signo, magnitud menor
  offline** → coherente con una cola de excitabilidad que es MÁXIMA durante el acto (online) y PERSISTE atenuada
  tras masticar (offline). Esto ancla empíricamente la "semilla común" de la cascada y, además, es un análisis
  NUEVO que cross-valida Paper 1 (que no reportó aperiódico). Script: `code/rev_paper1_aperiodic.py` →
  `outputs/stats/rev_paper1_aperiodic.npz`. *Caveat:* muestras/pipelines distintos → no comparar magnitudes
  absolutas; el punto robusto es el MISMO SIGNO significativo en ambos.

## 5. Cómo usarlo en el manuscrito (Discusión)
Convertir la disociación online/offline en un APORTE, no en un problema:
"La masticación siembra un estado de excitabilidad cortical central (aperiódico aplanado, vía trigémino-LC,
independiente de aferencia periférica; Espinoza 2025). Durante el acto, ese estado coexiste con un
acoplamiento corticomuscular broadband que enmascara la theta cognitiva e impone un costo motor-concurrente,
de modo que el beneficio conductual es sólo tendencia y la theta no es task-selectiva. Una vez cesado el acto
(Espinoza 2025), la cola de excitabilidad persiste sin el coupling motor y aflora como la theta frontoparietal
temprana y task-específica del paradigma offline. Los dos estudios son así ventanas complementarias de una
única cascada sensorimotor-cognitiva." Mantener "consistente con", no "demuestra" (estilo Wael).

## 5b. EXPLICACIÓN ALTERNATIVA (obligatoria por honestidad — auditoría Gemini)
La cascada de una-sola-theta-enmascarada es elegante PERO hay una alternativa más simple que NO debe omitirse:
**las ventanas medidas son fases cognitivas distintas** (offline 150–550 ms = encoding/atención temprana;
online 900–1300 ms = mantenimiento). Podría no ser "la misma theta enmascarada" sino **dinámicas de red
distintas** propias de cada fase, más el carácter dual-task (motor+cognitivo) del online. Es decir: la theta
temprana de encoding (offline) y la theta tardía de mantenimiento (online) podrían ser fenómenos parcialmente
diferentes, y la "cascada única" sería una sobre-interpretación. **Cómo decidir entre ambas:** (i) el test del
aperiódico offline (§4e) — si el aperiódico aplanado aparece en AMBOS, apoya la semilla común; (ii) reportar la
theta online en la MISMA ventana temprana (150–550 ms) que Paper 1 — si online la theta temprana también es
broadband/desacoplada, refuerza el enmascaramiento; si es limpia, apoya la alternativa de fases distintas.

## 6. Caveats honestos
- Correlacional; la cascada enlaza DOS muestras/diseños distintos (n=30 offline factorial vs n=31/15 online).
- **No afirmar driver "puramente central"** (la propiocepción profunda trigeminal sigue intacta; §3 MATIZ).
- **Alternativa viva:** ventanas = fases cognitivas distintas (§5b); la "cascada única" es hipótesis, no hecho.
- La CMC β online está atenuada por la limpieza de artefactos (límite inferior); la correlación β-PAC↔β-CMC
  (ρ=0.55) es la evidencia robusta del componente motor.
- El aperiódico no se midió en Paper 1 (test (e) lo cerraría).
- Ni offline ni online hay enlace individual theta↔conducta; el beneficio es a nivel de condición/estado.

Refs ancla: Espinoza2025 · Baddeley1986 · Omejc2025 · Gao2017 · Donoghue2020 · Miskovic2024 · Zelano2016/
Tort2018 (ritmo periférico→córtex) · Maurer2015 (fm-theta×desempeño) · Klimesch2018 (brain-body).
Fuentes: Paper 1 PMC12698685 (nature.com/articles/s41598-025-27606-5).
