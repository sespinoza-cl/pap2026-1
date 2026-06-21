# LIT_SUPPORT — Soporte de literatura + auditoría para el marco DUAL (θ-cognitivo + β-CMC)
Fecha: 2026-06-20 · Paper2 V1 (masticación rítmica × Two-Back × EEG/EMG, adultos SANOS, sin TMD)
Marco elegido por el usuario: **DUAL** — (A) θ-PAC cognitivo + (B) β-PAC como acoplamiento corticomuscular.
Fuentes: `opencite` (Semantic Scholar / Crossref / OpenAlex / PubMed), sólo peer-reviewed (excluidos
arXiv/bioRxiv/medRxiv/Zenodo/Research Square). Auditoría cruzada independiente: Gemini (bridge 001/002).

> Anclas empíricas canónicas que NO se pueden contradecir (RESULTS_CANONICAL.md):
> PAC broadband con β el más fuerte (zMI β 88–120 > α 22 > θ 13); mejora conductual dentro de casos
> SÓLIDA pero contraste diferencial Grupo×Bloque sólo TENDENCIA (p=0.054); PAC×conducta n.s.

---

## 0. VEREDICTO DE LA AUDITORÍA (Gemini, independiente) + RECONCILIACIÓN

Gemini auditó CLAIMS_SUMMARY + RESULTS_CANONICAL. Coincidimos en lo sustantivo; hay **dos flags HIGH**
que son guardarraíles de redacción, no errores de análisis. Resumen y mi reconciliación:

| # | Flag de Gemini | Sev. | ¿De acuerdo? | Acción |
|---|---|---|---|---|
| 1 | "El beneficio es mayor en quienes mastican" (Claim 1.2) es **stat_dubious**: p=0.054 + confound de orden | media | **SÍ** (ya es caveat canónico) | Redactar como "sugiere un beneficio adicional, aunque el contraste directo no alcanzó significancia (p=0.054)". |
| 2 | "Especificidad theta" es **internal_contradiction** si α/β también suben | **alta** | **Parcial** | El marco DUAL ya **abandona** la "exclusividad theta". Usar "incremento broadband con topografía frontal-medial theta-dominante en la interacción". La especificidad theta se reserva a TOPOGRAFÍA+CONDICIÓN+FOOOF, nunca a magnitud PAC. |
| 3 | "cambio de excitabilidad neural (E/I)" es **overgeneralization** (proxy inferido) | media | **SÍ** | Moderar: "interpretable como un posible desplazamiento del balance E/I (proxy del exponente 1/f; Gao 2017)". |
| 4 | **"β-PAC = CMC" es unsupported**: CMC clásico es coherencia fase-FASE; nuestro β-PAC es fase-AMPLITUD. Equipararlos sin justificación es overreach | **alta** | **SÍ — guardarraíl crítico** | NO escribir "β-PAC = CMC". Escribir "β-PAC como **firma de acoplamiento corticomuscular** sincronizada a la fase motora, **análoga/funcionalmente relacionada** al CMC", y justificarlo con (i) literatura emergente de *corticomuscular PAC* y (ii) el precedente de ritmos periféricos que acoplan al córtex (respiración). |

**ELEVAR AL USUARIO (decisión de redacción, no de análisis):** el flag #4. El marco dual es válido,
pero el rótulo debe ser **"β-PAC análogo a CMC"** y no una igualdad. La literatura abajo (§B) da la
cobertura para sostenerlo con honestidad.

Gemini lit-gap (002) señaló 3 vacíos —todos cubiertos abajo—: (i) *corticomuscular PAC* reciente como
extensión válida del CMC clásico; (ii) eje trigémino→LC→NE con evidencia pupilométrica; (iii) modulación
**dinámica** (no sólo rasgo/edad) del aperiódico durante WM. *Limitación:* Gemini devolvió punteros
temáticos, no DOIs duros; los DOIs los aporté vía opencite.

---

## A. θ-PAC COGNITIVO — theta frontal-medial, WM, pico genuino sobre 1/f, aperiódico=E/I

**A1. El theta frontal-medial es el sustrato canónico de la carga de WM y el control cognitivo.**
La potencia theta de línea media frontal escala con la carga de WM ([Onton 2005](https://doi.org/10.1016/j.neuroimage.2005.04.014);
[Gevins 2000](https://doi.org/10.1093/cercor/10.9.829)) y opera como mecanismo de control cognitivo cuya
frecuencia se desplaza hacia un óptimo según demanda ([Senoussi 2022](https://doi.org/10.1038/s41562-022-01335-5)).
Revisiones recientes consolidan el vínculo theta↔función cognitiva ([Tan 2024](https://doi.org/10.1016/j.dcn.2024.101404))
y el theta medio-frontal como índice transdiagnóstico de control ([McLoughlin 2021](https://doi.org/10.1016/j.biopsych.2021.08.020)).
Sostiene Claim 2.2 (topografía frontal-medial) y la lectura WM del θ-PAC.

**A2. El exponente aperiódico (1/f) es un proxy de excitabilidad E/I que se modula DINÁMICAMENTE en tareas.**
El aplanamiento del exponente indexa un desplazamiento hacia excitación ([Gao 2017]; seed), y su separación
del componente oscilatorio es metodológicamente necesaria ([Thuwal 2021](https://doi.org/10.1523/eneuro.0224-21.2021)).
La excitabilidad talamo-cortical guía el desempeño perceptual ([Kosciessa 2021](https://doi.org/10.1038/s41467-021-22511-7))
y la potencia aperiódica discrimina estados E/I ([Wilkinson 2021](https://doi.org/10.1186/s13229-021-00425-x)).
→ Cobertura del vacío #3 de Gemini: el cambio χ 1.110→0.735 (Claim 3.2) es interpretable como modulación
**aguda** del balance E/I por la masticación, no un rasgo. **Redacción moderada (flag #3).**

## B. β-PAC COMO ACOPLAMIENTO CORTICOMUSCULAR (no "=CMC") — guardarraíl + cobertura

**B1. La banda β es la firma canónica del drive cortico→muscular (CMC).** El CMC β refleja interacción
corticospinal efectiva durante salida motora sostenida ([Kristeva 2007]; seed) y la sincronía dentro del
sistema motor es banda-β por defecto ([van Wijk 2012](https://doi.org/10.3389/fnhum.2012.00252);
[Boonstra 2009](https://doi.org/10.1016/j.neulet.2009.07.043)). La potencia β cortical se modula con la
planificación/ejecución motora ([Zaepffel 2013](https://doi.org/10.1371/journal.pone.0060060);
[Solis-Escalante 2018](https://doi.org/10.1016/j.neuroimage.2018.12.045)).

**B2. El CMC es específicamente relevante a la musculatura masticatoria/oromotora.** Existe acoplamiento
EEG-EMG durante actividad de músculos mandibulares y de cuello en masticación ([Ishii 2016](https://doi.org/10.1016/j.physbeh.2016.03.023)),
control cerebral de la deglución ([Cheng 2022](https://doi.org/10.1016/j.jns.2022.120434)), y el CMC se
modula finamente con la demanda de la tarea motora ([Glories 2021](https://doi.org/10.1038/s41598-021-85851-w);
[Glories 2023](https://doi.org/10.1111/sms.14309); [Forman 2022](https://doi.org/10.1007/s00421-022-04938-y);
[Yang 2010](https://doi.org/10.1109/tnsre.2010.2047173)). → Sostiene que un β prominente acoplado a la fase
masticatoria es **esperable y fisiológico**, no un artefacto incómodo (Claim 4.3).

**B3. GUARDARRAÍL (flag #4): β-PAC ≠ CMC en método.** El CMC clásico es coherencia **fase-fase** EEG-EMG;
nuestro hallazgo es PAC **fase(periférica)-amplitud(cortical)**. Escribir "β-PAC = CMC" es un overreach.
*Cobertura honesta:* (i) la arquitectura de oscilaciones **cerebro-cuerpo** acopla ritmos periféricos y
corticales de forma jerárquica ([Klimesch 2018](https://doi.org/10.1111/ejn.14192)); (ii) un ritmo motor
periférico que **entrena/acopla** amplitud cortical es fenómeno establecido (ver §D respiración). Por tanto
el β-PAC se describe como "**firma de acoplamiento corticomuscular análoga al CMC**", reconociendo la
distinción metodológica. (El vacío #1 de Gemini —*corticomuscular PAC* reciente— apunta aquí; conviene una
búsqueda dirigida adicional de "corticomuscular phase-amplitude coupling" para 1–2 citas duras antes de la discusión.)

## C. METODOLOGÍA — PAC a un ritmo periférico + controles de no-sinusoidalidad

**C1. El PAC puede ser espurio por morfología no-sinusoidal; nuestros controles (AAFT + sharpness) son los
correctos.** Oscilaciones de morfología no-sinusoidal generan CFC espuria ([Lozano-Soldevilla 2016](https://doi.org/10.3389/fncom.2016.00087);
[Gerber 2016](https://doi.org/10.1371/journal.pone.0167351)); la forma de onda debe controlarse
([Cole/Jackson 2019](https://doi.org/10.1523/eneuro.0151-19.2019)). Hay criterios para discriminar PAC válido
de espurio ([Jensen 2016](https://doi.org/10.1523/eneuro.0334-16.2016); [Seymour 2017](https://doi.org/10.3389/fnins.2017.00487))
y métodos robustos ([Dupré la Tour 2017](https://doi.org/10.1371/journal.pcbi.1005893);
[Munia 2019](https://doi.org/10.1038/s41598-019-48870-2); [Siebenhühner 2020](https://doi.org/10.1371/journal.pbio.3000685)).
→ Sostiene que el doble null (circular-shift + **AAFT**, que rompe la forma de onda) es la defensa adecuada
(Claim 4.1) y que reportar el comodulograma como "descriptivo" es lo honesto (Claim 4.3).

**C2. Acoplar amplitud cortical a la fase de un ritmo MOTOR cortical es precedente directo.** La actividad
motora cortical humana está fase-acoplada de forma selectiva ([Miller 2012](https://doi.org/10.1371/journal.pcbi.1002655)).

## D. PRECEDENTE FUERTE — ritmos periféricos (respiración) que acoplan/entrenan el córtex y la cognición

Este es el respaldo más potente para "PAC a un ritmo periférico" como fisiología establecida (no artefacto):
la **respiración nasal entrena oscilaciones límbicas y modula la cognición/memoria** ([Zelano 2016](https://doi.org/10.1523/jneurosci.2586-16.2016));
los ritmos cerebrales entrenados por respiración son **globales** ([Tort 2018](https://doi.org/10.1016/j.tins.2018.01.007)
— nota: Tort es el autor del Modulation Index que usamos); la respiración como **ritmo fundamental** de la
función cerebral ([Heck 2017](https://doi.org/10.3389/fncir.2016.00115)); y los **ritmos interoceptivos**
estructuran la actividad cortical ([Engelen 2023](https://doi.org/10.1038/s41593-023-01425-1)).
→ Analogía directa: la masticación es otro ritmo oromotor periódico cuya fase puede acoplar amplitud cortical.
**Refuerza el marco entero** y desactiva el reflejo "PAC a periferia = artefacto".

## E. EJE TRIGÉMINO → LC-NE → EXCITABILIDAD CORTICAL (driver de arousal por masticación)

**E1. La masticación mejora la cognición vía aferencias trigeminales.** Entradas trigéminas/visceral/vestibular
mejoran la cognición ([De Cicco 2018](https://doi.org/10.3389/fnana.2017.00130); [De Cicco 2016](https://doi.org/10.1371/journal.pone.0148715));
masticar mantiene función hipocampal ([Chen 2015](https://doi.org/10.7150/ijms.11911)) y la deficiencia
masticatoria es factor de riesgo cognitivo ([Teixeira 2014](https://doi.org/10.7150/ijms.6801)). Masticar
mejora atención sostenida ([Hirano 2015](https://doi.org/10.1155/2015/367026); revisión [Allen 2015](https://doi.org/10.1155/2015/654806))
y es conducta de afrontamiento al estrés ([Kubo 2015](https://doi.org/10.1155/2015/876409)). Evidencia
reciente busca el enlace neural masticatorio-cognitivo directo ([Kang 2024](https://doi.org/10.3389/fncel.2024.1425645)).
→ Mecanismo blando para la cascada masticación→arousal→excitabilidad (aplanamiento 1/f). **Mantener como
'consistente con', no 'demuestra'** (estilo Wael; flag #3). *Vacío #2 de Gemini (pupilometría LC):* conviene
1 cita dura trigémino-LC-pupila si se quiere endurecer; pendiente de búsqueda dirigida.

---

## REFERENCIAS NUEVAS A INTEGRAR (peer-reviewed, con DOI) — staging para biblio.bib

> Se depositan en `outputs/reviewer/biblio_new.bib` (Task #2). NO se vuelcan aún a `biblio.bib` canónico:
> se mergean en Fase D sólo las que efectivamente cite la discusión, para no contaminar el .bib.

**β-CMC / oromotor (§B):** Glories 2021 `10.1038/s41598-021-85851-w` · Glories 2023 `10.1111/sms.14309` ·
Ishii 2016 `10.1016/j.physbeh.2016.03.023` · Boonstra 2009 `10.1016/j.neulet.2009.07.043` ·
van Wijk 2012 `10.3389/fnhum.2012.00252` · Yang 2010 `10.1109/tnsre.2010.2047173` ·
Zaepffel 2013 `10.1371/journal.pone.0060060` · Solis-Escalante 2018 `10.1016/j.neuroimage.2018.12.045` ·
Forman 2022 `10.1007/s00421-022-04938-y` · Cheng 2022 `10.1016/j.jns.2022.120434`

**Ritmo periférico → córtex (§D):** Zelano 2016 `10.1523/jneurosci.2586-16.2016` ·
Tort 2018 `10.1016/j.tins.2018.01.007` · Heck 2017 `10.3389/fncir.2016.00115` ·
Engelen 2023 `10.1038/s41593-023-01425-1` · Klimesch 2018 `10.1111/ejn.14192`

**θ-cognitivo / aperiódico-E/I (§A):** Onton 2005 `10.1016/j.neuroimage.2005.04.014` ·
Gevins 2000 `10.1093/cercor/10.9.829` · Senoussi 2022 `10.1038/s41562-022-01335-5` ·
Tan 2024 `10.1016/j.dcn.2024.101404` · McLoughlin 2021 `10.1016/j.biopsych.2021.08.020` ·
Thuwal 2021 `10.1523/eneuro.0224-21.2021` · Kosciessa 2021 `10.1038/s41467-021-22511-7` ·
Wilkinson 2021 `10.1186/s13229-021-00425-x`

**PAC validez / no-sinusoidalidad (§C):** Lozano-Soldevilla 2016 `10.3389/fncom.2016.00087` ·
Gerber 2016 `10.1371/journal.pone.0167351` · Jensen 2016 `10.1523/eneuro.0334-16.2016` ·
Cole/Jackson 2019 `10.1523/eneuro.0151-19.2019` · Siebenhühner 2020 `10.1371/journal.pbio.3000685` ·
Dupré la Tour 2017 `10.1371/journal.pcbi.1005893` · Seymour 2017 `10.3389/fnins.2017.00487` ·
Munia 2019 `10.1038/s41598-019-48870-2` · Miller 2012 `10.1371/journal.pcbi.1002655`

**Masticación-cognición / trigémino-LC (§E):** Chen 2015 `10.7150/ijms.11911` ·
Teixeira 2014 `10.7150/ijms.6801` · De Cicco 2018 `10.3389/fnana.2017.00130` ·
De Cicco 2016 `10.1371/journal.pone.0148715` · Hirano 2015 `10.1155/2015/367026` ·
Allen 2015 `10.1155/2015/654806` · Kubo 2015 `10.1155/2015/876409` · Kang 2024 `10.3389/fncel.2024.1425645`

---

## GAPS / PENDIENTES (búsquedas dirigidas antes de cerrar la discusión)
1. **"corticomuscular phase-amplitude coupling"** — 1–2 citas duras que legitimen PAC-a-periferia como
   extensión del CMC (cobertura directa del flag #4 y vacío #1 de Gemini).
2. **trigémino-LC-NE con pupilometría durante masticación** — 1 cita dura para endurecer §E (vacío #2).
3. Confirmar que las refs semilla de v2 (Conway1995, Omejc2025, Voytek2015, Gao2017, Donoghue2020,
   Miskovic2024, Jensen2010, Klimesch2012, Karlsson2023, Weisz2020, BuzsakiMoser2013, baddeley1986,
   espinoza2025) ya estén con DOI en biblio.bib/biblio2.bib (no re-resolver si ya existen).

## CÓMO ANCLA TÍTULO/DISCUSIÓN (Fases D/E)
- **Título:** evitar "theta-PAC específico"; abrazar lo dual ("masticación rítmica → theta frontal-medial +
  acoplamiento cortico-fase-masticatoria + mejora de WM en adultos sanos").
- **Discusión:** esqueleto de v2 (β-CMC + aperiódico-E/I + α-ERD + dual-stream + cascada temporal) PERO con
  los guardarraíles #2/#3/#4 y anclado a números canónicos. β-PAC = "análogo a CMC", no "=CMC".
