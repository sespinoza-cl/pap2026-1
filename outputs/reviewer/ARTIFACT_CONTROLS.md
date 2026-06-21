# ARTIFACT_CONTROLS — ¿El efecto theta/PAC es artefacto motor? (A1–A6)
Fecha: 2026-06-19 · θ=4–7 · ROI-18 · N=31 casos · seed=42
Fuentes: S6_artifact_controls.mat, S6b_A3_muscle.mat, v1_S4b_PAC_ROI.mat, FOOOF (F4), S4h_comodulogram_stats.mat

## VEREDICTO GLOBAL
El confound de artefacto motor **NO se sostiene como motor del efecto**. Convergen 5 controles
independientes; queda 1 caveat (la magnitud del PAC es broadband, no theta-específica → la
especificidad theta se sostiene en topografía + condición + FOOOF, no en magnitud del PAC).

---

## A1 — Disociación espacial ✅
Δtheta (Ch−Nc) por electrodo, t-test + FDR (de theta_topo Morlet).
- **Frontal-medial:** Fpz t=+4.22 pFDR=**0.013***, AFz t=+3.58 pFDR=**0.038***; Fz/AF3/AF4/F1/F2
  todos t>0 (p<0.03 sin sobrevivir FDR). Frontal-medial sig (FDR)=2/8, todos positivos.
- **Temporal (proxy muscular T7/T8/FT7/FT8/TP7/TP8):** **0/6 sig** (t=1.3–2.2, pFDR 0.14–0.28).
→ El efecto es frontal-medial; si fuera masetero/temporal, los temporales serían los máximos. No lo son.

## A2 — Disociación espectral (FOOOF/specparam) ✅
- Exponente aperiódico χ: Casos Nc=1.110 → Ch=**0.735** (Wilcoxon p<0.0001); Controles 1.065→1.030 (p=0.39 n.s.).
- Pico theta genuino sobre el aperiódico: Casos-Chew **97%** (CF=6.58 Hz, central en 4–7).
→ Hay una oscilación theta real (gaussiana sobre 1/f), no un artefacto broadband. El aplanamiento
  aperiódico es específico de casos (cambio de excitabilidad neuromodulatorio, no muscular).

## A3 — Proxy muscular ✅ (con el test correcto)
El primer intento (30–40 Hz en ROI frontal) dio ρ=+0.45 p=0.013, pero ese proxy mezcla gamma neural
con frontalis → ambiguo. Con proxies LIMPIOS (S6b):
- **(a) EMG masetero real (RMS 20–200 Hz, chew) × Δtheta: ρ=−0.221 p=0.232 (n.s.)** ← test directo.
- (b) Δ(30–40 Hz) TEMPORAL × Δtheta: ρ=+0.374 p=0.039 (sig)
- (ref) Δ(30–40 Hz) FRONTAL × Δtheta: ρ=+0.446 p=0.013 (sig)
→ **La cantidad real de actividad masetera NO predice el theta** (ρ negativo, n.s.). El 30–40 Hz EEG sí
  correlaciona, pero **igual de fuerte en frontal que en temporal** (no focal a sitios musculares) y
  **desacoplado del masetero real** → refleja el aplanamiento aperiódico compartido (A2), un cambio de
  excitabilidad neural, no spillover del músculo.

## A4 — No-sinusoidalidad (clave para C2) ✅
- **(i) Doble null:** el PAC theta supera AMBOS nulls — circular-shift (22/31) y **AAFT (21/31)**.
  Como AAFT preserva el espectro y rompe la forma de onda, superar AAFT descarta el artefacto de
  transiente no-sinusoidal (Kramer 2008; Aru 2015). zMI θ vs AAFT: M=12.85, Wilcoxon p=1.0×10⁻⁵.
- **(ii) Sharpness:** razón armónica de la envolvente EMG (2f/f, M=0.293±0.105) × zMI theta:
  **ρ=+0.165 p=0.375 (n.s.)** → el MI **no escala** con la no-sinusoidalidad de la onda EMG.
- **(iii) Comodulograma DESCRIPTIVO (S4h):** cluster 2D significativo (p<0.0001) pero **broadband**
  fa=[4–30 Hz] (θ+α+β), fp=[0.4–1.8 Hz]. Se reporta como descriptivo, no como prueba de PAC θ-exclusivo.

## A5 — Especificidad por condición ✅
zMI theta: Chew M=13.32 vs **NoChew M=1.08** (Wilcoxon p=**0.0001**); NoChew solo 7/31 sig.
→ El PAC está impulsado por la masticación; casi ausente sin ritmo masticatorio (la fase en f_chew±0.5
  en NoChew es ruido). Pilar mecanístico y control de condición.

## A6 — Anclaje conductual ⚪ (honesto: n.s.)
MI θ Late × RT ρ=−0.22 p=0.23; × IES ρ=−0.08 p=0.69; zMI × RT/IES n.s.
→ No hay correlación PAC–conducta significativa en esta muestra. El ρ=−0.377 histórico NO se reproduce.

---

## SÍNTESIS PARA EL REVISOR
El aumento de theta frontal durante la masticación y su acoplamiento a la fase masticatoria **no se
explican por artefacto motor**: (1) la topografía es frontal-medial, no temporal (A1); (2) la cantidad
real de EMG masetero no predice el theta (A3a); (3) el acoplamiento sobrevive surrogados AAFT y no
escala con la forma de onda EMG (A4); (4) existe un pico theta gaussiano genuino sobre el aperiódico
(A2); (5) el efecto es específico de la condición de masticación (A5). La covariación Δtheta–Δ(30–40 Hz)
refleja el aplanamiento aperiódico neural (A2), no spillover muscular (desacoplada del masetero real).

**Límites declarados (R5):** (a) el PAC es broadband (β≫α≫θ en magnitud) — la especificidad theta se
sostiene en la topografía de interacción + condición + FOOOF, NO en la magnitud del PAC; (b) sin
anclaje conductual (A6).

**Confound de orden/práctica (I6) — bien controlado por diseño.** Ambos grupos son sanos y comparten el
mismo orden fijo (B1 NoChew → B2); la ÚNICA diferencia es que el grupo "Casos" masticó en B2 y el grupo
"Controles" no (repitió la tarea sin masticar = control de práctica). La interacción Grupo×Bloque aísla
el efecto de masticación controlando práctica/repetición. **Manipulation check** (masetero 1–2.5 Hz,
S0c_chew_engagement.m): Casos Ch−Nc=+6.10 dB (31/31, p=1.2×10⁻⁶); Controles −0.72 dB (5/15, n.s.);
entre grupos p=7×10⁻⁸. Caveat: no hay brazo de sanos-masticando independiente, por lo que el contraste
mezcla masticación con asignación de grupo (ambos sanos); el control de práctica sí queda cubierto.
