# CLAIMS_SUMMARY — Qué sostiene el paper (canónico, 2026-06-19)
Diseño: 46 adultos SANOS. Manipulación entre-sujetos de masticación en el 2º bloque:
"Casos" (n=31) masticaron; "Controles" (n=15) repitieron la tarea sin masticar (control de práctica).
Ambos: B1=NoChew baseline, orden fijo. La interacción Grupo×Bloque aísla el efecto de masticar
controlando práctica. Manipulation check (masetero): Casos Ch−Nc=+6.1 dB (31/31); Controles ≈0 (n.s.).

Niveles: **SÓLIDO** (sobrevive controles/recompute) · **MODERADO** (sig. pero con matiz) · **TENDENCIA** · **NULO** (honesto).

---

## 1. CONDUCTA
**Claim 1.1 (SÓLIDO).** Masticar durante la tarea mejora el desempeño en quienes mastican.
- RT casos: 591→519 ms (Wilcoxon p=0.0001, d=−0.58); IES casos: 994→812 (p=0.0008, d=−0.49).
- Controles (práctica, sin masticar) NO mejoran: ΔRT=−13.5, ΔIES=+45.
- Sin diferencia basal entre grupos (LME Grupo p=0.86).

**Claim 1.2 (TENDENCIA).** El beneficio es mayor en los que mastican que en el control de práctica.
- LME RT Grupo×Bloque p=0.054; IES p=0.106. Es tendencia, no significativo.
- Honesto: el efecto dentro de casos es fuerte; el contraste diferencial vs práctica es solo tendencia.

**Caveat:** orden fijo NoChew→Chew, pero **controlado** por el grupo no-masticador (no mejora → la
mejora de casos no es práctica). No hay brazo de sanos-masticando independiente.

## 2. TIEMPO-FRECUENCIA (TF)
**Claim 2.1 (SÓLIDO).** Masticar aumenta la potencia cortical Chew>NoChew en los que mastican, no en controles.
- CBPT 2D fullband: Casos cluster p=0.005; Controles p=1.0 (n.s.); Interacción p=0.010.

**Claim 2.2 (SÓLIDO).** El efecto de interacción es theta frontal-medial.
- Topografía interacción FDR → ROI-18 frontal; Casos Ch>Nc FDR: **Fpz, AFz** (responde R1.M6).
- Disociación espacial (A1): frontal-medial sig, temporal 0/6 → no muscular.

**Claim 2.3 (MATIZ — importante).** La "especificidad theta" vive en la TOPOGRAFÍA de interacción,
NO en el CBPT temporal: en el eje tiempo, theta (p=0.037), alpha (p=0.0036) y beta (p=0.049) salen
significativas en casos. Decir "theta-específico" solo respecto a la topografía interacción.

## 3. FOOOF / specparam
**Claim 3.1 (SÓLIDO).** Hay una oscilación theta genuina, no un artefacto broadband.
- Pico theta sobre el aperiódico en 97% de casos-chew (CF=6.58 Hz, central en 4–7) → responde R1.M3.

**Claim 3.2 (SÓLIDO).** Masticar aplana el componente aperiódico (1/f) solo en los que mastican.
- Exponente χ casos: 1.110→0.735 (Wilcoxon p<0.0001); controles 1.065→1.030 (n.s.).
- Mecanismo: cambio de excitabilidad neural (no muscular). Explica la covariación Δtheta–Δ(30–40Hz) (A3).

## 4. PAC (acoplamiento fase masticatoria–amplitud cortical)
**Claim 4.1 (SÓLIDO).** Existe acoplamiento theta genuino a la fase masticatoria, sobre doble null.
- zMI theta supera circular-shift (22/31) Y AAFT (21/31); ambos 21/31; binom p≈0; Wilcoxon p≈1e-5.
- AAFT preserva espectro y rompe forma de onda → descarta artefacto no-sinusoidal (R7/A4i).

**Claim 4.2 (SÓLIDO).** El PAC está impulsado por la masticación (especificidad de condición).
- zMI theta Ch M=13.32 vs NoChew M=1.08 (Wilcoxon p=0.0001); NoChew solo 7/31 (A5).

**Claim 4.3 (MATIZ honesto).** El PAC NO es theta-específico en magnitud: es broadband.
- zMI: beta (88–120) > alpha (22) > theta (13); Friedman p<0.0001. Comodulograma 2D significativo pero
  broadband fa=[4–30 Hz] (S4h). Por eso se reporta como descriptivo; la especificidad theta se sostiene
  en topografía + condición + FOOOF, NO en magnitud del PAC.

**Claim 4.4 (NULO, honesto).** Sin anclaje conductual ni dosis-respuesta.
- PAC × RT/IES n.s. (A6). Intensidad de masticación × theta/PAC/mejora: todas n.s. (S7); el único
  nominal (f_chew×MI_Late ρ=−0.47 p=0.008) no sobrevive FDR. El efecto es categórico (Ch vs Nc), no graduado.

## 5. CONTROL NEGATIVO (PAC en controles, S8) — SÓLIDO
**Claim 5.1.** El PAC sigue a la masticación REAL, no a la etiqueta del bloque ni al grupo.
- Controles (no mastican) en su bloque "Ch": zMI theta M=0.52, Mdn=0.20, **1/15 sig**, Wilcoxon p=0.094 (n.s.).
- Completa el 2×2: Casos-Chew M_z=13.32 (22/31) ✅ ≫ Casos-NoChew 1.08 (7/31) ≈ Controles-"Chew" 0.52 (1/15).
- → El acoplamiento aparece sólo cuando hay ritmo masticatorio genuino; descarta que sea un efecto del
  diseño/bloque/grupo. Es el control negativo más directo del mecanismo.

---

## ESTRUCTURA NARRATIVA SUGERIDA PARA LA RESPUESTA AL REVISOR
1. Masticar mejora el desempeño (conducta sólida) y el control de práctica no mejora → no es práctica.
2. Masticar aumenta theta frontal-medial (TF + topografía), oscilación genuina sobre 1/f aplanado (FOOOF).
3. La amplitud theta se acopla a la fase masticatoria por encima de doble null (PAC), específico de la
   condición de masticación, y robusto a controles de artefacto motor (A1–A5).
4. Límites declarados: PAC broadband (especificidad por topografía/condición), sin anclaje conductual ni
   dosis-respuesta, sin brazo de sanos-masticando (no se aísla especificidad de grupo aparte de masticación).
