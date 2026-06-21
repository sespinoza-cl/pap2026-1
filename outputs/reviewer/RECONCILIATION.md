# RECONCILIATION (F8) — histórico vs recomputado canónico
Fecha: 2026-06-19 · θ=4–7 · ROI-18 · N=31 chewing + 15 non-chewing · seed=42
Fuentes históricas: BITACORA.md "NÚMEROS CLAVE", PAC/FOOOF_results_summary.txt (STALE), enunciado misión.
MATCH = igual (±redondeo) · STALE = cambió (recomputar) · CONTRADICHO = el hallazgo no se sostiene.

## Diseño
| Ítem | Histórico | Canónico | Estado |
|---|---|---|---|
| Grupos | "31 casos DC/TMD + 15 controles" | 31 chewing + 15 non-chewing, **ambos sanos** | **CONTRADICHO** (no hay TMD) |
| Manipulation check | — | Casos Ch−Nc=+6.1 dB 31/31; Ctrl ≈0 5/15; p=7e-8 | NUEVO |

## Conducta
| Métrica | Histórico | Canónico | Estado |
|---|---|---|---|
| RT basal Casos | 618.2 ms | 591.1±112 ms | STALE (menor; ¿mediana vs media?) |
| RT basal Controles | 578.4 ms | 584.2±98 ms | MATCH |
| Equivalencia basal | p_equiv=0.981 | LME Grupo p=0.86 (sin dif.) | MATCH (concepto) |
| LME RT Grupo×Bloque | F=3.81 p=0.054 | p=0.054 | **MATCH** |
| ΔRT Casos (B1→B2) | −74.5 ms | −72.4 ms | MATCH |
| ΔRT Controles | −30.5 ms | −13.5 ms | STALE |
| Mejora RT casos (test) | p=0.14 (t pareado) | **Wilcoxon p=0.0001** | **CONTRADICHO** (test equivocado; el efecto SÍ es sig.) |
| Mejora IES casos | p=0.66 (t) / d=−0.49 | **Wilcoxon p=0.0008** | **CONTRADICHO** (idem) |

## TF (CBPT)
| Métrica | Histórico | Canónico | Estado |
|---|---|---|---|
| Cluster Casos Ch>Nc | stat=15989 p=0.010 | stat=15988.9 **p=0.005** | MATCH (stat); STALE (p) |
| Cluster Controles | p=0.931 | p=1.0000 (n.s.) | MATCH |
| Interacción (2D fullband) | p=0.0002 (3D) | **p=0.010** (S2c 2D) | STALE (análisis distinto; canónico 0.010) |
| Topo Casos Ch>Nc (FDR) | Fpz, AFz | Fpz, AFz | **MATCH** |
| CBPT temporal por banda | — | θ p=0.037, α p=0.0036, β p=0.049 | NUEVO (especificidad NO en CBPT temporal) |
| ROI canónico | ambiguo (3/8/11/18; roi_cbpt=3) | **ROI-18** (S2d reproduce) | STALE→FIJADO (I4) |

## FOOOF / specparam
| Métrica | Histórico | Canónico | Estado |
|---|---|---|---|
| Exponente χ Casos Nc→Ch | "slope −8.65" (fit log-lineal) | **χ 1.110→0.735** (specparam) | CONTRADICHO (métrica vieja superseded) |
| χ Casos Ch vs Nc | LME F=26.31 p<0.001 | Wilcoxon p<0.0001 | MATCH (concepto) |
| χ Controles | — | 1.065→1.030 p=0.39 (n.s.) | MATCH (especif. grupo) |
| % pico theta Casos-Ch | 94% | 97% | STALE (menor) |

## PAC
| Métrica | Histórico | Canónico | Estado |
|---|---|---|---|
| Sujetos zMI>1.96 (theta) | **15/31** | **22/31** circ · 21/31 AAFT · 21/31 ambos | **CONTRADICHO/STALE** (C3) |
| Banda beta | 13–20 Hz | 13–30 Hz | STALE (C1) |
| Supresión beta Late | p=0.021 | no se reproduce (β B→L p=0.49) | **CONTRADICHO** (C1) |
| Rayleigh theta continuo | Z=13.24 R=0.653 | Z=13.47 R=0.659 | MATCH |
| MI θ Late × medRT | ρ=−0.377 p=0.037 | ρ=−0.22 p=0.23 (n.s.) | **CONTRADICHO** |
| MI θ × IES | ρ=−0.019 (n.s.) | ρ=−0.08 (n.s.) | MATCH |
| ROI usado para PAC | ROI_FC (4 elec) | ROI-18 | STALE (I1) |
| Null | solo circ-shift | **doble: circ + AAFT** | MEJORADO (R7) |
| zMI_nc | null de Ch (sesgo) | null propio de Nc | CORREGIDO (I3) |
| Proxy muscular × Δtheta | r=0.212 p=0.093 | EMG masetero RMS ρ=−0.22 n.s. | reformulado (A3 limpio) |
| PAC Ch vs Nc (condición) | — | Ch M=13.32 vs Nc 1.08, p=0.0001 | NUEVO (A5) |
| PAC control negativo (controles) | — | zMI θ M=0.52, 1/15 (n.s.) | NUEVO (S8) |
| Dosis-respuesta (engagement×efectos) | "dose-response" (rechazado R1.M4) | todas n.s. (S7) | CONFIRMADO ausente |

## Resumen
- **MATCH/robustos:** topo Fpz+AFz, Rayleigh theta, LME RT Grupo×Bloque, equivalencia basal, χ Ch<Nc.
- **STALE corregidos:** 15/31→22/31, β13-20→13-30, ROI→18, %peak 94→97, cluster p 0.010→0.005.
- **CONTRADICHOS (no se sostienen):** supresión beta Late, MI×medRT, "slope −8.65", marco TMD,
  y los tests t de conducta (el Wilcoxon SÍ da sig → el efecto conductual era subreportado).
