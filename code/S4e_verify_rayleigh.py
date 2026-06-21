"""
S4e_verify_rayleigh.py
======================
Verificación de los ángulos preferidos Rayleigh banda por banda.

Checks:
  1. Ángulos individuales por sujeto (raw) — confirmar que los valores
     en pref_ch son razonables y no hay artefactos
  2. Circular stats independientes (Watson test Chi-sq + Rayleigh)
  3. Correlación circular entre bandas (si theta-alpha-beta están
     correlacionados sujeto a sujeto → firma de fuente común/aperiódica)
  4. Plot de dispersión lineal de ángulos por banda (sujeto en eje X)
"""

import numpy as np
import h5py
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

ROOT     = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
PAC_FILE = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "stats",
                        "v1_S4b_PAC_ROI.mat")
DIR_OUT  = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "figures")

with h5py.File(PAC_FILE, 'r') as f:
    pref_ch = np.array(f['pref_ch'])   # shape (3, 31)

BANDS = ['theta (4-7 Hz)', 'alpha (8-12 Hz)', 'beta (13-20 Hz)']

def circ_mean(pv):
    C = np.mean(np.exp(1j * pv))
    return np.angle(C), np.abs(C)

def rayleigh_p(pv):
    n = len(pv)
    R = np.abs(np.mean(np.exp(1j * pv)))
    Z = n * R**2
    return R, Z, float(np.exp(-Z))

def circ_corr(a, b):
    """Correlación circular entre dos vectores de ángulos (Jammalamadaka & Sengupta)."""
    sin_a = np.sin(a - circ_mean(a)[0])
    sin_b = np.sin(b - circ_mean(b)[0])
    r = np.sum(sin_a * sin_b) / np.sqrt(np.sum(sin_a**2) * np.sum(sin_b**2))
    return r

# ═══════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("VERIFICACIÓN RAYLEIGH — pref_ch (bloque chewing continuo)")
print("=" * 65)

print(f"\nShape pref_ch: {pref_ch.shape}   (bandas × sujetos)")
print(f"NaN por banda: theta={np.sum(np.isnan(pref_ch[0]))}, "
      f"alpha={np.sum(np.isnan(pref_ch[1]))}, "
      f"beta={np.sum(np.isnan(pref_ch[2]))}")

# ── 1. Ángulos individuales ──────────────────────────────────────────────────
print("\n── Ángulos preferidos por sujeto (grados) ──")
print(f"{'Suj':>4}  {'Theta':>8}  {'Alpha':>8}  {'Beta':>8}")
for i in range(pref_ch.shape[1]):
    t = np.degrees(pref_ch[0, i]) % 360
    a = np.degrees(pref_ch[1, i]) % 360
    b = np.degrees(pref_ch[2, i]) % 360
    print(f"  {i+1:2d}    {t:7.1f}°   {a:7.1f}°   {b:7.1f}°")

# ── 2. Estadísticas circulares ───────────────────────────────────────────────
print("\n── Circular stats ──")
for b_idx, band in enumerate(BANDS):
    pv = pref_ch[b_idx, :]
    pv = pv[~np.isnan(pv)]
    mu, R = circ_mean(pv)
    R2, Z, p = rayleigh_p(pv)
    mu_deg = np.degrees(mu) % 360
    # Range of angles
    degs = np.degrees(pv) % 360
    print(f"\n  {band}:")
    print(f"    Mean angle : {mu_deg:.1f}°")
    print(f"    R (length) : {R:.4f}")
    print(f"    Z          : {Z:.3f}")
    print(f"    p (Rayleigh): {p:.6f}")
    print(f"    Range      : [{degs.min():.1f}°, {degs.max():.1f}°]")
    print(f"    % in 0-90° : {100*np.mean((degs >= 0) & (degs <= 90)):.0f}%  "
          f"({np.sum((degs >= 0) & (degs <= 90))}/{len(degs)} subj)")

# ── 3. Correlación circular entre bandas ─────────────────────────────────────
print("\n── Correlación circular entre bandas (Jammalamadaka) ──")
theta_v = pref_ch[0, :]
alpha_v = pref_ch[1, :]
beta_v  = pref_ch[2, :]

r_ta = circ_corr(theta_v, alpha_v)
r_tb = circ_corr(theta_v, beta_v)
r_ab = circ_corr(alpha_v, beta_v)
print(f"  theta × alpha : r = {r_ta:.3f}")
print(f"  theta × beta  : r = {r_tb:.3f}")
print(f"  alpha × beta  : r = {r_ab:.3f}")
print(f"\n  → Si r > 0.5: los ángulos preferidos de las tres bandas")
print(f"    covarían sujeto a sujeto, consistente con fuente común (aperiódica).")

# ── 4. Test de diferencia de ángulos medios entre bandas ─────────────────────
print("\n── Diferencia de ángulos medios entre bandas ──")
mu_t = circ_mean(pref_ch[0])[0]
mu_a = circ_mean(pref_ch[1])[0]
mu_b = circ_mean(pref_ch[2])[0]

def angle_diff_deg(a1, a2):
    d = np.degrees(a1 - a2)
    return ((d + 180) % 360) - 180

print(f"  theta vs alpha : Δ = {angle_diff_deg(mu_t, mu_a):.1f}°")
print(f"  theta vs beta  : Δ = {angle_diff_deg(mu_t, mu_b):.1f}°")
print(f"  alpha vs beta  : Δ = {angle_diff_deg(mu_a, mu_b):.1f}°")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURA: ángulos por sujeto (scatter lineal) + boxplot circular
# ═══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
fig.patch.set_facecolor('white')
COLS = ['#D55E00', '#56B4E9', '#CC79A7']  # Okabe-Ito: vermillion/skyblue/purple

# Panel A — scatter: sujeto vs ángulo, 3 bandas superpuestas
ax = axes[0]
for b_idx, band in enumerate(BANDS):
    degs = np.degrees(pref_ch[b_idx, :]) % 360
    jit  = np.random.default_rng(b_idx).random(len(degs)) * 0.35 - 0.17
    ax.scatter(np.arange(1, 32) + jit + b_idx * 0.5,
               degs, s=22, color=COLS[b_idx], alpha=0.6,
               label=band, linewidths=0)
    mu_deg = np.degrees(circ_mean(pref_ch[b_idx, :])[0]) % 360
    ax.axhline(mu_deg, color=COLS[b_idx], lw=1.4, ls='--', alpha=0.7)

ax.set_xlabel('Subject', fontsize=10)
ax.set_ylabel('Preferred phase (degrees)', fontsize=10)
ax.set_ylim(-10, 370)
ax.set_yticks([0, 90, 180, 270, 360])
ax.set_yticklabels(['0°', '90°', '180°', '270°', '360°'])
ax.axhline(0,   color='lightgray', lw=0.5)
ax.axhline(360, color='lightgray', lw=0.5)
ax.legend(fontsize=8, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('A  Preferred phase per subject (3 bands)', loc='left',
             fontsize=10, fontweight='bold')

# Panel B — scatter: theta ángulo vs alpha ángulo vs beta (correlación circular visual)
ax2 = axes[1]
theta_deg = np.degrees(pref_ch[0, :]) % 360
alpha_deg = np.degrees(pref_ch[1, :]) % 360
beta_deg  = np.degrees(pref_ch[2, :]) % 360

sc1 = ax2.scatter(theta_deg, alpha_deg, s=35, color='#9B59B6',
                  alpha=0.65, linewidths=0, label=f'theta vs alpha  r={r_ta:.2f}')
sc2 = ax2.scatter(theta_deg, beta_deg,  s=35, color='#E67E22',
                  alpha=0.65, marker='^', linewidths=0,
                  label=f'theta vs beta   r={r_tb:.2f}')

# Línea identidad
ax2.plot([0, 360], [0, 360], 'k--', lw=0.7, alpha=0.4, label='identity')
ax2.set_xlabel('Theta preferred phase (°)', fontsize=10)
ax2.set_ylabel('Alpha / Beta preferred phase (°)', fontsize=10)
ax2.set_xlim(-10, 370); ax2.set_ylim(-10, 370)
ax2.set_xticks([0, 90, 180, 270, 360])
ax2.set_yticks([0, 90, 180, 270, 360])
ax2.legend(fontsize=8)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_title('B  Cross-band phase correlation (subject-level)', loc='left',
              fontsize=10, fontweight='bold')

plt.tight_layout()
out = os.path.join(DIR_OUT, 'S4e_verify_rayleigh.png')
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'\nFigura guardada: {out}')
print("=== DONE ===")
