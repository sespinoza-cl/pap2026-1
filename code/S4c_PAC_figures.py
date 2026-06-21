"""
S4c_PAC_figures.py
==================
Figuras PAC individuales sin texto para montaje en Inkscape.

Salidas (Analysis_V1_Final/outputs/figures/):
  S4c_PAC_zMI_theta_vs_null.png        — box+scatter zMI theta
  S4c_PAC_behavior_theta_4to7Hz.png    — scatter MI_late × IES (theta)
  S4c_PAC_behavior_alpha_8to12Hz.png   — scatter MI_late × IES (alpha)
  S4c_PAC_behavior_beta_13to20Hz.png   — scatter MI_late × IES (beta)
"""

import numpy as np
import scipy.io as sio
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy import stats
from scipy.stats import theilslopes
import os

ROOT     = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
PAC_FILE = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "stats",
                        "v1_S4b_PAC_ROI.mat")
BEH_FILE = os.path.join(ROOT, "Analysis_V1_Final", "data", "computed",
                        "v1_S1_behavior_stats.mat")
DIR_OUT  = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "figures")
os.makedirs(DIR_OUT, exist_ok=True)

with h5py.File(PAC_FILE, 'r') as f:
    zMI_ch     = np.array(f['zMI_ch'])      # (3, 31)
    MI_late_ch = np.array(f['MI_late_ch'])  # (31, 3)

beh    = sio.loadmat(BEH_FILE)
IES_b2 = beh['ies_cas_b2'].flatten()

BANDS_LABEL  = ['theta_4to7Hz', 'alpha_8to12Hz', 'beta_13to20Hz']
BANDS_XLABEL = ['MI  θ (4–7 Hz)  [×10⁻³]',
                'MI  α (8–12 Hz)  [×10⁻³]',
                'MI  β (13–20 Hz)  [×10⁻³]']
COLS         = ['#D55E00', '#56B4E9', '#CC79A7']  # Okabe-Ito: vermillion/skyblue/purple

# ═══════════════════════════════════════════════════════════════════════════
# FIG 1 — zMI theta vs null  (sin texto)
# ═══════════════════════════════════════════════════════════════════════════
z   = zMI_ch[0, :]
col = COLS[0]

fig, ax = plt.subplots(figsize=(2.8, 4))
fig.patch.set_facecolor('white')

rng = np.random.default_rng(42)
jit = (rng.random(len(z)) - 0.5) * 0.22
ax.scatter(1 + jit, z, s=28, color=col, alpha=0.5, zorder=3, linewidths=0)

q1, q2, q3 = np.percentile(z, [25, 50, 75])
iqr = q3 - q1
wlo = z[z >= q1 - 1.5*iqr].min()
whi = z[z <= q3 + 1.5*iqr].max()
bw  = 0.28
rect = plt.Rectangle((1-bw/2, q1), bw, q3-q1,
                      facecolor=col, alpha=0.35, edgecolor=col, lw=1.6, zorder=2)
ax.add_patch(rect)
ax.plot([1-bw/2, 1+bw/2], [q2, q2], color=col, lw=2.2, zorder=4)
ax.plot([1, 1], [wlo, q1], color=col, lw=1.1, zorder=2)
ax.plot([1, 1], [q3, whi], color=col, lw=1.1, zorder=2)

ax.axhline(0,    color='k',    ls='--', lw=0.8, alpha=0.5)
ax.axhline(1.96, color='gray', ls=':',  lw=0.8)
ax.text(1.22, 1.96, 'z = 1.96', va='center', fontsize=8, color='gray')

# Escala symlog: lineal cerca de 0, logarítmica para valores altos
ax.set_yscale('symlog', linthresh=2)
ax.set_ylim(z.min() - 1, z.max() * 1.18)
ax.set_yticks([-1, 0, 2, 5, 10, 20, 50, 100])
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.set_xlim(0.55, 1.55)
ax.set_xticks([1])
ax.set_xticklabels([r'$\theta$-PAC  (4–7 Hz)'], fontsize=11)
ax.set_ylabel('z-score  (MI vs. shuffled null)  [log scale]', fontsize=10)
ax.tick_params(direction='out', labelsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.tight_layout(pad=0.3)
fname = 'S4c_PAC_zMI_theta_vs_null.png'
fig.savefig(os.path.join(DIR_OUT, fname),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Guardada: {fname}')

# ═══════════════════════════════════════════════════════════════════════════
# FIG 2 — PAC × IES, 3 figuras individuales (sin texto)
# ═══════════════════════════════════════════════════════════════════════════
for b in range(3):
    x   = MI_late_ch[:, b] * 1e3
    y   = IES_b2
    col = COLS[b]

    rho, p = stats.spearmanr(x, y)
    res    = theilslopes(y, x)
    xl     = np.array([x.min(), x.max()])

    fig, ax = plt.subplots(figsize=(3.2, 3.2))
    fig.patch.set_facecolor('white')

    ax.scatter(x, y, s=32, color=col, alpha=0.55, linewidths=0, zorder=3)
    ax.plot(xl, res.slope * xl + res.intercept,
            color=col, lw=1.8, alpha=0.85, zorder=2)

    # rho y p dentro del plot
    p_str  = f'p = {p:.3f}' if p >= 0.001 else 'p < 0.001'
    star   = '  *' if p < 0.05 else ''
    ax.text(0.97, 0.05,
            f'ρ = {rho:+.3f}\n{p_str}{star}',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color=col if p < 0.05 else '#888888',
            fontweight='bold' if p < 0.05 else 'normal')

    ax.set_xlabel(BANDS_XLABEL[b], fontsize=10)
    ax.set_ylabel('IES  (chewing block)', fontsize=10)
    ax.tick_params(direction='out', labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout(pad=0.3)
    fname = f'S4c_PAC_behavior_{BANDS_LABEL[b]}.png'
    fig.savefig(os.path.join(DIR_OUT, fname),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Guardada: {fname}  |  rho={rho:+.3f}  p={p:.3f}')
