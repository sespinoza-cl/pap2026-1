"""
S4d_Rayleigh_rose.py
====================
Rose histogram + vector resultante (negro), 3 bandas.
Figuras individuales sin texto para montaje en Inkscape.

Salida (Analysis_V1_Final/outputs/figures/):
  S4d_Rayleigh_rose_theta_4to7Hz.png
  S4d_Rayleigh_rose_alpha_8to12Hz.png
  S4d_Rayleigh_rose_beta_13to20Hz.png
"""

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

ROOT     = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
PAC_FILE = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "stats",
                        "v1_S4b_PAC_ROI.mat")
DIR_OUT  = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "figures")
os.makedirs(DIR_OUT, exist_ok=True)

with h5py.File(PAC_FILE, 'r') as f:
    pref_ch = np.array(f['pref_ch'])   # (3, 31)

BANDS_LABEL = ['theta_4to7Hz', 'alpha_8to12Hz', 'beta_13to20Hz']
COLS        = ['#D55E00', '#56B4E9', '#CC79A7']  # Okabe-Ito: vermillion/skyblue/purple
N_BINS      = 12

def rayleigh(pv):
    C  = np.mean(np.exp(1j * pv))
    R  = np.abs(C)
    mu = np.angle(C)
    return R, mu

bin_edges   = np.linspace(-np.pi, np.pi, N_BINS + 1)
bin_width   = 2 * np.pi / N_BINS
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

for b in range(3):
    pv  = pref_ch[b, :]
    pv  = pv[~np.isnan(pv)]
    col = COLS[b]

    counts, _ = np.histogram(pv, bins=bin_edges)
    R, mu     = rayleigh(pv)
    r_vec     = R * counts.max()

    fig, ax = plt.subplots(figsize=(4, 4),
                           subplot_kw=dict(projection='polar'))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Rose histogram
    bars = ax.bar(bin_centers, counts,
                  width=bin_width * 0.88,
                  color=col, alpha=0.60,
                  edgecolor='white', linewidth=0.7,
                  zorder=2)
    bars[np.argmax(counts)].set_alpha(0.85)

    # Círculo de referencia (radio = max_count)
    theta_ref = np.linspace(0, 2 * np.pi, 300)
    ax.plot(theta_ref, np.full(300, counts.max()),
            color='lightgray', lw=0.6, ls='--', zorder=1)

    # Vector resultante — NEGRO
    ax.annotate('',
                xy=(mu, r_vec),
                xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='black',
                                lw=2.2, mutation_scale=16))

    # Etiquetas de ángulo en notación π
    ax.set_xticks(np.radians([0, 90, 180, 270]))
    ax.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$'], fontsize=10)

    # Etiquetas radiales (conteo de sujetos)
    ymax = counts.max() + 1.5
    ax.set_ylim(0, ymax)
    yticks = np.arange(2, counts.max() + 1, 2)
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(y) for y in yticks], fontsize=7.5, color='gray')

    ax.spines['polar'].set_color('lightgray')
    ax.tick_params(pad=3)

    plt.tight_layout(pad=0.3)
    fname = f'S4d_Rayleigh_rose_{BANDS_LABEL[b]}.png'
    fig.savefig(os.path.join(DIR_OUT, fname),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Guardada: {fname}  |  R={R:.3f}  angle={np.degrees(mu)%360:.1f}°')
