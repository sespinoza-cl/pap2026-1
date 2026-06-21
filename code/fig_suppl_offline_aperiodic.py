"""
fig_suppl_offline_aperiodic.py — Figura suplementaria de convergencia cross-study:
exponente aperiódico que se aplana con masticación OFFLINE (Paper 1, re-análisis) y ONLINE (Paper 2).
Dos paneles cuadrados individuales (estándar _figstyle), PNG @300 dpi.
"""
import os, numpy as np, h5py
import matplotlib.pyplot as plt
from _figstyle import COL, DIR_OUT, DIR_STATS, new_square, finish, sig_bracket, pstars

# --- OFFLINE (Paper 1): npz exp_chew / exp_nochew (n=30) ---
z = np.load(os.path.join(DIR_STATS, "rev_paper1_aperiodic.npz"))
off_ch, off_nc = np.asarray(z["exp_chew"], float), np.asarray(z["exp_nochew"], float)

# --- ONLINE (Paper 2): FOOOF GR Cases exp_Ch / exp_Nc (n=31) ---
with h5py.File(os.path.join(DIR_STATS, "FOOOF_Workspace_V1.mat"), "r") as h:
    on_ch = np.asarray(h["GR"]["Cases"]["exp_Ch"]).ravel()
    on_nc = np.asarray(h["GR"]["Cases"]["exp_Nc"]).ravel()

def paired_panel(nc, ch, title, name, pstr):
    fig, ax = new_square()
    for i in range(len(nc)):
        ax.plot([0, 1], [nc[i], ch[i]], color=COL["CASE"], alpha=0.13, lw=0.7)
    m = [nc.mean(), ch.mean()]
    se = [nc.std(ddof=1)/np.sqrt(len(nc)), ch.std(ddof=1)/np.sqrt(len(ch))]
    ax.errorbar([0, 1], m, yerr=se, color=COL["CASE"], lw=2.4, marker="o", ms=7, capsize=3)
    d = ch.mean() - nc.mean()
    top = max(np.max(nc), np.max(ch))
    sig_bracket(ax, 0, 1, top*1.02, f"{pstr}  (Δχ={d:+.2f})")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["No-chew", "Chew"]); ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("Aperiodic exponent χ"); ax.set_title(title, pad=8)
    finish(fig, ax, name)

# OFFLINE: Wilcoxon p=0.005 (21/30 flatten)
paired_panel(off_nc, off_ch,
             "Offline (pre-task chewing)\nEspinoza 2025 re-analysis (n=30)",
             "fooofS_aperiodic_offline_paper1", "**")
# ONLINE: p<0.0001
paired_panel(on_nc, on_ch,
             "Online (chewing during task)\nPresent study, cases (n=31)",
             "fooofS_aperiodic_online_paper2", "***")
print("OFFLINE Δχ=%.3f (21/30, p=0.005) · ONLINE Δχ=%.3f (p<0.0001)"
      % (off_ch.mean()-off_nc.mean(), on_ch.mean()-on_nc.mean()))
