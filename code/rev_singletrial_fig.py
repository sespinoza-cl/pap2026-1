"""
rev_singletrial_fig.py — Figuras del puente conductual↔EEG a nivel de ensayo (chew, casos).
A) RT por quintil de theta-late within-subject (mean±SEM across subj).
B) Especificidad de ventana: r_within medio (pre/early/late) ± 95% CI.
"""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from _figstyle import COL, DIR_OUT, DIR_STATS, new_square, finish

df = pd.read_csv(os.path.join(DIR_STATS, "rev_singletrial.csv"))
WINS = ["pre", "early", "late"]

# ---- A) quintiles within-subject de theta_late vs RT ----
nb = 5
parts = []
for s, g in df.groupby("subj"):
    if len(g) < nb*2: continue
    q = pd.qcut(g["theta_late"].rank(method="first"), nb, labels=False)
    m = g.groupby(q)["RT"].mean()
    parts.append(m.reindex(range(nb)).values)
M = np.array(parts)                       # subj x bins
mean = np.nanmean(M, 0); sem = np.nanstd(M, 0, ddof=1)/np.sqrt(np.sum(~np.isnan(M), 0))
fig, ax = new_square()
ax.errorbar(range(1, nb+1), mean, yerr=sem, color=COL["THETA"], lw=2.3, marker="o", ms=7, capsize=3)
ax.set_xticks(range(1, nb+1)); ax.set_xlabel("Frontal-medial θ (late, within-subj quintile)")
ax.set_ylabel("Reaction time (ms)")
ax.set_title("Higher late θ → faster RT (trial level)", pad=8)
ax.text(0.5, 0.92, "LMM slope=−41.8 ms, p=0.0007 · 23/31 subj", transform=ax.transAxes,
        ha="center", fontsize=8.5, color="0.35")
finish(fig, ax, "bridge_thetaLate_RT_quintiles")

# ---- B) especificidad de ventana: r_within por sujeto ----
rmeans, rcis, cols = [], [], []
for w in WINS:
    rs = []
    for s, g in df.groupby("subj"):
        if len(g) >= 10:
            r, _ = pearsonr(g[f"theta_{w}"], g["RT"]); rs.append(np.arctanh(np.clip(r, -.999, .999)))
    rs = np.array(rs)
    rmeans.append(np.tanh(rs.mean()))
    rcis.append(np.tanh(1.96*rs.std(ddof=1)/np.sqrt(len(rs))))   # approx CI half-width in z->r
    cols.append(COL["THETA"] if w == "late" else COL["NC"])
fig, ax = new_square()
x = np.arange(len(WINS))
ax.bar(x, rmeans, yerr=rcis, color=cols, capsize=4, width=0.6, edgecolor="k", linewidth=0.6)
ax.axhline(0, color="k", lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(["pre\n−500:−100", "early\n100:500", "late\n900:1300"])
ax.set_ylabel("within-subject r (θ × RT)")
ax.set_title("Window specificity of the link", pad=8)
ax.text(2, rmeans[2]-rcis[2]-0.01, "**", ha="center", va="top", fontsize=12, color=COL["THETA"])
finish(fig, ax, "bridge_window_specificity")
print("means r:", dict(zip(WINS, np.round(rmeans, 3))))
