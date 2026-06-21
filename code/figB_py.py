"""
figB_py.py — Paneles SUPLEMENTARIOS (Python) Paper2 V1, estándar cuadrado/PNG/Arial.
B2 spatial frontal vs temporal (S6), B3 masseterRMS×Δθ (S6b), B6 dose-response (S7),
B7 MI por ventana (v1_S4b), B11/B12 comodulograma 2D cluster/zmap (S4h).
"""
import os, numpy as np, scipy.io as sio, h5py
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
from _figstyle import (COL, DIR_OUT, DIR_STATS, new_square, finish,
                       box_scatter, refline, robust_scatter)

def h5get(p, *keys):
    out = {}
    with h5py.File(os.path.join(DIR_STATS, p), "r") as h:
        for k in keys:
            out[k] = np.asarray(h[k]).ravel() if np.asarray(h[k]).ndim <= 2 else np.asarray(h[k])
    return out

def scatter_corr(ax, x, y, color, xlabel, ylabel, title, name):
    x, y = np.ravel(x).astype(float), np.ravel(y).astype(float)
    ax.scatter(x, y, s=22, color=color, alpha=0.6, edgecolor="none")
    lr = linregress(x, y)
    rho, pval = spearmanr(x, y)              # canónico: Spearman
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, lr.intercept + lr.slope*xs, color="k", lw=1.2, ls="--", alpha=0.7)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title, pad=8)
    tag = "n.s." if pval >= 0.05 else f"p={pval:.3f}"
    ax.text(0.04, 0.95, f"ρ={rho:+.2f}  {tag}  (Spearman)", transform=ax.transAxes,
            va="top", fontsize=9, color="0.25")
    return rho, pval

# ── B2 — disociación espacial frontal-medial vs temporal ─────────────────────
def B2():
    with h5py.File(os.path.join(DIR_STATS, "S6_artifact_controls.mat"), "r") as h:
        t_el = np.asarray(h["t_el"]).ravel()
        p_fdr = np.asarray(h["p_fdr"]).ravel()
        fm = np.asarray(h["fm_idx"]).ravel().astype(int) - 1
        tp = np.asarray(h["tp_idx"]).ravel().astype(int) - 1
    fig, ax = new_square()
    groups = [("Frontal-medial", fm, COL["THETA"]), ("Temporal", tp, COL["NC"])]
    for j, (lbl, idx, col) in enumerate(groups):
        tv = t_el[idx]; pv = p_fdr[idx]
        xj = np.random.default_rng(3+j).normal(j, 0.06, len(tv))
        sigm = pv < 0.05
        ax.scatter(xj[~sigm], tv[~sigm], s=34, facecolor="white", edgecolor=col, lw=1.3, zorder=2)
        ax.scatter(xj[sigm], tv[sigm], s=40, color=col, edgecolor="k", lw=0.6, zorder=3)
        ax.text(j, t_el[np.r_[fm, tp]].max()*1.04, f"{int(sigm.sum())}/{len(tv)} sig",
                ha="center", fontsize=8.5, color=col)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Frontal-\nmedial", "Temporal"])
    ax.set_xlim(-0.5, 1.5); ax.set_ylabel("t (Chew−No-chew, θ)")
    ax.set_title("Spatial dissociation (A1)", pad=8)
    finish(fig, ax, "artifactS_spatial_frontal_vs_temporal")

# ── B3 — EMG masetero real × Δθ (A3, n.s.) ───────────────────────────────────
def B3():
    with h5py.File(os.path.join(DIR_STATS, "S6b_A3_muscle.mat"), "r") as h:
        rms = np.asarray(h["emg_rms"]).ravel(); dth = np.asarray(h["dtheta"]).ravel()
    fig, ax = new_square()
    robust_scatter(ax, rms, dth, COL["BETA"], "Masseter EMG RMS", "Δθ power (Chew−Nc, dB)",
                   "Muscle artifact control (A3)")
    finish(fig, ax, "artifactS_masseterRMS_x_dtheta")

# ── B6 — dosis-respuesta: engagement × Δθ (n.s.) ─────────────────────────────
def B6():
    d = sio.loadmat(os.path.join(DIR_STATS, "S7_dose_response.mat"), squeeze_me=True)
    fig, ax = new_square()
    robust_scatter(ax, d["eng"], d["dth"], COL["CASE"], "Chewing engagement",
                   "Δθ power (Chew−Nc, dB)", "Dose-response (S7)")
    finish(fig, ax, "doseS_engagement_x_dtheta")

# ── B7 — MI θ por ventana (Base/Early/Late) ──────────────────────────────────
def B7():
    d = sio.loadmat(os.path.join(DIR_STATS, "v1_S4b_PAC_ROI.mat"), squeeze_me=True)
    MIw = np.asarray(d["MIw_ch"], float)        # (band, window, subj) esperado
    # localizar dim banda(=3) y window(=3): subj=31 es la última
    th = MIw[0]                                  # theta -> (window=3, subj=31)
    labels = ["Base", "Early", "Late"]
    med = np.median(th, 1)                        # mediana (canónico; robusta a outliers)
    q1 = np.percentile(th, 25, 1); q3 = np.percentile(th, 75, 1)
    fig, ax = new_square()
    ax.fill_between([0, 1, 2], q1, q3, color=COL["THETA"], alpha=0.18, lw=0)
    ax.plot([0, 1, 2], med, color=COL["THETA"], lw=2.3, marker="o", ms=7, zorder=3)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(labels); ax.set_xlim(-0.3, 2.3)
    ax.set_ylim(0, max(q3)*1.15)
    ax.set_ylabel("MI θ (Tort, median ± IQR)"); ax.set_title("θ-PAC by time window", pad=8)
    ax.text(0.5, 0.97, "Friedman p=0.0007 · Late>Early p=0.001", transform=ax.transAxes,
            ha="center", va="top", fontsize=8.5, color="0.35")
    finish(fig, ax, "pacS_MI_by_window")

# ── B11/B12 — comodulograma 2D estadístico ───────────────────────────────────
def _como2d(field, title, name, cmap, mask=True):
    with h5py.File(os.path.join(DIR_STATS, "S4h_comodulogram_stats.mat"), "r") as h:
        M = np.asarray(h[field]); fa = np.asarray(h["F_AMP"]).ravel()
        fp = np.asarray(h["F_PHASE"]).ravel(); sig = np.asarray(h["sig_mask"])
    # h5py lee (11,27) para (27,11) MATLAB -> orientar a (amp,phase)=(27,11)
    if M.shape == (len(fp), len(fa)): M = M.T
    if sig.shape == (len(fp), len(fa)): sig = sig.T
    fig, ax = new_square()
    pc = ax.pcolormesh(fp, fa, M, cmap=cmap, shading="auto")
    if mask and np.any(sig):
        ax.contour(fp, fa, sig.astype(float), levels=[0.5], colors="k", linewidths=1.2)
    ax.set_xlabel("Phase frequency (Hz)"); ax.set_ylabel("Amplitude frequency (Hz)")
    ax.set_title(title, pad=8)
    cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04); cb.ax.tick_params(labelsize=8)
    finish(fig, ax, name)

def B11(): _como2d("z_mean", "2D comodulogram — cluster", "comodulogramS_2D_cluster", "magma", True)
def B12(): _como2d("z_mean", "2D comodulogram — z map", "comodulogramS_2D_zmap", "viridis", False)

if __name__ == "__main__":
    print("== figB_py: suplementarios ==")
    for fn in (B2, B3, B6, B7, B11, B12):
        try: fn()
        except Exception as e: print(f"  [ERR] {fn.__name__}: {e}")
    print("Listo.")
