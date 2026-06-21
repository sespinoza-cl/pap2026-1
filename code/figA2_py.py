"""
figA2_py.py — Paneles principales (Python, datos de grupo precomputados) Paper2 V1.
A4-A6 mapas TF (S2c), A8-A10 FOOOF (FOOOF_Workspace_V1.mat GR), A16 comodulograma (S4g).
Todo cuadrado, Arial homogénea, PNG @300 dpi, etiquetas SIN TMD.
"""
import os, numpy as np, scipy.io as sio, h5py
import matplotlib.pyplot as plt
from _figstyle import COL, DIR_OUT, DIR_STATS, new_square, finish

# ── TF maps ──────────────────────────────────────────────────────────────────
TF = sio.loadmat(os.path.join(DIR_STATS, "S2c_TF_GroupFigure.mat"), squeeze_me=True)
frex = np.asarray(TF["FREX_TF"], float); t = np.asarray(TF["times_anal"], float)

def _tf_panel(mapk, maskk, title, name):
    M = np.asarray(TF[mapk], float)
    vmax = np.nanpercentile(np.abs(M), 98)
    fig, ax = new_square()
    pc = ax.pcolormesh(t, frex, M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    if maskk in TF and np.any(TF[maskk]):
        ax.contour(t, frex, np.asarray(TF[maskk], float), levels=[0.5],
                   colors="k", linewidths=1.0)
    ax.set_yscale("log"); ax.set_ylim(frex.min(), 30)
    ax.set_yticks([4, 7, 13, 30]); ax.set_yticklabels([4, 7, 13, 30])
    ax.set_xlabel("Time (ms)"); ax.set_ylabel("Frequency (Hz)")
    ax.set_title(title, pad=8)
    cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Chew − No-chew (dB)", fontsize=9); cb.ax.tick_params(labelsize=8)
    ax.axvline(0, color="k", lw=0.7, ls=":")
    finish(fig, ax, name)

def A4(): _tf_panel("map_cas", "mask_cas", "Cases: Chew − No-chew", "tf_cases_chew_minus_nochew")
def A5(): _tf_panel("map_ctr", "mask_ctr", "Controls: Chew − No-chew", "tf_controls_chew_minus_nochew")
def A6(): _tf_panel("map_int", "mask_int", "Interaction (Group × Cond)", "tf_interaction")

# ── FOOOF (GR struct, v7.3) ──────────────────────────────────────────────────
def _gr(h, grp, fld):
    return np.asarray(h["GR"][grp][fld])

def _load_fooof():
    h = h5py.File(os.path.join(DIR_STATS, "FOOOF_Workspace_V1.mat"), "r")
    return h

def A8():
    h = _load_fooof()
    f = _gr(h, "Cases", "f").ravel()
    fig, ax = new_square()
    for grp, col, lbl in [("Cases", COL["CASE"], "Cases"), ("Controls", COL["CTRL"], "Controls")]:
        psd = _gr(h, grp, "PSD_log_Ch").mean(0)
        ap  = _gr(h, grp, "AP_Ch").mean(0)
        ax.plot(f, psd, color=col, lw=2.0, label=f"{lbl} (chew)")
        ax.plot(f, ap, color=col, lw=1.2, ls="--", alpha=0.8)
    ax.set_xscale("log"); ax.set_xlim(2, 40)
    ax.set_xticks([2, 4, 7, 13, 30]); ax.set_xticklabels([2, 4, 7, 13, 30])
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("log power (a.u.)")
    ax.set_title("PSD + aperiodic fit", pad=8)
    ax.legend(loc="upper right", frameon=False, handlelength=1.4)
    finish(fig, ax, "fooof_PSD_aperiodic_fit_groups")
    h.close()

def _paired(ax, nc, ch, col, label):
    nc, ch = nc.ravel(), ch.ravel()
    for i in range(len(nc)):
        ax.plot([0, 1], [nc[i], ch[i]], color=col, alpha=0.13, lw=0.7)
    m = [nc.mean(), ch.mean()]
    se = [nc.std(ddof=1)/np.sqrt(len(nc)), ch.std(ddof=1)/np.sqrt(len(ch))]
    ax.errorbar([0, 1], m, yerr=se, color=col, lw=2.2, marker="o", ms=6, capsize=3, label=label)

def A9():
    h = _load_fooof()
    fig, ax = new_square()
    _paired(ax, _gr(h, "Controls", "exp_Nc"), _gr(h, "Controls", "exp_Ch"), COL["CTRL"], "Controls")
    _paired(ax, _gr(h, "Cases", "exp_Nc"), _gr(h, "Cases", "exp_Ch"), COL["CASE"], "Cases")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["No-chew", "Chew"]); ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("Aperiodic exponent χ")
    ax.set_title("Aperiodic flattening", pad=8)
    ax.text(0.5, 0.04, "Cases p<0.0001 · Controls n.s.", transform=ax.transAxes,
            ha="center", fontsize=8.5, color="0.35")
    ax.legend(loc="upper right", frameon=False, handlelength=1.2)
    finish(fig, ax, "fooof_exponent_cases_vs_controls")
    h.close()

def A10():
    h = _load_fooof()
    f = _gr(h, "Cases", "f").ravel()
    res = _gr(h, "Cases", "Res_Ch")
    m = res.mean(0); se = res.std(0, ddof=1)/np.sqrt(res.shape[0])
    fig, ax = new_square()
    ax.axvspan(4, 7, color=COL["THETA"], alpha=0.12)
    ax.plot(f, m, color=COL["THETA"], lw=2.0)
    ax.fill_between(f, m-se, m+se, color=COL["THETA"], alpha=0.22, lw=0)
    ax.set_xlim(2, 30); ax.set_xticks([2, 4, 7, 13, 30]); ax.set_xticklabels([2, 4, 7, 13, 30])
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Periodic power (a.u.)")
    ax.set_title("Genuine θ peak (Cases, chew)", pad=8)
    ax.text(0.97, 0.93, "peak 97%\nCF≈6.6 Hz", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color=COL["THETA"])
    ax.axhline(0, color="k", lw=0.7, ls=":")
    finish(fig, ax, "fooof_theta_peak_overlay")
    h.close()

# ── Comodulograma descriptivo (S4g) ──────────────────────────────────────────
def A16():
    g = sio.loadmat(os.path.join(DIR_STATS, "S4g_comodulogram.mat"), squeeze_me=True)
    MI = np.asarray(g["MI_mean"], float)        # (amp 27, phase 11) ? -> shape (11,27) per dump
    fp = np.asarray(g["F_PHASE"], float); fa = np.asarray(g["F_AMP"], float)
    # MI_mean dump was (11,27)=(phase,amp); transpose so amp on y
    if MI.shape == (len(fp), len(fa)):
        MI = MI.T
    fig, ax = new_square()
    pc = ax.pcolormesh(fp, fa, MI, cmap="magma", shading="auto")
    ax.set_xlabel("Phase frequency (Hz)"); ax.set_ylabel("Amplitude frequency (Hz)")
    ax.set_title("Comodulogram (descriptive)", pad=8)
    cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("MI", fontsize=9); cb.ax.tick_params(labelsize=8)
    ax.text(0.5, 1.0, "broadband fa 4–30 Hz", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8, color="0.4")
    finish(fig, ax, "comodulogram_descriptive_MI")

if __name__ == "__main__":
    print("== figA2_py: TF + FOOOF + comodulograma ==")
    for fn in (A4, A5, A6, A8, A9, A10, A16):
        try:
            fn();
        except Exception as e:
            print(f"  [ERR] {fn.__name__}: {e}")
    print("Listo.")
