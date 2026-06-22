"""
figC_repro.py — Reproduce la estética figures_ok para FOOOF y TF, en CUADRADO.
FOOOF: PSD 4-condiciones (SEM + ajuste aperiódico punteado + bandas θ/α/β) · residual periódico · exponente.
TF: 3 mapas (casos / controles / interacción) en JET ±1.5, eje y lineal, contorno blanco del cluster,
líneas de referencia (t=0 dashed, ventanas/bandas dotted), colorbar. Todos 3.5×3.5 in PNG@300.
"""
import os, numpy as np, scipy.io as sio, h5py
import matplotlib.pyplot as plt
from _figstyle import COL, DIR_OUT, DIR_STATS, new_square, finish

TH, AL, BE = (4, 7), (8, 13), (13, 30)

# ════════════════════════════ FOOOF (4 condiciones) ════════════════════════
def _gr(h, grp, fld): return np.asarray(h["GR"][grp][fld])
def _ms(M):  # mean, sem sobre sujetos (M: subj x freq)
    return M.mean(0), M.std(0, ddof=1)/np.sqrt(M.shape[0])

def fooof_psd():
    with h5py.File(os.path.join(DIR_STATS, "FOOOF_Workspace_V1.mat"), "r") as h:
        f = _gr(h, "Cases", "f").ravel()
        fr = np.asarray(h["GR"]["FIT_RANGE"]).ravel() if "FIT_RANGE" in h["GR"] else np.array([3, 40])
        d = {}
        for grp in ("Cases", "Controls"):
            for c in ("Ch", "Nc"):
                d[(grp, c, "psd")] = 10*_gr(h, grp, f"PSD_log_{c}")
                d[(grp, c, "ap")]  = 10*_gr(h, grp, f"AP_{c}")
    m = (f >= fr[0]) & (f <= fr[1]); fv = f[m]
    fig, ax = new_square()
    ylo, yhi = None, None
    series = [("Controls", "Nc", COL["CTRL"], "-", "Controls – No-chew"),
              ("Controls", "Ch", COL["CTRL"], "--", "Controls – Chew"),
              ("Cases", "Nc", COL["NC"], "-", "Cases – No-chew"),
              ("Cases", "Ch", COL["CASE"], "--", "Cases – Chew")]
    means = [d[(g, c, "psd")][:, m].mean(0) for g, c, *_ in series]
    ylo, yhi = np.floor(min(x.min() for x in means)-1), np.ceil(max(x.max() for x in means)+1)
    for (lo, hi), lbl in [(TH, "θ"), (AL, "α"), (BE, "β")]:
        ax.axvspan(lo, hi, color="0.92", alpha=0.7, zorder=0)
        ax.text((lo+hi)/2, yhi-0.6, lbl, ha="center", fontsize=9, color="0.5")
    for g, c, col, ls, lbl in series:
        mu, se = _ms(d[(g, c, "psd")][:, m])
        ax.fill_between(fv, mu-se, mu+se, color=col, alpha=0.15, lw=0, zorder=1)
        ax.plot(fv, mu, ls=ls, color=col, lw=1.8, label=lbl, zorder=3)
        apmu = d[(g, c, "ap")][:, m].mean(0)
        ax.plot(fv, apmu, ls=":", color=col, lw=1.1, alpha=0.8, zorder=2)
    ax.plot([], [], "k:", lw=1.1, label="aperiodic fit")
    ax.set_xlim(fr[0], fr[1]); ax.set_ylim(ylo, yhi)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Power (dB)")
    ax.set_title("Power spectral density", pad=8)
    ax.legend(loc="lower left", fontsize=7, frameon=False, handlelength=1.6)
    finish(fig, ax, "fooof_PSD_aperiodic_fit_groups")

def fooof_residual():
    with h5py.File(os.path.join(DIR_STATS, "FOOOF_Workspace_V1.mat"), "r") as h:
        f = _gr(h, "Cases", "f").ravel()
        d = {(g, c): _gr(h, g, f"Res_{c}") for g in ("Cases", "Controls") for c in ("Ch", "Nc")}
    m = (f >= 2) & (f <= 30); fv = f[m]
    fig, ax = new_square()
    ax.axvspan(*TH, color="0.92", alpha=0.7, zorder=0)  # gray, consistent with fooof_psd
    for (g, c), col, ls, lbl in [(("Controls", "Nc"), COL["CTRL"], "-", "Controls – No-chew"),
                                 (("Controls", "Ch"), COL["CTRL"], "--", "Controls – Chew"),
                                 (("Cases", "Nc"), COL["NC"], "-", "Cases – No-chew"),
                                 (("Cases", "Ch"), COL["CASE"], "--", "Cases – Chew")]:
        mu, se = _ms(d[(g, c)][:, m])
        ax.fill_between(fv, mu-se, mu+se, color=col, alpha=0.15, lw=0)
        ax.plot(fv, mu, ls=ls, color=col, lw=1.8, label=lbl)
    ax.axhline(0, color="k", lw=0.7, ls=":")
    ax.set_xlim(2, 30); ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Periodic power (a.u.)")
    ax.set_title("Periodic component (FOOOF residual)", pad=8)
    ax.legend(loc="upper right", fontsize=7, frameon=False, handlelength=1.6)
    finish(fig, ax, "fooof_periodic_residual")

# ════════════════════════════ TF (jet, cuadrado) ═══════════════════════════
def tf_panels():
    TF = sio.loadmat(os.path.join(DIR_STATS, "S2c_TF_GroupFigure.mat"), squeeze_me=True)
    frex = np.asarray(TF["FREX_TF"], float); t = np.asarray(TF["times_anal"], float)
    panels = [("map_cas", "mask_cas", "tf_cases_chew_minus_nochew"),
              ("map_ctr", "mask_ctr", "tf_controls_chew_minus_nochew"),
              ("map_int", "mask_int", "tf_interaction")]
    for mapk, maskk, name in panels:
        M = np.asarray(TF[mapk], float)
        fig, ax = new_square()
        levels = np.linspace(-1.5, 1.5, 61)
        # sin extend -> colorbar con extremos PLANOS (como MATLAB), datos recortados a ±1.5
        cf = ax.contourf(t, frex, np.clip(M, -1.4999, 1.4999), levels=levels, cmap="jet", extend="neither")
        if maskk in TF and np.any(TF[maskk]):
            ax.contour(t, frex, np.asarray(TF[maskk], float), levels=[0.5], colors="white", linewidths=2.0)
        ax.axvline(0, color="k", lw=1.4, ls="--")
        for xv in (100, 900, 1300):
            ax.axvline(xv, color="white", lw=0.8, ls=":", alpha=0.8)
        for yv in (7, 13):
            ax.axhline(yv, color="white", lw=0.8, ls=":", alpha=0.8)
        ax.set_ylim(2, 40); ax.set_yticks([5, 10, 15, 20, 25, 30, 35, 40])
        ax.set_xlim(t.min(), t.max())
        ax.set_xlabel("Time (ms)"); ax.set_ylabel("Frequency (Hz)")
        cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, ticks=[-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
        cb.set_label("Chew − No-chew (dB)", fontsize=9); cb.ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(True)   # TF lleva marco completo (estilo figures_ok)
        fig.tight_layout(pad=0.6)
        out = os.path.join(DIR_OUT, name+".png")
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
        print(f"  saved {name}.png")

if __name__ == "__main__":
    print("== figC_repro: FOOOF + TF (estética figures_ok, cuadrado) ==")
    fooof_psd(); fooof_residual(); tf_panels()
    print("Listo.")
