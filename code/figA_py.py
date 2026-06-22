"""
figA_py.py — Paneles principales (Python) Paper2 V1, marco DUAL θ+β-CMC.
Genera PNG cuadrados individuales en outputs/figures/ desde .mat CANÓNICOS.
A1 RT, A2 IES, A3 masseter manip-check, A11 zMI θ doble null, A12 2x2 mecanismo,
A13 rose θ, A14 band comparison (θ/α/β), A15 rose β.
"""
import os, numpy as np, scipy.io as sio
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
from _figstyle import (COL, LBL_CASE, LBL_CTRL, DIR_OUT, DIR_STATS, DIR_DATA,
                       new_square, finish, pstars, sig_bracket,
                       box_scatter, refline, symlog_zmi, robust_scatter, paired_box, rose_v2)

BEH = sio.loadmat(os.path.join(DIR_DATA, "v1_S1_behavior_stats.mat"), squeeze_me=True)
CHW = sio.loadmat(os.path.join(DIR_STATS, "chew_engagement_check.mat"), squeeze_me=True)
PAC = sio.loadmat(os.path.join(DIR_STATS, "v1_S4b_PAC_ROI.mat"), squeeze_me=True)
S8  = sio.loadmat(os.path.join(DIR_STATS, "S8_controls_PAC.mat"), squeeze_me=True)
TH, AL, BE = 0, 1, 2  # columnas banda

def _paired_group(ax, b1, b2, color, label):
    x = np.array([0, 1.0])
    for i in range(len(b1)):
        ax.plot(x, [b1[i], b2[i]], color=color, alpha=0.13, lw=0.7, zorder=1)
    m = [np.mean(b1), np.mean(b2)]
    se = [np.std(b1, ddof=1)/np.sqrt(len(b1)), np.std(b2, ddof=1)/np.sqrt(len(b2))]
    ax.errorbar(x, m, yerr=se, color=color, lw=2.2, marker="o", ms=6,
                capsize=3, zorder=3, label=label)
    return max(np.max(b1), np.max(b2))

# ── A1 — RT por bloque × grupo ───────────────────────────────────────────────
def A1():
    fig, ax = new_square()
    paired_box(ax, [
        {"label": "Controls", "color": COL["CTRL"], "nc": BEH["rt_ctr_b1"], "ch": BEH["rt_ctr_b2"], "stars": None},
        {"label": "Cases", "color": COL["CASE"], "nc": BEH["rt_cas_b1"], "ch": BEH["rt_cas_b2"],
         "stars": pstars(float(BEH["p_rt_sign"]))},
    ], ylabel="Median RT (ms)", title="Reaction time")
    # Reviewer-requested: n.s. bar between NoChew baselines (Controls x=0, Cases x=3)
    _, p_bl = mannwhitneyu(BEH["rt_ctr_b1"], BEH["rt_cas_b1"], alternative="two-sided")
    all_rt = np.concatenate([BEH["rt_ctr_b1"], BEH["rt_cas_b1"],
                              BEH["rt_ctr_b2"], BEH["rt_cas_b2"]])
    ytop = float(np.max(all_rt))
    sig_bracket(ax, 0, 3, ytop * 1.05, pstars(p_bl))
    ax.set_ylim(top=ytop * 1.20)
    finish(fig, ax, "behavior_RT_by_block_group")

# ── A2 — IES por bloque × grupo ──────────────────────────────────────────────
def A2():
    fig, ax = new_square()
    paired_box(ax, [
        {"label": "Controls", "color": COL["CTRL"], "nc": BEH["ies_ctr_b1"], "ch": BEH["ies_ctr_b2"], "stars": None},
        {"label": "Cases", "color": COL["CASE"], "nc": BEH["ies_cas_b1"], "ch": BEH["ies_cas_b2"],
         "stars": pstars(float(BEH["p_ies_sign"]))},
    ], ylabel="IES (ms)", title="Inverse efficiency score")
    # Reviewer-requested: n.s. bar between NoChew baselines (Controls x=0, Cases x=3)
    _, p_bl = mannwhitneyu(BEH["ies_ctr_b1"], BEH["ies_cas_b1"], alternative="two-sided")
    all_ies = np.concatenate([BEH["ies_ctr_b1"], BEH["ies_cas_b1"],
                               BEH["ies_ctr_b2"], BEH["ies_cas_b2"]])
    ytop = float(np.max(all_ies))
    sig_bracket(ax, 0, 3, ytop * 1.05, pstars(p_bl))
    ax.set_ylim(top=ytop * 1.20)
    finish(fig, ax, "behavior_IES_by_block_group")

# ── A3 — Manipulation check: masetero Ch−Nc por grupo ────────────────────────
def A3():
    cas, ctr = np.asarray(CHW["casd"], float), np.asarray(CHW["ctrd"], float)
    fig, ax = new_square()
    box_scatter(ax, [ctr, cas], [COL["CTRL"], COL["CASE"]],
                [f"Control\n(n={len(ctr)})", f"Cases\n(n={len(cas)})"], s=26)
    refline(ax, 0, style="--", color="k", lw=0.8, alpha=0.5)
    ax.set_ylabel("Masseter Ch−Nc (dB)")
    ax.set_title("Jaw-muscle chewing rhythm (masseter SNR)", pad=8)
    ax.text(1, np.max(cas)*1.04, "***", ha="center", va="bottom", fontsize=13, color=COL["CASE"])
    ax.text(0, np.max(cas)*0.55, "n.s.", ha="center", color=COL["CTRL"], fontsize=9)
    finish(fig, ax, "manipcheck_masseter_SNR_ChNc_group")

def _zbox(ax, series, colors, labels, thr=1.96, symlog=True):
    """Box+scatter estilo figures_ok, con symlog para zMI de rango amplio.
    Significancia visual = puntos sobre la línea z=1.96 (conteos en el caption, no en el plot)."""
    box_scatter(ax, series, colors, labels, s=26)
    allv = np.concatenate([np.asarray(v, float) for v in series])
    if symlog:
        symlog_zmi(ax, allv, thr=2)
    refline(ax, thr, label=f"z = {thr}", style="--", color=COL["SIG"], alpha=0.8)

# ── A11 — zMI θ vs doble null (circ + AAFT) ──────────────────────────────────
def A11():
    fig, ax = new_square()
    s = [PAC["zC_ch"][:, TH], PAC["zA_ch"][:, TH]]
    _zbox(ax, s, [COL["THETA"], COL["THETA"]],
          ["circular-shift\nnull", "AAFT\nnull"])
    ax.set_ylabel("zMI θ (z-score, symlog)")
    ax.set_title("θ-PAC over double null", pad=8)
    finish(fig, ax, "pac_zMI_theta_doublenull")

# ── A12 — 2x2 mecanismo: Ch vs Nc + control negativo ─────────────────────────
def A12():
    fig, ax = new_square()
    s = [PAC["zC_ch"][:, TH], PAC["zC_nc"][:, TH], S8["zC"][:, TH]]
    _zbox(ax, s, [COL["CASE"], COL["NC"], COL["CTRL"]],
          ["Cases\nChew", "Cases\nNo-chew", "Controls\nNo-chew"])
    ax.set_ylabel("zMI θ (z-score, symlog)")
    ax.set_title("θ-PAC by condition", pad=8)
    ytop = max(np.asarray(s[0]).max(), np.asarray(s[1]).max())   # bracket por encima de TODO
    sig_bracket(ax, 0, 1, ytop*1.10, "***")
    finish(fig, ax, "pac_zMI_theta_chew_vs_nochew_negctrl")

# ── A14 — comparación de bandas (θ/α/β) doble null ───────────────────────────
def A14():
    fig, ax = new_square()
    s = [PAC["zC_ch"][:, TH], PAC["zC_ch"][:, AL], PAC["zC_ch"][:, BE]]
    _zbox(ax, s, [COL["THETA"], COL["ALPHA"], COL["BETA"]],
          ["θ\n4–7", "α\n8–13", "β\n13–30"])
    ax.set_ylabel("zMI (z-score, symlog)")
    ax.set_title("Broadband PAC", pad=8)
    finish(fig, ax, "pac_zMI_band_comparison_doublenull")

# ── A13 / A15 — Rayleigh rose ────────────────────────────────────────────────
def _rose(band_idx, color, band_lbl, name):
    ang = np.asarray(PAC["pref_ch"][:, band_idx], float)
    R = float(np.asarray(PAC["rayl_R"])[band_idx]); p = float(np.asarray(PAC["rayl_p"])[band_idx])
    rose_v2(ang, R, color, band_lbl, pstars(p), name)

def A13():    _rose(TH, COL["THETA"], r"$\theta$", "pac_rayleigh_theta_rose")
def Aalpha(): _rose(AL, COL["ALPHA"], r"$\alpha$", "pac_rayleigh_alpha_rose")
def A15():    _rose(BE, COL["BETA"],  r"$\beta$",  "pac_rayleigh_beta_rose")

if __name__ == "__main__":
    print("== figA_py: paneles principales ==")
    for fn in (A1, A2, A3, A11, A12, A14, A13, Aalpha, A15):
        fn()
    print("Listo. PNGs en outputs/figures/")
