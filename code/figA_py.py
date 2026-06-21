"""
figA_py.py — Paneles principales (Python) Paper2 V1, marco DUAL θ+β-CMC.
Genera PNG cuadrados individuales en outputs/figures/ desde .mat CANÓNICOS.
A1 RT, A2 IES, A3 masseter manip-check, A11 zMI θ doble null, A12 2x2 mecanismo,
A13 rose θ, A14 band comparison (θ/α/β), A15 rose β.
"""
import os, numpy as np, scipy.io as sio
import matplotlib.pyplot as plt
from _figstyle import (COL, LBL_CASE, LBL_CTRL, DIR_OUT, DIR_STATS, DIR_DATA,
                       new_square, finish, pstars, sig_bracket,
                       box_scatter, refline, symlog_zmi, robust_scatter)

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
    _paired_group(ax, BEH["rt_ctr_b1"], BEH["rt_ctr_b2"], COL["CTRL"], LBL_CTRL)
    top = _paired_group(ax, BEH["rt_cas_b1"], BEH["rt_cas_b2"], COL["CASE"], LBL_CASE)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["No-chew\n(B1)", "Chew\n(B2)"])
    ax.set_ylabel("Median RT (ms)"); ax.set_xlim(-0.35, 1.35)
    p = float(BEH["p_rt_sign"]); d = float(BEH["d_rt"])
    sig_bracket(ax, 0, 1, top*1.04, f"{pstars(p)}  (d={d:.2f})")
    ax.set_title("Reaction time", pad=8)
    ax.legend(loc="upper right", frameon=False, handlelength=1.2)
    finish(fig, ax, "behavior_RT_by_block_group")

# ── A2 — IES por bloque × grupo ──────────────────────────────────────────────
def A2():
    fig, ax = new_square()
    _paired_group(ax, BEH["ies_ctr_b1"], BEH["ies_ctr_b2"], COL["CTRL"], LBL_CTRL)
    top = _paired_group(ax, BEH["ies_cas_b1"], BEH["ies_cas_b2"], COL["CASE"], LBL_CASE)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["No-chew\n(B1)", "Chew\n(B2)"])
    ax.set_ylabel("IES (ms)"); ax.set_xlim(-0.35, 1.35)
    p = float(BEH["p_ies_sign"]); d = float(BEH["d_ies"])
    sig_bracket(ax, 0, 1, top*1.04, f"{pstars(p)}  (d={d:.2f})")
    ax.set_title("Inverse efficiency score", pad=8)
    ax.legend(loc="upper right", frameon=False, handlelength=1.2)
    finish(fig, ax, "behavior_IES_by_block_group")

# ── A3 — Manipulation check: masetero Ch−Nc por grupo ────────────────────────
def A3():
    cas, ctr = np.asarray(CHW["casd"], float), np.asarray(CHW["ctrd"], float)
    fig, ax = new_square()
    box_scatter(ax, [ctr, cas], [COL["CTRL"], COL["CASE"]],
                [f"Control\n(n={len(ctr)})", f"Cases\n(n={len(cas)})"], s=26)
    refline(ax, 0, style="--", color="k", lw=0.8, alpha=0.5)
    npos = int(np.sum(cas > 0))
    ax.set_ylabel("Masseter Ch−Nc (dB)")
    ax.set_title("Mastication manipulation check", pad=8)
    ax.text(1, np.max(cas)*1.02, f"+{np.mean(cas):.1f} dB · {npos}/{len(cas)} ***",
            ha="center", va="bottom", fontsize=9, color=COL["CASE"])
    ax.text(0, np.max(cas)*0.55, "n.s.", ha="center", color=COL["CTRL"], fontsize=9)
    finish(fig, ax, "manipcheck_masseter_SNR_ChNc_group")

def _zbox(ax, series, colors, labels, thr=1.96, symlog=True):
    """Box+scatter estilo figures_ok, con symlog para zMI de rango amplio y k/N por banda."""
    box_scatter(ax, series, colors, labels, s=26)
    allv = np.concatenate([np.asarray(v, float) for v in series])
    if symlog:
        symlog_zmi(ax, allv, thr=2)
    refline(ax, thr, label=f"z = {thr}", style="--", color=COL["SIG"], alpha=0.8)
    ytop = np.max(allv) * (1.05 if not symlog else 1.0)
    for j, dat in enumerate(series):
        dat = np.asarray(dat, float); k = int(np.sum(dat > thr))
        ax.annotate(f"{k}/{len(dat)}", (j, ytop), ha="center", va="bottom",
                    fontsize=8.5, color=colors[j], annotation_clip=False)

# ── A11 — zMI θ vs doble null (circ + AAFT) ──────────────────────────────────
def A11():
    fig, ax = new_square()
    s = [PAC["zC_ch"][:, TH], PAC["zA_ch"][:, TH]]
    _zbox(ax, s, [COL["THETA"], COL["THETA"]],
          ["circular-shift\nnull", "AAFT\nnull"])
    ax.set_ylabel("zMI θ (vs null)")
    ax.set_title("θ-PAC over double null", pad=8)
    ax.text(0.5, -0.02, "k = #subjects z>1.96", transform=ax.transAxes,
            ha="center", va="top", fontsize=8, color="0.4")
    finish(fig, ax, "pac_zMI_theta_doublenull")

# ── A12 — 2x2 mecanismo: Ch vs Nc + control negativo ─────────────────────────
def A12():
    fig, ax = new_square()
    s = [PAC["zC_ch"][:, TH], PAC["zC_nc"][:, TH], S8["zC"][:, TH]]
    _zbox(ax, s, [COL["CASE"], COL["NC"], COL["CTRL"]],
          ["Cases\nChew", "Cases\nNo-chew", "Controls\n\"Chew\""])
    ax.set_ylabel("zMI θ (circular null)")
    ax.set_title("θ-PAC driven by mastication", pad=8)
    top = np.percentile(PAC["zC_ch"][:, TH], 95)
    sig_bracket(ax, 0, 1, top*1.10, "***")
    finish(fig, ax, "pac_zMI_theta_chew_vs_nochew_negctrl")

# ── A14 — comparación de bandas (θ/α/β) doble null ───────────────────────────
def A14():
    fig, ax = new_square()
    s = [PAC["zC_ch"][:, TH], PAC["zC_ch"][:, AL], PAC["zC_ch"][:, BE]]
    _zbox(ax, s, [COL["THETA"], COL["ALPHA"], COL["BETA"]],
          ["θ\n4–7", "α\n8–13", "β\n13–30"])
    ax.set_ylabel("zMI (circular null)")
    ax.set_title("Broadband PAC — β strongest", pad=8)
    ax.text(0.5, 0.92, "Friedman p<0.0001", transform=ax.transAxes,
            ha="center", fontsize=8.5, color="0.35")
    finish(fig, ax, "pac_zMI_band_comparison_doublenull")

# ── A13 / A15 — Rayleigh rose ────────────────────────────────────────────────
def _rose(band_idx, color, title, name):
    ang = np.asarray(PAC["pref_ch"][:, band_idx], float)
    Z = float(np.asarray(PAC["rayl_Z"])[band_idx]); R = float(np.asarray(PAC["rayl_R"])[band_idx])
    fig = plt.figure(figsize=(3.5, 3.5))
    ax = fig.add_subplot(111, projection="polar")
    nb = 16; edges = np.linspace(-np.pi, np.pi, nb+1)
    cnt, _ = np.histogram(ang, bins=edges)
    width = np.diff(edges); centers = edges[:-1] + width/2
    ax.bar(centers, cnt, width=width, color=color, alpha=0.85,
           edgecolor="white", linewidth=0.6)
    mu = np.angle(np.mean(np.exp(1j*ang)))
    ax.annotate("", xy=(mu, np.max(cnt)*R + 0.1), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.6))
    ax.set_theta_zero_location("E"); ax.set_yticklabels([])
    ax.set_title(f"{title}\nZ={Z:.2f}, R={R:.2f}", fontsize=11, pad=12)
    ax.tick_params(labelsize=9)
    out = os.path.join(DIR_OUT, name + ".png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print(f"  saved {name}.png")

def A13(): _rose(TH, COL["THETA"], "θ preferred phase", "pac_rayleigh_theta_rose")
def A15(): _rose(BE, COL["BETA"],  "β preferred phase", "pac_rayleigh_beta_rose")

if __name__ == "__main__":
    print("== figA_py: paneles principales ==")
    for fn in (A1, A2, A3, A11, A12, A14, A13, A15):
        fn()
    print("Listo. PNGs en outputs/figures/")
