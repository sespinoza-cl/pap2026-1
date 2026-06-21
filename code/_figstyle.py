"""
_figstyle.py — Estándar único de figuras Paper2 V1 (marco DUAL θ+β-CMC).
Cuadradas 3.5x3.5 in · SOLO PNG @300 dpi · Arial título 12 / label 11 / ticks 9 / leyenda 9.
Paleta Okabe-Ito (= S0_config.m). Etiquetas SIN TMD. Salida -> outputs/figures/.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ── Okabe-Ito (idéntico a S0_config.m) ───────────────────────────────────────
COL = {
    "CASE":  "#009E73",  # bluish-green
    "CTRL":  "#0072B2",  # blue
    "NC":    "#888888",  # neutral grey
    "THETA": "#D55E00",  # vermillion
    "ALPHA": "#56B4E9",  # sky blue
    "BETA":  "#CC79A7",  # reddish-purple
    "SIG":   "#E63232",  # red marker
}
LBL_CASE = "Chewing group\n(Cases, n=31)"
LBL_CTRL = "Control / no-chew\n(n=15)"

# ── Tipografía homogénea ─────────────────────────────────────────────────────
_FS_TITLE, _FS_LABEL, _FS_TICK, _FS_LEG = 12, 11, 9, 9
def _pick_font():
    for fam in ("Arial", "Liberation Sans", "DejaVu Sans"):
        try:
            font_manager.findfont(fam, fallback_to_default=False)
            return fam
        except Exception:
            continue
    return "DejaVu Sans"
_FAM = _pick_font()
plt.rcParams.update({
    "font.family": _FAM, "font.size": _FS_TICK,
    "axes.titlesize": _FS_TITLE, "axes.labelsize": _FS_LABEL,
    "xtick.labelsize": _FS_TICK, "ytick.labelsize": _FS_TICK,
    "legend.fontsize": _FS_LEG, "axes.linewidth": 0.9,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "savefig.dpi": 300, "figure.dpi": 300,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
if _FAM != "Arial":
    print(f"[_figstyle] WARNING: Arial no encontrada, usando '{_FAM}'.")

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT  = os.path.dirname(_HERE)
DIR_OUT   = os.path.join(ROOT, "outputs", "figures")
DIR_STATS = os.path.join(ROOT, "outputs", "stats")
DIR_DATA  = os.path.join(ROOT, "data", "computed")
os.makedirs(DIR_OUT, exist_ok=True)

def new_square():
    """Figura cuadrada estándar 3.5x3.5 in."""
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    return fig, ax

def finish(fig, ax, name):
    """Aplica límites/estilo común y guarda SOLO PNG @300 dpi."""
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.6)
    out = os.path.join(DIR_OUT, name if name.endswith(".png") else name + ".png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved {os.path.basename(out)}")
    return out

def pstars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."

def sig_bracket(ax, x1, x2, y, label, lw=1.0):
    ax.plot([x1, x1, x2, x2], [y, y*1.0, y*1.0, y], lw=lw, c="k", clip_on=False)
    ax.text((x1+x2)/2, y, label, ha="center", va="bottom", fontsize=_FS_TICK)

# ── Estética figures_ok: box+scatter, refline, regresión robusta ─────────────
import matplotlib.ticker as _mticker
from matplotlib.patches import Rectangle as _Rect

def box_scatter(ax, series, colors, labels, bw=0.30, jitter=0.10, s=30, seed=42):
    """Box(IQR, mediana, bigotes) + scatter jittered translucido — estilo figures_ok.
    series: lista de arrays. Devuelve posiciones x (0..n-1)."""
    rng = np.random.default_rng(seed)
    allv = np.concatenate([np.asarray(v, float) for v in series])
    for j, (v, col) in enumerate(zip(series, colors)):
        v = np.asarray(v, float)
        q1, q2, q3 = np.percentile(v, [25, 50, 75]); iqr = q3 - q1
        wlo = v[v >= q1 - 1.5*iqr].min() if np.any(v >= q1-1.5*iqr) else v.min()
        whi = v[v <= q3 + 1.5*iqr].max() if np.any(v <= q3+1.5*iqr) else v.max()
        ax.add_patch(_Rect((j-bw/2, q1), bw, q3-q1, facecolor=col, alpha=0.32,
                           edgecolor=col, lw=1.5, zorder=2))
        ax.plot([j-bw/2, j+bw/2], [q2, q2], color=col, lw=2.4, zorder=4)
        ax.plot([j, j], [wlo, q1], color=col, lw=1.1, zorder=2)
        ax.plot([j, j], [q3, whi], color=col, lw=1.1, zorder=2)
        jit = (rng.random(len(v))-0.5)*2*jitter
        ax.scatter(j+jit, v, s=s, color=col, alpha=0.50, linewidths=0, zorder=3)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
    ax.set_xlim(-0.6, len(labels)-0.4)
    return np.arange(len(labels))

def refline(ax, y, label=None, style=":", color="gray", lw=0.9, alpha=1.0, lx=None):
    ax.axhline(y, color=color, ls=style, lw=lw, alpha=alpha, zorder=1)
    if label:
        x = lx if lx is not None else ax.get_xlim()[1]
        ax.text(x, y, label, va="center", ha="right", fontsize=_FS_TICK-1, color=color)

def symlog_zmi(ax, vals, thr=2):
    ax.set_yscale("symlog", linthresh=thr)
    ax.set_yticks([0, 2, 5, 10, 20, 50, 100])
    ax.yaxis.set_major_formatter(_mticker.ScalarFormatter())
    ax.set_ylim(min(np.min(vals), -1), np.max(vals)*1.18)

def robust_scatter(ax, x, y, color, xlabel, ylabel, title=None, s=34):
    """Scatter + recta robusta (Theil-Sen) + anotación rho/p (Spearman) — estilo figures_ok."""
    from scipy.stats import spearmanr, theilslopes
    x = np.asarray(x, float); y = np.asarray(y, float)
    ax.scatter(x, y, s=s, color=color, alpha=0.55, linewidths=0, zorder=3)
    res = theilslopes(y, x); xl = np.array([x.min(), x.max()])
    ax.plot(xl, res[1] + res[0]*xl, color=color, lw=1.8, alpha=0.85, zorder=2)
    rho, p = spearmanr(x, y)
    pstr = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    star = "  *" if p < 0.05 else ""
    ax.text(0.96, 0.05, f"$\\rho$ = {rho:+.2f}\n{pstr}{star}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=_FS_TICK,
            color=color if p < 0.05 else "#888888",
            fontweight="bold" if p < 0.05 else "normal")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title: ax.set_title(title, pad=8)
    return rho, p
