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
LBL_CASE = "Cases"
LBL_CTRL = "Controls"

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
        # etiqueta POR ENCIMA de la línea (no la cruza)
        ax.text(x, y, label, va="bottom", ha="right", fontsize=_FS_TICK-1, color=color)

def rose_v2(phi, R, color, band_lbl, stars, name, nb=12):
    """Rose estilo v2: sectores rellenos + círculos de referencia graduados con etiqueta
    de escala radial (n sujetos por bin) + vector resultante. Cuadrado, PNG@300."""
    phi = np.asarray(phi, float); phi = phi[~np.isnan(phi)]
    fig = plt.figure(figsize=(3.5, 3.5)); ax = fig.add_axes([0.02, 0.02, 0.96, 0.92])
    ax.set_aspect("equal"); ax.axis("off")
    edges = np.linspace(-np.pi, np.pi, nb+1)
    counts, _ = np.histogram(phi, edges)
    n_outer = int(counts.max())
    # Round reference values for scale circles
    n_mid = max(1, round(n_outer / 2))
    rmax = n_outer * 1.30 + 0.5
    tt = np.linspace(0, 2*np.pi, 120)
    # Draw sectors
    for b in range(nb):
        if counts[b] == 0: continue
        th = np.linspace(edges[b], edges[b+1], 20)
        xs = np.concatenate([[0], counts[b]*np.sin(th), [0]])
        ys = np.concatenate([[0], counts[b]*np.cos(th), [0]])
        ax.fill(xs, ys, color=color, alpha=0.85, edgecolor="white", lw=0.5, zorder=2)
    # Outer frame circle (gray)
    ax.plot(rmax*np.sin(tt), rmax*np.cos(tt), "-", color="0.75", lw=0.8, zorder=1)
    # Radial scale: mid circle (dashed) + max circle (dashed)
    for nr, ls in [(n_mid, ":"), (n_outer, "--")]:
        ax.plot(nr*np.sin(tt), nr*np.cos(tt), ls, color="0.65", lw=0.7, zorder=1)
        # Label at ~45° (upper-right quadrant)
        ang45 = np.pi / 4
        ax.text(nr*np.sin(ang45)*1.06, nr*np.cos(ang45)*1.06,
                str(nr), ha="left", va="bottom", fontsize=7.5, color="0.40")
    # Cross lines
    ax.plot([-rmax, rmax], [0, 0], "-", color="0.85", lw=0.8, zorder=0)
    ax.plot([0, 0], [-rmax, rmax], "-", color="0.85", lw=0.8, zorder=0)
    # Mean resultant vector
    mu = np.angle(np.mean(np.exp(1j*phi)))
    ax.annotate("", xy=(rmax*0.88*R*np.sin(mu), rmax*0.88*R*np.cos(mu)), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=2.4), zorder=4)
    # Cardinal labels and radius unit note
    ax.text(0, rmax*1.08, "0", ha="center", va="bottom", fontsize=9)
    ax.text(0, -rmax*1.08, "$\\pi$", ha="center", va="top", fontsize=10)
    ax.text(0, -rmax*1.28, "radius = subjects per bin", ha="center", va="top",
            fontsize=7, color="0.50")
    ax.set_xlim(-rmax*1.20, rmax*1.20); ax.set_ylim(-rmax*1.35, rmax*1.15)
    ax.set_title(f"{band_lbl}   R={R:.2f}  {stars}", fontsize=_FS_TITLE, fontweight="normal")
    out = os.path.join(DIR_OUT, name+".png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print(f"  saved {name}.png")

def symlog_zmi(ax, vals, thr=2):
    ax.set_yscale("symlog", linthresh=thr)
    ax.set_yticks([0, 2, 5, 10, 20, 50, 100])
    ax.yaxis.set_major_formatter(_mticker.ScalarFormatter())
    ax.set_ylim(min(np.min(vals), -1), np.max(vals)*1.18)

def paired_box(ax, groups, ylabel, title=None, bw=0.32, gap=1.0):
    """Boxplot (IQR+mediana+bigotes) por condición con puntos UNIDOS (pareados) y
    significancia SOLO con asteriscos. groups: lista de dict(label,color,nc,ch,stars).
    Sin n ni p en el plot."""
    pos, labs = [], []
    x = 0.0; ymax = -np.inf
    for g in groups:
        nc, ch = np.asarray(g["nc"], float), np.asarray(g["ch"], float)
        xp = [x, x+1]; ymax = max(ymax, nc.max(), ch.max())
        for i in range(len(nc)):                      # líneas pareadas (sin jitter)
            ax.plot(xp, [nc[i], ch[i]], color=g["color"], alpha=0.10, lw=0.6, zorder=1)
        for xi, dat in zip(xp, [nc, ch]):
            q1, q2, q3 = np.percentile(dat, [25, 50, 75]); iqr = q3-q1
            wlo = dat[dat >= q1-1.5*iqr].min(); whi = dat[dat <= q3+1.5*iqr].max()
            ax.add_patch(_Rect((xi-bw/2, q1), bw, q3-q1, facecolor=g["color"], alpha=0.30,
                               edgecolor=g["color"], lw=1.4, zorder=2))
            ax.plot([xi-bw/2, xi+bw/2], [q2, q2], color=g["color"], lw=2.3, zorder=4)
            ax.plot([xi, xi], [wlo, q1], color=g["color"], lw=1.0, zorder=2)
            ax.plot([xi, xi], [q3, whi], color=g["color"], lw=1.0, zorder=2)
            ax.scatter([xi]*len(dat), dat, s=14, color=g["color"], alpha=0.45,
                       linewidths=0, zorder=3)
        ax.text(x+0.5, 1.0, g["label"], transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=_FS_TICK, color=g["color"])
        pos += xp; labs += ["No-chew", "Chew"]
        if g.get("stars"):
            yb = max(nc.max(), ch.max())
            sig_bracket(ax, xp[0], xp[1], yb*1.03 if yb > 0 else yb+0.03*abs(yb), g["stars"])
        x += 2 + gap
    ax.set_xticks(pos); ax.set_xticklabels(labs)
    ax.set_xlim(-0.6, pos[-1]+0.6); ax.set_ylabel(ylabel)
    if title: ax.set_title(title, pad=10)

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
