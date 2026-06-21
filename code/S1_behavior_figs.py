"""
S1_behavior_figs.py
===================
Figuras conductuales individuales para montaje en Inkscape.

Salidas (Analysis_V1_Final/outputs/figures/):
  S1a_RT_baseline.png        — RT bloque 1: Controls vs Cases (equivalencia, p=0.981)
  S1b_IES_Cases.png          — IES NoChew vs Chew, Cases (paired, ***)
  S1c_IES_Controls.png       — IES NoChew vs Chew, Controls (paired, NS)
  S1d_RT_Cases.png           — RT NoChew vs Chew, Cases (paired, ***)
  S1e_RT_Controls.png        — RT NoChew vs Chew, Controls (paired, NS)
  S1f_IES_delta.png          — ΔIES (Chew-NoChew) Controls vs Cases + sig vs 0
  S1g_RT_delta.png           — ΔRT  (Chew-NoChew) Controls vs Cases + sig vs 0
"""

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import statsmodels.formula.api as smf
import os

ROOT    = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
BEH_RAW = os.path.join(ROOT, "Analysis_V1_Final", "data", "computed",
                        "data_beh_tb.mat")
DIR_OUT = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "figures")
os.makedirs(DIR_OUT, exist_ok=True)

# ── Cargar datos desde fuente canónica ───────────────────────────────────
with h5py.File(BEH_RAW, 'r') as f:
    rt_cas_b1  = np.array(f['tb_data/casos/nochew/med']).flatten()
    rt_cas_b2  = np.array(f['tb_data/casos/chew/med']).flatten()
    rt_ctr_b1  = np.array(f['tb_data/controles/nochew/med']).flatten()
    rt_ctr_b2  = np.array(f['tb_data/controles/chew/med']).flatten()
    ies_cas_b1 = np.array(f['tb_data/casos/nochew/ies']).flatten()
    ies_cas_b2 = np.array(f['tb_data/casos/chew/ies']).flatten()
    ies_ctr_b1 = np.array(f['tb_data/controles/nochew/ies']).flatten()
    ies_ctr_b2 = np.array(f['tb_data/controles/chew/ies']).flatten()

# ── Estadísticas ─────────────────────────────────────────────────────────
# Equivalencia baseline (B1)
_, p_rt_b1  = stats.mannwhitneyu(rt_cas_b1,  rt_ctr_b1,  alternative='two-sided')
_, p_ies_b1 = stats.mannwhitneyu(ies_cas_b1, ies_ctr_b1, alternative='two-sided')

# Efecto within (Wilcoxon signed-rank vs 0)
_, p_rt_sign  = stats.wilcoxon(rt_cas_b2  - rt_cas_b1)
_, p_ies_sign = stats.wilcoxon(ies_cas_b2 - ies_cas_b1)
_, p_rt_ctr   = stats.wilcoxon(rt_ctr_b2  - rt_ctr_b1)
_, p_ies_ctr  = stats.wilcoxon(ies_ctr_b2 - ies_ctr_b1)

# LME Group × Block (equivalente al fitlme de MATLAB)
n_cas, n_ctr = len(rt_cas_b1), len(rt_ctr_b1)
subj = (list(range(n_cas)) * 2 + list(range(n_cas, n_cas + n_ctr)) * 2)
grp  = (['Cases'] * n_cas + ['Cases'] * n_cas +
        ['Controls'] * n_ctr + ['Controls'] * n_ctr)
blk  = ([1] * n_cas + [2] * n_cas + [1] * n_ctr + [2] * n_ctr)

df_rt  = pd.DataFrame({'RT':  np.concatenate([rt_cas_b1,  rt_cas_b2,  rt_ctr_b1,  rt_ctr_b2]),
                        'IES': np.concatenate([ies_cas_b1, ies_cas_b2, ies_ctr_b1, ies_ctr_b2]),
                        'Group': grp, 'Block': blk, 'Subject': subj})
df_rt['Block']   = df_rt['Block'].astype('category')
df_rt['Group']   = df_rt['Group'].astype('category')

lme_rt  = smf.mixedlm('RT  ~ C(Group) * C(Block)', df_rt,
                        groups=df_rt['Subject']).fit(reml=False)
lme_ies = smf.mixedlm('IES ~ C(Group) * C(Block)', df_rt,
                        groups=df_rt['Subject']).fit(reml=False)

# p del término de interacción
int_key_rt  = [k for k in lme_rt.pvalues.index  if 'Group' in k and 'Block' in k][0]
int_key_ies = [k for k in lme_ies.pvalues.index if 'Group' in k and 'Block' in k][0]
p_lme_int_rt  = lme_rt.pvalues[int_key_rt]
p_lme_int_ies = lme_ies.pvalues[int_key_ies]

print(f'\n=== Stats desde data_beh_tb.mat ===')
print(f'Baseline RT  Cases vs Controls: p = {p_rt_b1:.3f}')
print(f'Baseline IES Cases vs Controls: p = {p_ies_b1:.3f}')
print(f'Wilcoxon IES Cases (vs 0): p = {p_ies_sign:.4f}')
print(f'Wilcoxon RT  Cases (vs 0): p = {p_rt_sign:.4f}')
print(f'Wilcoxon IES Controls (vs 0): p = {p_ies_ctr:.4f}')
print(f'Wilcoxon RT  Controls (vs 0): p = {p_rt_ctr:.4f}')
print(f'LME Group×Block RT:  p = {p_lme_int_rt:.4f}')
print(f'LME Group×Block IES: p = {p_lme_int_ies:.4f}')

COL_CAS = '#009E73'   # Okabe-Ito bluish-green  (Cases)
COL_CTR = '#0072B2'   # Okabe-Ito blue         (Controls)
COL_NC  = '#888888'   # gris neutro             (NoChew)


def pstar(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'n.s.'


def fmt_p(p):
    if p < 0.001: return 'p < 0.001'
    return f'p = {p:.3f}'


# ═══════════════════════════════════════════════════════════════════════════
# FIG A — RT baseline (B1): Controls vs Cases  [raincloud style]
# ═══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(4, 5))
fig.patch.set_facecolor('white')

rng = np.random.default_rng(42)
groups = [rt_ctr_b1, rt_cas_b1]
cols   = [COL_CTR, COL_CAS]
xlabs  = ['Controls', 'Cases']
xpos   = [1, 2]

for i, (data, col, lab) in enumerate(zip(groups, cols, xlabs)):
    x = xpos[i]
    # Jitter
    jit = (rng.random(len(data)) - 0.5) * 0.22
    ax.scatter(x - 0.3 + jit, data, s=28, color=col, alpha=0.55,
               linewidths=0, zorder=3)
    # Boxplot manual
    q1, q2, q3 = np.percentile(data, [25, 50, 75])
    iqr = q3 - q1
    wlo = data[data >= q1 - 1.5*iqr].min()
    whi = data[data <= q3 + 1.5*iqr].max()
    bw = 0.22
    rect = plt.Rectangle((x - bw/2, q1), bw, q3-q1,
                          facecolor=col, alpha=0.40, edgecolor=col,
                          lw=1.5, zorder=2)
    ax.add_patch(rect)
    ax.plot([x - bw/2, x + bw/2], [q2, q2], color='white', lw=2, zorder=4)
    ax.plot([x, x], [wlo, q1], color=col, lw=1.1, zorder=2)
    ax.plot([x, x], [q3, whi], color=col, lw=1.1, zorder=2)
    # Mean triangle
    ax.scatter(x, np.mean(data), marker='^', s=55, color='black',
               zorder=5, linewidths=0)

# Significance bracket
ymax = max(rt_ctr_b1.max(), rt_cas_b1.max())
ytop = ymax * 1.08
ax.plot([1, 1, 2, 2], [ytop, ytop*1.02, ytop*1.02, ytop],
        color='black', lw=1.2)
ps = pstar(p_rt_b1)
ax.text(1.5, ytop * 1.03, f'{fmt_p(p_rt_b1)}', ha='center', va='bottom',
        fontsize=9, color='black')

ax.set_xlim(0.5, 2.8)
ax.set_xticks(xpos)
ax.set_xticklabels(xlabs, fontsize=11)
ax.set_ylabel('Baseline RT  (ms)', fontsize=11)
ax.set_ylim(250, ytop * 1.12)
ax.tick_params(direction='out', labelsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=0.4)
fname = 'S1a_RT_baseline.png'
fig.savefig(os.path.join(DIR_OUT, fname),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Guardada: {fname}  |  p={p_rt_b1:.3f}')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURAS B-E — Paired boxplot NoChew vs Chew (estilo FOOOF exponent)
# ═══════════════════════════════════════════════════════════════════════════
def paired_boxplot(data_nc, data_ch, col, p_val, ylabel, fname):
    fig, ax = plt.subplots(figsize=(3.5, 5))
    fig.patch.set_facecolor('white')

    xpos = [1, 2]
    xlabs = ['No Chew', 'Chew']
    col_nc = '#888888'

    # Líneas individuales
    for i in range(len(data_nc)):
        ax.plot([1, 2], [data_nc[i], data_ch[i]],
                color='gray', lw=0.7, alpha=0.45, zorder=1)

    for j, (data, col_box, x) in enumerate(
            zip([data_nc, data_ch], [col_nc, col], xpos)):
        rng_l = np.random.default_rng(j + 10)
        jit = (rng_l.random(len(data)) - 0.5) * 0.18
        ax.scatter(x + jit, data, s=26, color=col_box, alpha=0.55,
                   linewidths=0, zorder=3)
        q1, q2, q3 = np.percentile(data, [25, 50, 75])
        iqr = q3 - q1
        wlo = data[data >= q1 - 1.5*iqr].min()
        whi = data[data <= q3 + 1.5*iqr].max()
        bw = 0.26
        rect = plt.Rectangle((x - bw/2, q1), bw, q3-q1,
                              facecolor=col_box, alpha=0.38,
                              edgecolor=col_box, lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.plot([x - bw/2, x + bw/2], [q2, q2],
                color='white', lw=2.2, zorder=4)
        ax.plot([x, x], [wlo, q1], color=col_box, lw=1.1, zorder=2)
        ax.plot([x, x], [q3, whi], color=col_box, lw=1.1, zorder=2)
        ax.scatter(x, np.mean(data), marker='^', s=55, color='black',
                   zorder=5, linewidths=0)

    # Bracket significancia
    yall = np.concatenate([data_nc, data_ch])
    ytop = yall.max() * 1.10
    ax.plot([1, 1, 2, 2], [ytop, ytop*1.02, ytop*1.02, ytop],
            color='black', lw=1.2)
    star = pstar(p_val)
    txt = star if p_val < 0.05 else fmt_p(p_val)
    ax.text(1.5, ytop * 1.035, txt, ha='center', va='bottom',
            fontsize=11 if p_val < 0.05 else 9,
            fontweight='bold' if p_val < 0.05 else 'normal')

    ax.set_xlim(0.55, 2.55)
    ax.set_xticks(xpos)
    ax.set_xticklabels(xlabs, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_ylim(max(0, yall.min() * 0.88), ytop * 1.13)
    ax.tick_params(direction='out', labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout(pad=0.4)
    fig.savefig(os.path.join(DIR_OUT, fname),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Guardada: {fname}  |  p={p_val:.4f}  ({pstar(p_val)})')


paired_boxplot(ies_cas_b1, ies_cas_b2, COL_CAS, p_ies_sign,
               'IES  (ms)', 'S1b_IES_Cases.png')

paired_boxplot(ies_ctr_b1, ies_ctr_b2, COL_CTR, p_ies_ctr,
               'IES  (ms)', 'S1c_IES_Controls.png')

paired_boxplot(rt_cas_b1, rt_cas_b2, COL_CAS, p_rt_sign,
               'Median RT  (ms)', 'S1d_RT_Cases.png')

paired_boxplot(rt_ctr_b1, rt_ctr_b2, COL_CTR, p_rt_ctr,
               'Median RT  (ms)', 'S1e_RT_Controls.png')


# ═══════════════════════════════════════════════════════════════════════════
# FIGURAS F-G — Delta (Chew - NoChew): Controls vs Cases, sig vs 0
# ═══════════════════════════════════════════════════════════════════════════
ies_delta_cas = ies_cas_b2 - ies_cas_b1
ies_delta_ctr = ies_ctr_b2 - ies_ctr_b1
rt_delta_cas  = rt_cas_b2  - rt_cas_b1
rt_delta_ctr  = rt_ctr_b2  - rt_ctr_b1


def delta_figure(delta_ctr, delta_cas, p_ctr_vs0, p_cas_vs0, p_between,
                 ylabel, fname):
    fig, ax = plt.subplots(figsize=(4, 5))
    fig.patch.set_facecolor('white')

    datasets = [delta_ctr, delta_cas]
    cols     = [COL_CTR, COL_CAS]
    xlabs    = ['Controls', 'Cases']
    xpos     = [1, 2]
    p_vs0    = [p_ctr_vs0, p_cas_vs0]

    rng2 = np.random.default_rng(99)
    for i, (data, col, p0) in enumerate(zip(datasets, cols, p_vs0)):
        x = xpos[i]
        jit = (rng2.random(len(data)) - 0.5) * 0.22
        ax.scatter(x + jit, data, s=28, color=col, alpha=0.55,
                   linewidths=0, zorder=3)
        # Boxplot manual
        q1, q2, q3 = np.percentile(data, [25, 50, 75])
        iqr = q3 - q1
        wlo = data[data >= q1 - 1.5*iqr].min()
        whi = data[data <= q3 + 1.5*iqr].max()
        bw = 0.28
        rect = plt.Rectangle((x - bw/2, q1), bw, q3-q1,
                              facecolor=col, alpha=0.38, edgecolor=col,
                              lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.plot([x - bw/2, x + bw/2], [q2, q2],
                color='white', lw=2.2, zorder=4)
        ax.plot([x, x], [wlo, q1], color=col, lw=1.1, zorder=2)
        ax.plot([x, x], [q3, whi], color=col, lw=1.1, zorder=2)
        ax.scatter(x, np.mean(data), marker='^', s=55, color='black',
                   zorder=5, linewidths=0)

        # Marca de sig vs 0 encima del whisker superior
        star = pstar(p0)
        y_ann = whi + (whi - wlo) * 0.06
        ax.text(x, y_ann, star, ha='center', va='bottom',
                fontsize=12 if p0 < 0.05 else 9,
                color=col if p0 < 0.05 else '#888888',
                fontweight='bold' if p0 < 0.05 else 'normal')

    # Línea en 0
    ax.axhline(0, color='black', ls='--', lw=0.9, alpha=0.6)

    # Bracket entre grupos (significancia del delta entre grupos)
    all_data = np.concatenate([delta_ctr, delta_cas])
    ymin_all = all_data.min()
    ymax_all = all_data.max()
    span = ymax_all - ymin_all
    ytop = ymax_all + span * 0.28
    ax.plot([1, 1, 2, 2], [ytop, ytop + span*0.03,
                            ytop + span*0.03, ytop],
            color='black', lw=1.2)
    ps_b = pstar(p_between)
    txt_b = ps_b if p_between < 0.05 else fmt_p(p_between)
    ax.text(1.5, ytop + span*0.04, txt_b, ha='center', va='bottom',
            fontsize=10 if p_between < 0.05 else 8,
            fontweight='bold' if p_between < 0.05 else 'normal')

    ax.set_xlim(0.5, 2.8)
    ax.set_xticks(xpos)
    ax.set_xticklabels(xlabs, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.tick_params(direction='out', labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout(pad=0.4)
    fig.savefig(os.path.join(DIR_OUT, fname),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Guardada: {fname}  |  Cases p={p_cas_vs0:.4f}  Ctrl p={p_ctr_vs0:.4f}  '
          f'between p={p_between:.4f}')


delta_figure(ies_delta_ctr, ies_delta_cas,
             p_ies_ctr, p_ies_sign, p_lme_int_ies,
             'ΔIES  Chew − No Chew  (ms)', 'S1f_IES_delta.png')

delta_figure(rt_delta_ctr, rt_delta_cas,
             p_rt_ctr, p_rt_sign, p_lme_int_rt,
             'ΔRT  Chew − No Chew  (ms)', 'S1g_RT_delta.png')
