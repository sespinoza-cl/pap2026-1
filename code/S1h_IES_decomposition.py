"""
S1h_IES_decomposition.py
========================
IES Shift Decomposition: descompone ΔIES (Chew - NoChew) en
contribución de RT y contribución de ACC, por grupo.

IES = mean_RT / accuracy  (por sujeto)
RT contrib_i  = (RT_ch_i - RT_nc_i) / acc_nc_i
ACC contrib_i = RT_ch_i * (1/acc_ch_i - 1/acc_nc_i)
Verificación: RT + ACC = IES_ch_i - IES_nc_i

Salida: S1h_IES_shift_decomposition.png
"""

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os

ROOT    = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
BEH_RAW = os.path.join(ROOT, "Analysis_V1_Final", "data", "computed",
                        "data_beh_tb.mat")
DIR_OUT = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "figures")
os.makedirs(DIR_OUT, exist_ok=True)

with h5py.File(BEH_RAW, 'r') as f:
    rt_cas_nc  = np.array(f['tb_data/casos/nochew/mean']).flatten()
    rt_cas_ch  = np.array(f['tb_data/casos/chew/mean']).flatten()
    rt_ctr_nc  = np.array(f['tb_data/controles/nochew/mean']).flatten()
    rt_ctr_ch  = np.array(f['tb_data/controles/chew/mean']).flatten()
    acc_cas_nc = np.array(f['tb_data/casos/nochew/acc']).flatten()
    acc_cas_ch = np.array(f['tb_data/casos/chew/acc']).flatten()
    acc_ctr_nc = np.array(f['tb_data/controles/nochew/acc']).flatten()
    acc_ctr_ch = np.array(f['tb_data/controles/chew/acc']).flatten()
    ies_cas_nc = np.array(f['tb_data/casos/nochew/ies']).flatten()
    ies_cas_ch = np.array(f['tb_data/casos/chew/ies']).flatten()
    ies_ctr_nc = np.array(f['tb_data/controles/nochew/ies']).flatten()
    ies_ctr_ch = np.array(f['tb_data/controles/chew/ies']).flatten()

# ── Descomposición por sujeto ─────────────────────────────────────────────
# RT contribution (holding ACC at baseline)
rt_contrib_cas  = (rt_cas_ch - rt_cas_nc) / acc_cas_nc
acc_contrib_cas = rt_cas_ch * (1.0/acc_cas_ch - 1.0/acc_cas_nc)

rt_contrib_ctr  = (rt_ctr_ch - rt_ctr_nc) / acc_ctr_nc
acc_contrib_ctr = rt_ctr_ch * (1.0/acc_ctr_ch - 1.0/acc_ctr_nc)

# Como % de IES_nc (usando IES almacenado como denominador)
rt_pct_cas  = rt_contrib_cas  / ies_cas_nc * 100
acc_pct_cas = acc_contrib_cas / ies_cas_nc * 100
tot_pct_cas = (ies_cas_ch - ies_cas_nc) / ies_cas_nc * 100

rt_pct_ctr  = rt_contrib_ctr  / ies_ctr_nc * 100
acc_pct_ctr = acc_contrib_ctr / ies_ctr_nc * 100
tot_pct_ctr = (ies_ctr_ch - ies_ctr_nc) / ies_ctr_nc * 100

# ── Estadísticas (Wilcoxon vs 0) ─────────────────────────────────────────
_, p_cas_total = stats.wilcoxon(tot_pct_cas)
_, p_ctr_total = stats.wilcoxon(tot_pct_ctr)
_, p_cas_rt    = stats.wilcoxon(rt_pct_cas)
_, p_ctr_rt    = stats.wilcoxon(rt_pct_ctr)

print('=== IES Shift Decomposition ===')
print(f'Controls  RT: {np.mean(rt_pct_ctr):.1f}%  '
      f'ACC: {np.mean(acc_pct_ctr):.1f}%  '
      f'Total: {np.mean(tot_pct_ctr):.1f}%  (p={p_ctr_total:.3f})')
print(f'Cases     RT: {np.mean(rt_pct_cas):.1f}%  '
      f'ACC: {np.mean(acc_pct_cas):.1f}%  '
      f'Total: {np.mean(tot_pct_cas):.1f}%  (p={p_cas_total:.4f})')

# ── Medias de grupo ───────────────────────────────────────────────────────
rt_means  = [np.mean(rt_pct_ctr),  np.mean(rt_pct_cas)]
acc_means = [np.mean(acc_pct_ctr), np.mean(acc_pct_cas)]
tot_means = [np.mean(tot_pct_ctr), np.mean(tot_pct_cas)]

# ── Figura ────────────────────────────────────────────────────────────────
# Paleta Okabe-Ito: Controls=azul #0072B2, Cases=verde #009E73
# RT contrib = color saturado del grupo; ACC contrib = versión clara
COLS_RT  = ['#0072B2', '#009E73']   # [Controls, Cases]
COLS_ACC = ['#7BBBDC', '#70CCB0']   # versión clara  [Controls, Cases]
COL_BLK  = 'black'

fig, ax = plt.subplots(figsize=(5, 5.5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

x = np.array([0.0, 1.0])
w = 0.55

# Barras apiladas: cada grupo usa su propio color Okabe-Ito
for i, xi in enumerate(x):
    ax.bar(xi, rt_means[i],  width=w, color=COLS_RT[i],  zorder=3,
           label='RT Contribution'  if i == 0 else '_nolegend_')
    ax.bar(xi, acc_means[i], width=w, bottom=rt_means[i], color=COLS_ACC[i], zorder=3,
           label='ACC Contribution' if i == 0 else '_nolegend_')

# Diamante = total ΔIES
ax.scatter(x, tot_means, marker='D', s=90, color=COL_BLK, zorder=5,
           label='Total Δ IES (%)')

# Línea cero
ax.axhline(0, color='black', lw=1.0, zorder=2)

# Grid horizontal
ax.yaxis.grid(True, ls='--', lw=0.6, color='#BBBBBB', zorder=0)
ax.set_axisbelow(True)

# Etiquetas de significancia sobre/bajo las barras
def pstar(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'n.s.'

y_offset = 0.8  # % arriba del total
for i, (tot, p) in enumerate(zip(tot_means, [p_ctr_total, p_cas_total])):
    ax.text(x[i], tot - y_offset, pstar(p),
            ha='center', va='top', fontsize=13, fontweight='bold', color='black')

# Ejes y etiquetas
ax.set_xticks(x)
ax.set_xticklabels(['Controls', 'Cases'], fontsize=14, fontweight='bold')
ax.set_ylabel('Δ IES (%)  (Chew − NoChew)', fontsize=12, fontweight='bold')
ax.set_title('IES Shift Decomposition', fontsize=14, fontweight='bold')

# Límites Y — un poco de margen abajo del mínimo
ymin_all = min(min(rt_means[i] + acc_means[i] for i in range(2)),
               min(tot_means))
ax.set_ylim(ymin_all * 1.25, ax.get_ylim()[1])

ax.legend(frameon=False, fontsize=10, loc='lower left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(pad=0.4)
fname = 'S1h_IES_shift_decomposition.png'
fig.savefig(os.path.join(DIR_OUT, fname),
            dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Guardada: {fname}')
