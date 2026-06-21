"""
rev_dissociation2.py — Enfoque "elegante v1": frecuencia-PICO por banda en el ROI +
amplitud MEDIA/MEDIANA del ROI por banda, correlacionadas con la conducta (intra-casos, n=31).
Fuente: FOOOF_Workspace_V1.mat GR.Cases.PSD_Ch/PSD_Nc (espectro ROI por sujeto) + behavior deltas.
"""
import os, numpy as np, scipy.io as sio, h5py
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
ST = os.path.join(ROOT, "outputs", "stats"); DC = os.path.join(ROOT, "data", "computed")

beh = sio.loadmat(os.path.join(DC, "v1_S1_behavior_stats.mat"), squeeze_me=True)
drt = np.asarray(beh["rt_delta_cas"], float); dies = np.asarray(beh["ies_delta_cas"], float)

with h5py.File(os.path.join(ST, "FOOOF_Workspace_V1.mat"), "r") as h:
    f = np.asarray(h["GR"]["Cases"]["f"]).ravel()
    psd_ch = np.asarray(h["GR"]["Cases"]["PSD_Ch"])   # (31,193) o (193,31)
    psd_nc = np.asarray(h["GR"]["Cases"]["PSD_Nc"])
if psd_ch.shape[1] != len(f): psd_ch = psd_ch.T; psd_nc = psd_nc.T  # -> (subj, freq)
n = psd_ch.shape[0]
bands = {"theta": (4, 7), "alpha": (8, 13), "beta": (13, 30)}

def band_feats(psd, lo, hi):
    m = (f >= lo) & (f <= hi); fb = f[m]
    pk = np.array([fb[np.argmax(p[m])] for p in psd])   # peak freq
    mn = psd[:, m].mean(1); md = np.median(psd[:, m], 1)  # mean/median amp
    return pk, mn, md

def corr(x, y, lbl, out):
    rho, p = spearmanr(x, y)
    sig = "*" if p < 0.05 else " "
    out.append((lbl, rho, p, sig))
    return p

print(f"n={n}  ΔRT={drt.mean():+.1f}  ΔIES={dies.mean():+.1f}  (mejora = negativo)")
results = []
for b, (lo, hi) in bands.items():
    pk_c, mn_c, md_c = band_feats(psd_ch, lo, hi)
    pk_n, mn_n, md_n = band_feats(psd_nc, lo, hi)
    feats = {
        f"{b}_peakHz_chew": pk_c, f"{b}_dpeakHz": pk_c - pk_n,
        f"{b}_meanAmp_chew": mn_c, f"{b}_dmeanAmp": mn_c - mn_n,
        f"{b}_medAmp_chew": md_c, f"{b}_dmedAmp": md_c - md_n,
    }
    for name, val in feats.items():
        corr(val, drt, f"{name:22} × ΔRT", results)
        corr(val, dies, f"{name:22} × ΔIES", results)

# imprimir, marcar los nominalmente sig y aplicar Bonferroni
ps = [r[2] for r in results]
mtot = len(ps); bonf = 0.05 / mtot
print(f"\n=== Correlaciones neural(pico/amp por banda) × conducta — {mtot} tests, Bonferroni α={bonf:.4f} ===")
for lbl, rho, p, sig in results:
    flag = "  <-- nominal" if p < 0.05 else ""
    surv = "  [SOBREVIVE Bonferroni]" if p < bonf else ""
    print(f"  {lbl}  rho={rho:+.3f} p={p:.3f}{flag}{surv}")
nsig = sum(1 for _, _, p, _ in results if p < 0.05)
nsurv = sum(1 for _, _, p, _ in results if p < bonf)
print(f"\nnominal p<0.05: {nsig}/{mtot} (esperados por azar ~{0.05*mtot:.1f}) · sobreviven Bonferroni: {nsurv}")
