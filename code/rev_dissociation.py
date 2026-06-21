"""
rev_dissociation.py — Test disociativo θ/aperiódico ↔ conducta (intra-casos, n=31).
Responde al issue CRITICO del panel (DA-1): ¿algún cambio neural predice la mejora conductual
a nivel individual, aunque el PAC no lo haga?
"""
import os, numpy as np, scipy.io as sio, h5py
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
ST = os.path.join(ROOT, "outputs", "stats"); DC = os.path.join(ROOT, "data", "computed")

# ---- conducta (deltas chew-nochew por caso) ----
beh = sio.loadmat(os.path.join(DC, "v1_S1_behavior_stats.mat"), squeeze_me=True)
drt = np.asarray(beh["rt_delta_cas"], float)    # 31
dies = np.asarray(beh["ies_delta_cas"], float)  # 31

# ---- theta-late power por canal/sujeto (Morlet NPL dB, WIN_LATE) ----
def h5arr(path, *keys):
    out = {}
    with h5py.File(path, "r") as h:
        for k in keys:
            out[k] = np.asarray(h[k])
    return out
T = h5arr(os.path.join(DC, "v1_S2_TF_data.mat"), "theta_topo_ch", "theta_topo_nc")
def orient(X, n=64):
    return X if X.shape[0] == n else X.T
tch = orient(T["theta_topo_ch"]); tnc = orient(T["theta_topo_nc"])  # 64 x 31

# ROI-18 indices
R = h5arr(os.path.join(ST, "ROI_canonical.mat"), "ROI_IDX")
roi = np.asarray(R["ROI_IDX"]).ravel().astype(int) - 1
dtheta_late = (tch - tnc)[roi, :].mean(0)        # 31, Δθ-late ROI

# ---- aperiodic exponent por caso (FOOOF GR) ----
with h5py.File(os.path.join(ST, "FOOOF_Workspace_V1.mat"), "r") as h:
    exp_ch = np.asarray(h["GR"]["Cases"]["exp_Ch"]).ravel()
    exp_nc = np.asarray(h["GR"]["Cases"]["exp_Nc"]).ravel()
dexp = exp_ch - exp_nc                            # 31, negativo = aplanamiento

# ---- PAC zMI theta (para comparar) ----
pac = sio.loadmat(os.path.join(ST, "v1_S4b_PAC_ROI.mat"), squeeze_me=True)
zth = np.asarray(pac["zC_ch"], float)[:, 0]

print(f"n: drt={len(drt)} dies={len(dies)} dtheta={len(dtheta_late)} dexp={len(dexp)} zth={len(zth)}")
print(f"means: Δθ_late={dtheta_late.mean():+.3f} dB · Δexp={dexp.mean():+.3f} (flattening if <0) · "
      f"ΔRT={drt.mean():+.1f} ms · ΔIES={dies.mean():+.1f}")

def corr(x, y, lbl):
    n = min(len(x), len(y))
    rho, p = spearmanr(x[:n], y[:n])
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
    print(f"  {lbl:34} rho={rho:+.3f}  p={p:.3f}  {sig}")
    return rho, p

print("\n=== DISOCIATIVO: cambio neural × mejora conductual (intra-casos, Spearman) ===")
print("(mejora = ΔRT/ΔIES negativos; aplanamiento = Δexp negativo)")
corr(dtheta_late, drt, "Δθ-late  × ΔRT")
corr(dtheta_late, dies, "Δθ-late  × ΔIES")
corr(dexp, drt, "Δexponent × ΔRT")
corr(dexp, dies, "Δexponent × ΔIES")
print("\n--- referencia (lo ya conocido n.s.) ---")
corr(zth, drt, "zMI θ (PAC) × ΔRT")
corr(zth, dies, "zMI θ (PAC) × ΔIES")
print("\n--- absolutos en chew (theta-late chew × IES improvement) ---")
corr(tch[roi,:].mean(0), dies, "θ-late(chew) × ΔIES")
corr(dexp, dtheta_late, "Δexponent × Δθ-late (¿broadband?)")
