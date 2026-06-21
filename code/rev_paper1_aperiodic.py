"""
rev_paper1_aperiodic.py — Test (e) del relato online/offline:
¿el exponente aperiódico 1/f se APLANA post-masticación (chew vs no-chew) en el Paper 1 (OFFLINE),
igual que online (Paper 2, χ 1.110→0.735)? Si sí → semilla de excitabilidad central común a ambos.
Datos: D:/Exp1/Exp1/Revision_Paper/epochs (S1..S30; tb[a/n][ct/nct]). FOOOF 1.1, ROI frontocentral.
"""
import os, glob, warnings, numpy as np
warnings.filterwarnings("ignore")
import mne
from fooof import FOOOF
from scipy.stats import wilcoxon
mne.set_log_level("ERROR")

EPDIR = r"D:\Exp1\Exp1\Revision_Paper\epochs"
ROI = ["AFz", "F1", "Fz", "F2", "FC1", "FCz", "FC2"]
CHEW = ["tbnct", "tbact"]      # chew (no-anes + anes)
NOCHEW = ["tbnnct", "tbanct"]  # no-chew (no-anes + anes)
FIT = [3, 40]

def load_psd(subj, conds):
    psds = []
    for c in conds:
        f = os.path.join(EPDIR, f"{subj}_{c}.set")
        if not os.path.exists(f):
            continue
        ep = mne.io.read_epochs_eeglab(f)
        chn = [c for c in ROI if c in ep.ch_names]
        if len(chn) < 3:
            # fallback: frontocentral por prefijo
            chn = [c for c in ep.ch_names if c.upper().startswith(("FZ","FCZ","AFZ","F1","F2","FC1","FC2"))][:7]
        ep.pick(chn)
        sp = ep.compute_psd(method="welch", fmin=1, fmax=45, verbose=False)
        p = sp.get_data().mean(axis=(0, 1))   # mean over epochs+channels
        psds.append((p, sp.freqs))
    if not psds:
        return None, None
    # promedio ponderado simple de las condiciones (chew/anes pooled)
    f0 = psds[0][1]
    P = np.mean([p for p, _ in psds], axis=0)
    return P, f0

def expo(P, freqs):
    fm = FOOOF(peak_width_limits=[0.5, 12], max_n_peaks=8, min_peak_height=0.01,
               peak_threshold=1.5, aperiodic_mode="fixed", verbose=False)
    fm.fit(freqs, P, FIT)
    return fm.aperiodic_params_[-1], fm.r_squared_

subs = sorted({os.path.basename(f).split("_")[0] for f in glob.glob(os.path.join(EPDIR, "S*_tb*.set"))},
              key=lambda s: int(s[1:]))
print(f"Sujetos: {len(subs)}  ROI={ROI}")
ech, enc, r2c, r2n, used = [], [], [], [], []
for s in subs:
    Pc, fc = load_psd(s, CHEW); Pn, fn = load_psd(s, NOCHEW)
    if Pc is None or Pn is None:
        print(f"  {s}: faltan sets"); continue
    xc, rc = expo(Pc, fc); xn, rn = expo(Pn, fn)
    ech.append(xc); enc.append(xn); r2c.append(rc); r2n.append(rn); used.append(s)
    print(f"  {s}: χ_chew={xc:.3f} (R²={rc:.2f})  χ_nochew={xn:.3f} (R²={rn:.2f})  Δ={xc-xn:+.3f}")

ech, enc = np.array(ech), np.array(enc)
d = ech - enc
nflat = int(np.sum(d < 0))
print(f"\n=== APERIÓDICO OFFLINE (Paper 1, n={len(ech)}) ===")
print(f"χ chew  : media={ech.mean():.3f}  mediana={np.median(ech):.3f}")
print(f"χ nochew: media={enc.mean():.3f}  mediana={np.median(enc):.3f}")
print(f"Δ(chew-nochew): media={d.mean():+.3f}  (negativo = aplanamiento)  aplanan {nflat}/{len(d)}")
W, p2 = wilcoxon(ech, enc)                       # two-sided
try: Wg, p1 = wilcoxon(ech, enc, alternative="less")  # chew < nochew (flatter)
except Exception: p1 = np.nan
print(f"Wilcoxon chew vs nochew: two-sided p={p2:.4f} | one-sided(chew<nochew) p={p1:.4f}")
print(f"R² medio ajuste FOOOF: chew={np.mean(r2c):.2f} nochew={np.mean(r2n):.2f}")
print(f"\nCOMPARA con ONLINE (Paper 2): χ 1.110→0.735 (Δ=-0.375, p<0.0001).")
print("Si offline también aplana (Δ<0 sig) -> semilla de excitabilidad central COMÚN (cascada apoyada).")
np.savez(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
         "outputs", "stats", "rev_paper1_aperiodic.npz"),
         exp_chew=ech, exp_nochew=enc, subs=used)
