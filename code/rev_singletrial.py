"""
rev_singletrial.py — Puente conductual↔EEG a nivel de ENSAYO (chew block, casos).
Para cada ensayo HIT (target chew=40 con respuesta=45): RT = t(45)-t(estímulo) y theta-power
frontal-medial (ROI-18, 4-7 Hz) en ventanas pre/early/late. LMM RT ~ theta_within + (1|sujeto).
Mucho más potente que las Δ-correlaciones entre-sujetos.
"""
import warnings, os, glob, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import mne; mne.set_log_level("ERROR")
import statsmodels.formula.api as smf
from scipy.signal import hilbert

ROOT = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
DATA = os.path.join(ROOT, "Data_PAC")
ROI18 = ["Fp1","AF7","AF3","F1","F3","F5","F7","FC1","Fpz","Fp2","AF8","AF4","AFz","Fz","F2","F4","FC2","FCz"]
WINS = {"pre": (-0.5,-0.1), "early": (0.1,0.5), "late": (0.9,1.3)}
TH = (4,7)

# casos = E3S* con _Ch_clean_emg.set
files = sorted(glob.glob(os.path.join(DATA, "E3S*_Ch_clean_emg.set")))
rows = []
for f in files:
    sid = os.path.basename(f).split("_")[0]
    raw = mne.io.read_raw_eeglab(f, preload=True)
    fs = raw.info["sfreq"]
    picks = [c for c in ROI18 if c in raw.ch_names]
    if len(picks) < 8: print(f"  {sid}: ROI<8, skip"); continue
    raw.pick(picks)
    raw.filter(TH[0], TH[1], verbose=False)
    amp = np.abs(hilbert(raw.get_data(), axis=1))      # ch x time, theta envelope
    power = (amp**2).mean(0)                            # ROI-mean theta power (time,)
    ev, evid = mne.events_from_annotations(raw, verbose=False); inv = {v:k for k,v in evid.items()}
    codes = np.array([inv[c] for c in ev[:,2]]); lat = ev[:,0]/fs
    tgt = np.where(codes == "40")[0]                   # targets in chew block
    resp_t = lat[codes == "45"]
    n_hit = 0
    for i in tgt:
        ts = lat[i]
        # respuesta 45 en (ts, ts+1.2]
        r = resp_t[(resp_t > ts) & (resp_t <= ts+1.2)]
        if len(r) == 0: continue                       # miss -> sin RT
        rt = (r[0]-ts)*1000
        if rt < 150 or rt > 1200: continue             # RT plausible
        feats = {}
        ok = True
        for w,(a,b) in WINS.items():
            s0,s1 = int((ts+a)*fs), int((ts+b)*fs)
            if s0 < 0 or s1 > len(power): ok=False; break
            feats[w] = np.log10(power[s0:s1].mean()+1e-20)
        if not ok: continue
        rows.append({"subj":sid, "RT":rt, **{f"theta_{w}":feats[w] for w in WINS}})
        n_hit += 1
    print(f"  {sid}: hits con RT={n_hit}")

df = pd.DataFrame(rows)
print(f"\nTotal ensayos: {len(df)} de {df.subj.nunique()} sujetos (media {len(df)/df.subj.nunique():.0f}/suj)")

# centrar theta DENTRO de sujeto -> efecto puramente within-subject
for w in WINS:
    df[f"theta_{w}_c"] = df.groupby("subj")[f"theta_{w}"].transform(lambda x: x - x.mean())

print("\n=== LMM  RT ~ theta_within + (1|sujeto)  (theta centrada por sujeto) ===")
for w in WINS:
    md = smf.mixedlm(f"RT ~ theta_{w}_c", df, groups=df["subj"])
    r = md.fit(method="lbfgs")
    b = r.params[f"theta_{w}_c"]; p = r.pvalues[f"theta_{w}_c"]; se = r.bse[f"theta_{w}_c"]
    sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "n.s."
    print(f"  theta {w:5} ({WINS[w][0]*1000:.0f}-{WINS[w][1]*1000:.0f} ms): "
          f"slope={b:+.1f} ms/log-unit (SE {se:.1f})  p={p:.4f}  {sig}")

# correlación within-subject media (Fisher) como referencia
print("\n=== correlación within-subject (media de r por sujeto, Fisher-z) ===")
from scipy.stats import pearsonr
for w in WINS:
    rs=[]
    for s,g in df.groupby("subj"):
        if len(g)>5:
            rr,_=pearsonr(g[f"theta_{w}"], g["RT"]); rs.append(np.arctanh(np.clip(rr,-.999,.999)))
    rbar=np.tanh(np.mean(rs))
    print(f"  theta {w:5}: r_within medio={rbar:+.3f} (n={len(rs)} suj)")
df.to_csv(os.path.join(ROOT,"Analysis_V1_Final","outputs","stats","rev_singletrial.csv"), index=False)
print("\nguardado outputs/stats/rev_singletrial.csv")
