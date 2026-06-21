"""
S5_PLI_theta.py  —  Phase Lag Index (PLI) theta en WIN_LATE
================================================================
Calcula PLI entre todos los pares de canales EEG en la banda theta
(4-7 Hz) durante WIN_LATE (900-1300 ms), comparando Ch vs Nc.

Basado en Paper 1 (Espinoza et al., Sci Rep 2025) que encontró
PLI theta frontoparietal aumentado con masticación.

Inputs:
  Data_PAC/Epochs/E3S*_Ch_ep.set  (EEGLAB, 66 ch, 256 Hz)
  Data_PAC/Epochs/E3S*_Nc_ep.set

Outputs:
  Analysis_V1_Final/outputs/stats/v1_S5_PLI.mat
    PLI_ch_cas    (31, 64, 64) — PLI theta WIN_LATE, bloque Ch, Cases
    PLI_nc_cas    (31, 64, 64) — PLI theta WIN_LATE, bloque Nc, Cases
    PLI_ch_ctr    (15, 64, 64) — Controls Ch
    PLI_nc_ctr    (15, 64, 64) — Controls Nc
    delta_PLI_cas (31, 64, 64) — Ch - Nc, Cases
    delta_PLI_ctr (15, 64, 64) — Ch - Nc, Controls
    t_cas         (64, 64)     — t-stat paired t-test Ch vs Nc, Cases
    p_cas         (64, 64)     — p-valor (raw)
    p_cas_fdr     (64, 64)     — p-valor FDR-corrected
    ch_names      (64,)        — nombres de canales EEG

Método PLI:
  PLI(i,j) = |mean_t,ep(sign(imag(analytic_i(t) * conj(analytic_j(t)))))|
  = |mean_t,ep(sign(sin(phi_i(t) - phi_j(t))))|
  Robusto a volumen de conducción (descarta phase-lag=0)

Wael/Stam2007 reference: PLI ∈ [0,1], 0 = no acoplamiento, 1 = acoplamiento perfecto
"""

import numpy as np
import scipy.io as sio
from scipy.signal import butter, filtfilt, hilbert
from scipy import stats
import mne
import os, sys

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT       = r"C:\Users\Pc - Casa\Desktop\Proyectos_Claude\Phd\Paper2\P2V1"
DIR_EPOCHS = os.path.join(ROOT, "Data_PAC", "Epochs")
OUT_STATS  = os.path.join(ROOT, "Analysis_V1_Final", "outputs", "stats")
os.makedirs(OUT_STATS, exist_ok=True)

# ── Sujetos ────────────────────────────────────────────────────────────────
CASES    = ['E3S1','E3S2','E3S4','E3S6','E3S7','E3S8','E3S9','E3S10',
            'E3S11','E3S12','E3S13','E3S14','E3S15','E3S16','E3S17','E3S18',
            'E3S19','E3S20','E3S21','E3S22','E3S23','E3S24','E3S25','E3S26',
            'E3S27','E3S28','E3S29','E3S30','E3S31','E3S32','E3S33']
CONTROLS = ['E3C1','E3C2','E3C3','E3C4','E3C5','E3C6','E3C7','E3C8',
            'E3C9','E3C10','E3C11','E3C12','E3C13','E3C14','E3C15']

# ── Parámetros ─────────────────────────────────────────────────────────────
FBAND    = (4.0, 7.0)     # theta Hz
WIN_LATE = (0.900, 1.300) # segundos
N_EEG    = 64             # primeros 64 canales son EEG

# ── Funciones ──────────────────────────────────────────────────────────────

def bandpass_theta(data, sfreq, fband=FBAND):
    """Butterworth bandpass 4th order, zero-phase."""
    nyq = sfreq / 2
    b, a = butter(4, [fband[0]/nyq, fband[1]/nyq], btype='band')
    return filtfilt(b, a, data, axis=-1)

def compute_pli(epochs_data, times, win=(0.900, 1.300), sfreq=256.0, fband=FBAND):
    """
    Calcula PLI para todos los pares de canales EEG.

    epochs_data : (n_epochs, n_ch, n_times)  float, ya en microvolt
    times       : (n_times,)  en segundos
    Devuelve    : (n_ch, n_ch) PLI simétrico, diagonal = 0
    """
    n_ep, n_ch, n_t = epochs_data.shape

    # Filtrar theta
    filt = bandpass_theta(epochs_data, sfreq, fband)   # (n_ep, n_ch, n_t)

    # Ventana temporal
    t_idx = (times >= win[0]) & (times <= win[1])
    filt_late = filt[:, :, t_idx]                       # (n_ep, n_ch, n_late)

    # Señal analítica via Hilbert (sobre el eje de tiempo)
    analytic = hilbert(filt_late, axis=2)               # (n_ep, n_ch, n_late)

    # Reshape → (n_ch, n_ep × n_late) para vectorización
    n_late = np.sum(t_idx)
    A = analytic.transpose(1, 0, 2).reshape(n_ch, -1)   # (n_ch, N)
    N = A.shape[1]

    # PLI vectorizado: imag(A[i] * conj(A[j])) per sample, luego sign, luego mean
    # Para evitar OOM: procesamos por bloques de i
    PLI = np.zeros((n_ch, n_ch), dtype=np.float32)

    for i in range(n_ch):
        # cross: (n_ch, N) = A * conj(A[i])
        cross_imag = np.imag(A * np.conj(A[i:i+1, :]))  # (n_ch, N)
        PLI[i, :] = np.abs(np.mean(np.sign(cross_imag), axis=1))

    np.fill_diagonal(PLI, 0)
    # Simetrizar (por construcción ya es simétrico, pero asegurar)
    PLI = (PLI + PLI.T) / 2
    return PLI

def load_epochs(subj, cond, dir_epochs):
    """Carga épocas EEGLAB y devuelve (data, times, ch_names, sfreq)."""
    fname = os.path.join(dir_epochs, f"{subj}_{cond}_ep.set")
    ep = mne.io.read_epochs_eeglab(fname, verbose=False)
    data   = ep.get_data()[:, :N_EEG, :]   # (n_ep, 64, n_t)
    times  = ep.times                        # segundos
    ch_names = ep.ch_names[:N_EEG]
    sfreq  = ep.info['sfreq']
    return data, times, ch_names, sfreq

# ── Loop principal ──────────────────────────────────────────────────────────

def run_group(subjects, label):
    n = len(subjects)
    PLI_ch_all = []
    PLI_nc_all = []
    failures = []

    for k, subj in enumerate(subjects):
        sys.stdout.write(f"\r  {label} {k+1}/{n}  {subj}    ")
        sys.stdout.flush()
        try:
            data_ch, times, ch_names, sfreq = load_epochs(subj, 'Ch', DIR_EPOCHS)
            data_nc, _,     _,        _     = load_epochs(subj, 'Nc', DIR_EPOCHS)
            PLI_ch = compute_pli(data_ch, times, sfreq=sfreq)
            PLI_nc = compute_pli(data_nc, times, sfreq=sfreq)
            PLI_ch_all.append(PLI_ch)
            PLI_nc_all.append(PLI_nc)
        except Exception as e:
            print(f"\n    ERROR {subj}: {e}")
            failures.append(subj)
            PLI_ch_all.append(np.full((N_EEG, N_EEG), np.nan, dtype=np.float32))
            PLI_nc_all.append(np.full((N_EEG, N_EEG), np.nan, dtype=np.float32))

    print(f"\n  Done {label}. Failures: {failures if failures else 'none'}")
    return np.stack(PLI_ch_all), np.stack(PLI_nc_all), ch_names

print("=== S5_PLI_theta — PLI theta WIN_LATE ===")
print("Banda: %.0f-%.0f Hz  |  WIN_LATE: %.0f-%.0f ms" % (
      FBAND[0], FBAND[1], WIN_LATE[0]*1000, WIN_LATE[1]*1000))
print()

print("[1/2] Procesando Cases (n=31)...")
PLI_ch_cas, PLI_nc_cas, ch_names = run_group(CASES, 'Cases')

print("[2/2] Procesando Controls (n=15)...")
PLI_ch_ctr, PLI_nc_ctr, _ = run_group(CONTROLS, 'Controls')

# ── Estadísticas ────────────────────────────────────────────────────────────
print("\nComputando estadísticas...")

delta_cas = PLI_ch_cas - PLI_nc_cas   # (31, 64, 64)
delta_ctr = PLI_ch_ctr - PLI_nc_ctr  # (15, 64, 64)

# Paired t-test per connection (Cases: Ch vs Nc)
t_cas = np.zeros((N_EEG, N_EEG))
p_cas = np.ones((N_EEG, N_EEG))

for i in range(N_EEG):
    for j in range(i+1, N_EEG):
        t_val, p_val = stats.ttest_rel(PLI_ch_cas[:, i, j],
                                        PLI_nc_cas[:, i, j],
                                        nan_policy='omit')
        t_cas[i, j] = t_val
        t_cas[j, i] = t_val
        p_cas[i, j] = p_val
        p_cas[j, i] = p_val

# FDR correction (Benjamini-Hochberg) sobre triángulo superior
from scipy.stats import false_discovery_control
triu_idx = np.triu_indices(N_EEG, k=1)
p_upper  = p_cas[triu_idx]
p_fdr    = false_discovery_control(p_upper, method='bh')
p_cas_fdr = np.ones((N_EEG, N_EEG))
p_cas_fdr[triu_idx] = p_fdr
p_cas_fdr[triu_idx[1], triu_idx[0]] = p_fdr

# Resumen
sig_raw = np.sum(p_cas[triu_idx] < 0.05)
sig_fdr = np.sum(p_fdr < 0.05)
print(f"  Conexiones Ch>Nc significativas (p<0.05 raw): {sig_raw}/{len(p_upper)}")
print(f"  Conexiones Ch>Nc significativas (FDR q<0.05): {sig_fdr}/{len(p_upper)}")

# Wilcoxon sobre mean PLI (ROI = promedio de todos los pares)
mean_PLI_ch = np.nanmean(PLI_ch_cas[:, triu_idx[0], triu_idx[1]], axis=1)  # (31,)
mean_PLI_nc = np.nanmean(PLI_nc_cas[:, triu_idx[0], triu_idx[1]], axis=1)
w_stat, p_wilcoxon = stats.wilcoxon(mean_PLI_ch, mean_PLI_nc)
print(f"  Mean global PLI Ch={mean_PLI_ch.mean():.4f}  Nc={mean_PLI_nc.mean():.4f}")
print(f"  Wilcoxon mean PLI (Ch vs Nc, Cases): p={p_wilcoxon:.4f}")

# Same for Controls
mean_PLI_ch_ctr = np.nanmean(PLI_ch_ctr[:, triu_idx[0], triu_idx[1]], axis=1)
mean_PLI_nc_ctr = np.nanmean(PLI_nc_ctr[:, triu_idx[0], triu_idx[1]], axis=1)
if len(mean_PLI_ch_ctr) > 2:
    _, p_wilcoxon_ctr = stats.wilcoxon(mean_PLI_ch_ctr, mean_PLI_nc_ctr)
    print(f"  Controls: Ch={mean_PLI_ch_ctr.mean():.4f}  Nc={mean_PLI_nc_ctr.mean():.4f}  p={p_wilcoxon_ctr:.4f}")

# ── Guardar ─────────────────────────────────────────────────────────────────
out_file = os.path.join(OUT_STATS, 'v1_S5_PLI.mat')
sio.savemat(out_file, {
    'PLI_ch_cas':  PLI_ch_cas.astype(np.float32),
    'PLI_nc_cas':  PLI_nc_cas.astype(np.float32),
    'PLI_ch_ctr':  PLI_ch_ctr.astype(np.float32),
    'PLI_nc_ctr':  PLI_nc_ctr.astype(np.float32),
    'delta_PLI_cas': delta_cas.astype(np.float32),
    'delta_PLI_ctr': delta_ctr.astype(np.float32),
    't_cas':       t_cas.astype(np.float32),
    'p_cas':       p_cas.astype(np.float32),
    'p_cas_fdr':   p_cas_fdr.astype(np.float32),
    'ch_names':    np.array(ch_names),
    'mean_PLI_ch_cas': mean_PLI_ch,
    'mean_PLI_nc_cas': mean_PLI_nc,
    'mean_PLI_ch_ctr': mean_PLI_ch_ctr,
    'mean_PLI_nc_ctr': mean_PLI_nc_ctr,
    'p_wilcoxon_cas': p_wilcoxon,
    'FBAND':       np.array(FBAND),
    'WIN_LATE_S':  np.array(WIN_LATE),
})
print(f"\nGuardado: {out_file}")
print("=== DONE ===")
