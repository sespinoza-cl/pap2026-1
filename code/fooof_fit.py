"""
fooof_fit.py — FOOOF/specparam batch fitting para Paper 2 V3
Llamado desde S6_FOOOF_V3.m via system().

Uso:
    python fooof_fit.py --input <psd_mat> --output <results_mat>

Input .mat (v7.3, HDF5):
    psds      [nFreq × nSubj]   potencia lineal (µV²/Hz), double
    freqs     [nFreq × 1]       vector de frecuencias (Hz)
    fit_range [1 × 2]           [f_low, f_high] Hz para el ajuste

Output .mat:
    exponents [nSubj × 1]       exponente χ (pendiente aperiódica, positivo)
    offsets   [nSubj × 1]       offset b
    ap_fits   [nFreq × nSubj]   curva aperiódica (log10 µV²/Hz)
    residuals [nFreq × nSubj]   residual periódico (log10)
    r_squared [nSubj × 1]       R² del ajuste
    psd_log   [nFreq × nSubj]   log10(psd) original
    freqs_out [nFreq × 1]       vector de frecuencias (idéntico al de entrada)
"""

import sys
import argparse
import numpy as np
import scipy.io as sio

def load_mat(path):
    """Carga .mat v7.3 (HDF5) o v5 automáticamente."""
    try:
        import h5py
        with h5py.File(path, 'r') as f:
            data = {}
            for k in f.keys():
                v = f[k][()]
                if v.dtype.kind in ('O',):       # object/cell arrays
                    data[k] = np.array(v).flatten()
                else:
                    data[k] = np.array(v, dtype=float)
        return data
    except Exception:
        return sio.loadmat(path, squeeze_me=True)

def fit_fooof_subject(freqs, psd_linear, fit_range, max_n_peaks=8,
                       min_peak_height=0.05, peak_threshold=2.0):
    """Ajusta FOOOF (specparam) para un sujeto."""
    # Intentar con specparam (nuevo) o fooof (viejo)
    try:
        from specparam import SpectralModel as FOOOF_cls
    except ImportError:
        try:
            from fooof import FOOOF as FOOOF_cls
        except ImportError:
            raise ImportError(
                "Ningún paquete FOOOF encontrado.\n"
                "Instalar: pip install specparam  o  pip install fooof"
            )

    fm = FOOOF_cls(
        peak_width_limits=[0.5, 12.0],
        max_n_peaks=max_n_peaks,
        min_peak_height=min_peak_height,
        peak_threshold=peak_threshold,
        aperiodic_mode='fixed',
        verbose=False
    )
    fm.fit(freqs, psd_linear, fit_range)

    exp     = fm.aperiodic_params_[1]
    offset  = fm.aperiodic_params_[0]
    ap_log  = fm._ap_fit          # log10 aperiodic component (fit range only)
    res     = fm._spectrum_flat   # flattened spectrum (residual) over fit range

    # Construir curvas sobre el vector completo de frecuencias
    ap_full  = np.full(len(freqs), np.nan)
    res_full = np.full(len(freqs), np.nan)
    psd_log  = np.where(psd_linear > 0, np.log10(psd_linear), np.nan)

    f_mask = (freqs >= fit_range[0]) & (freqs <= fit_range[1])
    ap_full[f_mask]  = ap_log
    res_full[f_mask] = res

    # R² sobre el rango de ajuste
    y_true = psd_log[f_mask]
    y_pred = ap_log + res
    ss_res = np.nansum((y_true - y_pred)**2)
    ss_tot = np.nansum((y_true - np.nanmean(y_true))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return exp, offset, ap_full, res_full, psd_log, r2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input',  required=True,  help='Ruta al .mat de entrada')
    parser.add_argument('--output', required=True,  help='Ruta al .mat de salida')
    parser.add_argument('--fit_low',  type=float, default=3.0,  help='Límite inferior Hz')
    parser.add_argument('--fit_high', type=float, default=35.0, help='Límite superior Hz')
    args = parser.parse_args()

    fit_range = [args.fit_low, args.fit_high]
    print(f"  Rango ajuste: {fit_range[0]}–{fit_range[1]} Hz")

    data  = load_mat(args.input)
    psds  = np.array(data['psds'])        # [nFreq × nSubj]
    freqs = np.array(data['freqs']).flatten()

    # MATLAB guarda en column-major; si viene transpuesto, corregir
    if psds.ndim == 1:
        psds = psds[:, np.newaxis]
    if psds.shape[0] != len(freqs):
        psds = psds.T

    nFreq, nSubj = psds.shape
    print(f"  PSDs: {nFreq} frecuencias × {nSubj} sujetos")

    exponents = np.full(nSubj, np.nan)
    offsets   = np.full(nSubj, np.nan)
    r_squared = np.full(nSubj, np.nan)
    ap_fits   = np.full((nFreq, nSubj), np.nan)
    residuals = np.full((nFreq, nSubj), np.nan)
    psd_log   = np.full((nFreq, nSubj), np.nan)

    for si in range(nSubj):
        psd_s = psds[:, si]
        if np.any(np.isnan(psd_s)) or np.all(psd_s <= 0):
            print(f"  [SKIP] sujeto {si+1}: PSD inválida")
            continue
        try:
            exp, off, ap, res, pl, r2 = fit_fooof_subject(
                freqs, psd_s, fit_range)
            exponents[si] = exp
            offsets[si]   = off
            r_squared[si] = r2
            ap_fits[:, si]   = ap
            residuals[:, si] = res
            psd_log[:, si]   = pl
            print(f"  Sujeto {si+1:3d}: exp={exp:.3f}  b={off:.3f}  R2={r2:.4f}")
        except Exception as e:
            print(f"  [ERROR] sujeto {si+1}: {e}")

    sio.savemat(args.output, {
        'exponents': exponents,
        'offsets':   offsets,
        'r_squared': r_squared,
        'ap_fits':   ap_fits,
        'residuals': residuals,
        'psd_log':   psd_log,
        'freqs_out': freqs,
    })
    print(f"  Guardado: {args.output}")
    print(f"  exp medio: {np.nanmean(exponents):.3f} +/- {np.nanstd(exponents):.3f}")


if __name__ == '__main__':
    main()
