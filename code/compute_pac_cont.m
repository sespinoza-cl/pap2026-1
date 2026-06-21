function R = compute_pac_cont(phase, eeg_roi, bands, fs, n_bins, n_surr, min_sh, stream)
% compute_pac_cont  PAC continuo con DOBLE NULL (R7) para UNA condición.
%   phase   : [1×T] fase EMG (Hilbert de f_chew±0.5 Hz)
%   eeg_roi : [1×T] EEG ROI-18 promediado (crudo; se filtra por banda aquí)
%   bands   : cell de [lo hi] Hz (theta/alpha/beta)
%   fs      : Hz
%   n_bins  : bins de fase Tort (18)
%   n_surr  : nº surrogados por null (500)
%   min_sh  : shift mínimo en MUESTRAS (5 s · fs)
%   stream  : RandStream (reproducibilidad por sujeto)
%
%   R.mi  [nB]  MI observado
%   R.zc  [nB]  z vs null CIRCULAR-SHIFT  (rompe timing, preserva todo)
%   R.za  [nB]  z vs null AAFT            (preserva espectro, rompe forma de onda)
%   R.ncm/R.ncs, R.nam/R.nas : media/sd de cada null
%   R.pref[nB]  fase preferida (ángulo del vector resultante de amplitud por bin)

    nB = numel(bands);
    T  = numel(eeg_roi);
    be = linspace(-pi, pi, n_bins+1);
    bc = (be(1:end-1) + be(2:end)) / 2;
    [R.mi, R.zc, R.za, R.ncm, R.ncs, R.nam, R.nas, R.pref] = deal(nan(nB,1));

    ph_col = phase(:);
    range_sh = max(1, T - 2*min_sh);

    for b = 1:nB
        [bb, ab] = butter(4, bands{b}/(fs/2), 'bandpass');
        sigf = filtfilt(bb, ab, eeg_roi);
        amp  = abs(hilbert(sigf));
        amp_col = amp(:);
        R.mi(b) = tort_mi(ph_col, amp_col, n_bins);

        % Fase preferida
        ad = zeros(n_bins,1);
        for k = 1:n_bins
            ik = phase >= be(k) & phase < be(k+1);
            if any(ik), ad(k) = mean(amp(ik)); end
        end
        ad = ad/(sum(ad)+eps);
        R.pref(b) = angle(sum(ad .* exp(1i*bc')));

        % NULL A — circular-shift de la amplitud (min_sh muestras)
        mc = nan(n_surr,1);
        for k = 1:n_surr
            sh = min_sh + randi(stream, range_sh) - 1;
            mc(k) = tort_mi(ph_col, circshift(amp_col, sh), n_bins);
        end
        R.ncm(b) = mean(mc,'omitnan');
        R.ncs(b) = std(mc,'omitnan');
        R.zc(b)  = (R.mi(b) - R.ncm(b)) / (R.ncs(b)+eps);

        % NULL B — AAFT de la señal filtrada (preserva |FFT| y distribución)
        ma = nan(n_surr,1);
        for k = 1:n_surr
            sigs = aaft(sigf, stream);
            ma(k) = tort_mi(ph_col, abs(hilbert(sigs))', n_bins);
        end
        R.nam(b) = mean(ma,'omitnan');
        R.nas(b) = std(ma,'omitnan');
        R.za(b)  = (R.mi(b) - R.nam(b)) / (R.nas(b)+eps);
    end
end
