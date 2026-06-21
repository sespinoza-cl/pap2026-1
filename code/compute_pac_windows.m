function MIw = compute_pac_windows(eeg_roi, phase, events, wS, bands, fs, n_bins, wins_ms)
% compute_pac_windows  PAC por ventana, INDUCIDO sin el bug de épocas solapadas (F5).
%   Epoca PRIMERO, resta la media-por-época (ERP) a cada época, LUEGO filtra
%   cada época y agrupa (fase,amp) por sub-ventana → MI Tort.
%
%   eeg_roi : [1×T] EEG ROI-18 crudo continuo
%   phase   : [1×T] fase EMG continua
%   events  : [1×nev] latencias (muestras) de onset de trial válidas
%   wS      : [start end] offsets en muestras relativos al onset (de WIN_EPOCH)
%   bands   : cell {[lo hi]...}
%   fs, n_bins
%   wins_ms : 3×2 [base; early; late] en ms relativos al onset
%   MIw     : [nB × 3] MI por banda × {base, early, late}

    nB   = numel(bands);
    nev  = numel(events);
    Lwin = diff(wS) + 1;
    MIw  = nan(nB, 3);
    if nev < 5, return; end

    % Epocar señal cruda + fase
    ep_raw = nan(Lwin, nev);
    ep_ph  = nan(Lwin, nev);
    for k = 1:nev
        ix = events(k)+wS(1) : events(k)+wS(2);
        ep_raw(:,k) = eeg_roi(ix)';
        ep_ph(:,k)  = phase(ix)';
    end

    % Inducido: restar ERP (media a través de épocas) a CADA época
    erp    = mean(ep_raw, 2, 'omitnan');
    ep_ind = ep_raw - erp;

    t_ms = ((0:Lwin-1)/fs*1000) + (wS(1)/fs*1000);

    for b = 1:nB
        [bb,ab] = butter(4, bands{b}/(fs/2), 'bandpass');
        amp_ep  = nan(Lwin, nev);
        for k = 1:nev
            amp_ep(:,k) = abs(hilbert(filtfilt(bb,ab, ep_ind(:,k))));
        end
        for w = 1:3
            iw = t_ms >= wins_ms(w,1) & t_ms <= wins_ms(w,2);
            if ~any(iw), continue; end
            ph_p = reshape(ep_ph(iw,:),  [], 1);
            am_p = reshape(amp_ep(iw,:), [], 1);
            ok   = ~isnan(ph_p) & ~isnan(am_p);
            MIw(b,w) = tort_mi(ph_p(ok), am_p(ok), n_bins);
        end
    end
end
