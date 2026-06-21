function [emg, used, snr_used] = emg_bilateral(data, chs, fs)
% emg_bilateral  Criterio ÚNICO de canal EMG (R6): promedio bilateral 65+66,
%   con fallback al canal vivo si el otro está muerto. Idéntico en Ch y Nc.
%   data : matriz [nchan × T]
%   chs  : canales EMG candidatos (p.ej. [65 66]) ya validados <= nbchan
%   fs   : Hz
%   emg  : [1×T] señal EMG combinada (detrended)
%   used : 0=ambos · 65/66=solo ese canal · -1=fallback ciego (ambos sospechosos)
%   snr_used : SNR masticatorio (dB, in[1-2.5]/out) del/los canal(es) usados
%
% "Muerto" = std≈0 o SNR no finito. SNR masticatorio: potencia de la envolvente
% rectificada en [1,2.5] Hz vs [0.1,1]∪[2.5,5] Hz.

    [bb,ab] = butter(4, [20 min(200, fs/2-1)]/(fs/2), 'bandpass');
    good = []; snrs = []; sigs = {};
    for c = chs(:)'
        sig = detrend(double(data(c,:)));
        if std(sig) < 1e-6, continue; end          % canal plano = muerto
        env = detrend(abs(filtfilt(bb,ab,sig)));
        [px,fx] = pwelch(env, round(fs*2), round(fs), [], fs);
        inb  = fx>=1 & fx<=2.5;
        outb = (fx>=0.1 & fx<1) | (fx>2.5 & fx<=5);
        snr = NaN;
        if any(inb) && any(outb)
            snr = 10*log10(mean(px(inb))/mean(px(outb)));
        end
        if isfinite(snr)
            good(end+1) = c;  snrs(end+1) = snr;  sigs{end+1} = sig; %#ok<AGROW>
        end
    end

    if isempty(good)
        % fallback ciego: promedio crudo de los candidatos
        emg = detrend(mean(double(data(chs,:)),1));
        used = -1; snr_used = NaN;
    elseif numel(good) == 2
        emg = (sigs{1} + sigs{2}) / 2;
        used = 0; snr_used = mean(snrs);
    else
        emg = sigs{1};
        used = good(1); snr_used = snrs(1);
    end
end
