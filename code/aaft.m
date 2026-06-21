function s = aaft(x, stream)
% aaft  Amplitude Adjusted Fourier Transform surrogate (Theiler et al. 1992).
%   s = aaft(x, stream)
%   x      : vector real (columna o fila)
%   stream : RandStream para reproducibilidad (opcional; default global)
%   s      : surrogado con el MISMO espectro de potencia y la MISMA distribución
%            de amplitud que x, pero con fase aleatorizada → rompe la forma de
%            onda no-sinusoidal / no-linealidad. Control para artefacto PAC (C2/A4).
%
% Uso en PAC: se aplica a la señal EEG filtrada por banda; su envolvente ya no
% se acopla a la fase EMG por estructura de forma, sólo por azar espectral.

    if nargin < 2 || isempty(stream), stream = RandStream.getGlobalStream; end
    x  = x(:);
    n  = numel(x);

    % 1) Gaussianizar preservando el rango ordinal de x
    [~, ix] = sort(x);
    rg      = sort(randn(stream, n, 1));
    g       = zeros(n, 1);
    g(ix)   = rg;

    % 2) Aleatorizar fase de la versión gaussianizada (preserva |FFT|)
    gp = phase_randomize(g, stream);

    % 3) Re-imponer la distribución de amplitud original según el rango de gp
    [~, igp] = sort(gp);
    xs       = sort(x);
    s        = zeros(n, 1);
    s(igp)   = xs;
end

function y = phase_randomize(x, stream)
    x    = x(:);
    n    = numel(x);
    X    = fft(x);
    half = floor(n/2);
    ph   = rand(stream, half-1, 1) * 2*pi - pi;     % fases para bins 2..half
    Y        = X;
    Y(2:half) = abs(X(2:half)) .* exp(1i*ph);
    Y(n:-1:n-half+2) = conj(Y(2:half));             % simetría hermitiana
    if mod(n,2) == 0
        Y(half+1) = abs(X(half+1));                 % bin de Nyquist real
    end
    y = real(ifft(Y));
end
