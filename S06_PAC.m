%% ============================================================================
%  S6_PAC_4Groups.m — PAC fase-EMG × amplitud-EEG (4 GRUPOS)
%  ----------------------------------------------------------------------------
%  Extiende S5b_PAC_Continuous.m al diseño 2×2 completo:
%
%    Grupos    : Casos (E3S) | Controles (E3C)
%    Condiciones: Ch (marca 40) | Nc (marca 30)
%
%  LÓGICA DEL CONTROL:
%    Casos/Ch  → EMG rítmico a ~1-2 Hz: fase EMG tiene estructura biológica.
%    Casos/Nc  → EMG quieto: fase extraída es ruido → PAC ≈ 0 (verificación).
%    Controles/Ch  → No masticaron, EMG basal: fase = ruido a ~1-2 Hz.
%    Controles/Nc  → Ídem.
%    Si Rayleigh es sig solo en Casos/Ch → el acoplamiento requiere masticación.
%
%  Para Controles se usa peak_hz = mediana de Casos (no tienen frecuencia
%  individual). El filtro extrae actividad basal del masetero → fase aleatoria.
%
%  IDÉNTICO a S5b_PAC_Continuous en:
%    - Filtrado de señal continua completa (sin edge effects)
%    - Surrogados por permutación de trials (N=200)
%    - Tort MI, n_bins=18
%    - Ventanas: Early[0-300] Late[300-900] Active[200-700] ms
%    - Bandas: Theta[4-7] Alpha[8-13] Beta[14-30] Hz
%    - ROI EEG: F1 + FC1 (fallback ch 12+16)
%
%  OUTPUTS:
%    EEG/PAC/FourGroups/PAC_4Groups_Workspace.mat
%    EEG/PAC/FourGroups/Reports/Reporte_PAC_4Groups.txt
%    EEG/PAC/FourGroups/Plots/   (figuras por grupo)
%
%  Sebastián, 2026-05
%  ============================================================================
clear; clc; close all; rng(2026); tic;

%% ── 1. CONFIGURACIÓN ────────────────────────────────────────────────────────
S = S0_paths();                 % rutas + parámetros UNIFICADOS (fuente única)
P = struct();

P.dir_cont    = S.dir_cont;     % *_final.set (referenciado en C:\Desktop\Exp2, solo lectura)
P.dir_lists   = S.dir_data;     % incluidos45.mat (copiado al proyecto)
P.dir_out     = S.fig04;        % workspace + reportes + figuras → proyecto/outputs
P.dir_reports = S.fig04;
P.dir_plots   = S.fig04;
P.eeglab_path = S.eeglab_path;
P.f_chewfreq  = S.file_chew;
P.f_behavior  = S.file_beh;

% ROI EEG: buscar por label, fallback por índice (desde config)
P.roi_labels  = S.roi_labels;
P.roi_idx_fb  = S.roi_ch;

% EMG: labels posibles en S5_Final (buscar en orden)
P.emg_labels  = {'EMG_L','EMG_R','Mas_L','Mas_R','EXG7','EXG8'};
P.emg_idx_fb  = [67, 68];   % fallback: índices en S5_Final (68 ch total)

% Markers
P.trig_ch     = 40;   % target bloque chewing
P.trig_nc     = 30;   % target bloque no-chewing

% Ventanas temporales (ms relativas al marker)
P.wins.labels = {'Early','Late','Active'};
P.wins.ms     = [0 300; 300 900; 200 700];

% Bandas amplitud EEG (definición ÚNICA desde S0_paths)
P.amp_freqs   = 4:2:60;
P.bands       = S.bands;
P.bands_hz    = S.bands_hz;     % θ[4-7] α[8-12] β[13-30]

% Parámetros PAC
P.bw_pac      = 0.5;    % semi-ancho Hz del filtro fase EMG (peak ± bw)
P.n_bins      = 18;     % bins Tort MI
P.n_surr      = 200;    % surrogados permutación de trials
P.filt_order  = 4;      % orden Butterworth

% Sujetos excluidos (desde config)
P.excluded    = S.excluded;

% Crear directorios
for d = {P.dir_out, P.dir_reports, P.dir_plots}
    if ~exist(d{1},'dir'); mkdir(d{1}); end
end

% EEGLAB
if exist(P.eeglab_path,'dir')
    addpath(P.eeglab_path);
    addpath(genpath(fullfile(P.eeglab_path,'functions')));
    eeglab nogui;
end

sep = repmat('=',1,82);
fprintf('\n%s\n  S6_PAC_4Groups.m — PAC 2×2 (Casos|Controles × Ch|Nc)\n  %s\n%s\n', ...
    sep, datestr(now), sep);
fprintf('  Dir datos : %s\n', P.dir_cont);
fprintf('  Ventanas  : Early[0-300] Late[300-900] Active[200-700] ms\n');
fprintf('  Bandas    : Theta[%d-%d] Alpha[%d-%d] Beta[%d-%d] Hz\n', ...
    P.bands_hz{1}, P.bands_hz{2}, P.bands_hz{3});
fprintf('  Surrogados: permutación trials (N=%d)\n\n', P.n_surr);


%% ── 2. LISTAS DE SUJETOS ────────────────────────────────────────────────────
tmp      = load(fullfile(P.dir_lists, 'incluidos45.mat'));
fn       = fieldnames(tmp); todos = tmp.(fn{1})(:);

casos      = sort(setdiff(todos(startsWith(todos, 'E3S')), P.excluded, 'stable'));
controles  = sort(todos(startsWith(todos, 'E3C')));
nCas = numel(casos); nCtr = numel(controles);

fprintf('  Casos     : %d sujetos  (%s)\n', nCas, strjoin(casos(1:min(3,end)),', '));
fprintf('  Controles : %d sujetos  (%s)\n', nCtr, strjoin(controles(1:min(3,end)),', '));
fprintf('  Excluidos : %s\n\n', strjoin(P.excluded,', '));

GROUPS = {'Casos','Controles'};
SUBJ   = struct('Casos', {casos}, 'Controles', {controles});
NCAS   = struct('Casos', nCas,    'Controles', nCtr);
CONDS  = {'Ch','Nc'};


%% ── 3. FRECUENCIAS MASTICATORIAS ────────────────────────────────────────────
% Casos: individual desde chew_metrics.mat
% Controles: mediana de Casos (no mastican)

if ~exist(P.f_chewfreq,'file')
    error('No encontrado: %s\n→ Verificar ruta chew_metrics.mat', P.f_chewfreq);
end
load(P.f_chewfreq, 'T_freq');

peak_hz_cas = nan(nCas, 1);
for j = 1:nCas
    idx = strcmp(T_freq.Sujeto, casos{j});
    if any(idx)
        peak_hz_cas(j) = mean([T_freq.Freq_Left(idx), T_freq.Freq_Right(idx)], 'omitnan');
    end
end
peak_med = median(peak_hz_cas, 'omitnan');
peak_hz_ctr = repmat(peak_med, nCtr, 1);   % Controles: mediana de Casos

fprintf('  Peak EMG Casos: mediana=%.3f Hz | rango=[%.3f–%.3f Hz]\n', ...
    peak_med, min(peak_hz_cas,[],'omitnan'), max(peak_hz_cas,[],'omitnan'));
fprintf('  Peak EMG Controles (fijo): %.3f Hz (mediana de Casos)\n\n', peak_med);

PEAK = struct('Casos', peak_hz_cas, 'Controles', peak_hz_ctr);


%% ── 4. PRE-ALLOCACIÓN ───────────────────────────────────────────────────────
nW = size(P.wins.ms, 1);
nF = numel(P.amp_freqs);

R = struct();
for g = 1:2
    grp = GROUPS{g};
    nS  = NCAS.(grp);
    for c = 1:2
        cn = CONDS{c};
        R.(grp).(cn).zMI        = nan(nS, nF, nW);
        R.(grp).(cn).MI         = nan(nS, nF, nW);
        R.(grp).(cn).MI_surr_mu = nan(nS, nF, nW);
        R.(grp).(cn).MI_surr_sd = nan(nS, nF, nW);
        R.(grp).(cn).pref_phase = nan(nS, nF, nW);
        R.(grp).(cn).n_trials   = zeros(nS, 1);
        R.(grp).(cn).n_valid    = zeros(nS, nW);
        R.(grp).(cn).peak_used  = nan(nS, 1);
    end
end


%% ── 5. LOOP PRINCIPAL: GRUPOS × SUJETOS (parfor) × CONDICIONES ─────────────
fprintf('%s\n  LOOP PRINCIPAL (parfor sobre sujetos)\n%s\n', repmat('-',1,82), repmat('-',1,82));

% Pool de workers (Ryzen 5900X: 12 núcleos físicos). Cada sujeto es independiente.
if isempty(gcp('nocreate')); parpool('local', min(12, feature('numcores'))); end

for g = 1:2
    grp  = GROUPS{g};
    subs = SUBJ.(grp);
    nS   = NCAS.(grp);
    pk_v = PEAK.(grp);

    fprintf('\n>>> GRUPO: %s  (%d sujetos)\n', upper(grp), nS);

    for j = 1:nS
        suj  = subs{j};
        pk   = pk_v(j);
        t_sj = tic;

        fprintf('\n[%s %2d/%d] %-10s  peak_hz=%.3f Hz\n', grp, j, nS, suj, pk);

        if isnan(pk)
            fprintf('  [SKIP] Sin frecuencia válida.\n'); continue;
        end

        f_set = fullfile(P.dir_cont, [suj '_final.set']);
        if ~exist(f_set,'file')
            fprintf('  [SKIP] %s_final.set no encontrado.\n', suj); continue;
        end

        % ── Carga ──────────────────────────────────────────────────────────
        try
            EEG = pop_loadset('filename',[suj '_final.set'], 'filepath',P.dir_cont);
        catch ME
            fprintf('  [ERROR] Carga: %s\n', ME.message); continue;
        end
        fs = EEG.srate;
        fprintf('  Cargado: %d ch × %d pts | fs=%d Hz | %d eventos\n', ...
            EEG.nbchan, EEG.pnts, fs, numel(EEG.event));

        % ── Detectar canales por label ──────────────────────────────────────
        all_lbl = {EEG.chanlocs.labels};

        emg_idx = [];
        for lbl = P.emg_labels
            h = find(strcmpi(all_lbl, lbl{1}));
            if ~isempty(h) && ~ismember(h(1),emg_idx); emg_idx(end+1)=h(1); end
            if numel(emg_idx)==2; break; end
        end
        if numel(emg_idx)<2
            fb = P.emg_idx_fb(P.emg_idx_fb<=EEG.nbchan);
            emg_idx = fb(1:min(2,end));
            fprintf('  [AVISO] EMG fallback → ch[%s]\n', num2str(emg_idx));
        end

        roi_idx = [];
        for lbl = P.roi_labels
            h = find(strcmpi(all_lbl, lbl{1}));
            if ~isempty(h); roi_idx(end+1)=h(1); end
        end
        if isempty(roi_idx)
            roi_idx = P.roi_idx_fb(P.roi_idx_fb<=EEG.nbchan);
            fprintf('  [AVISO] ROI fallback → ch[%s]\n', num2str(roi_idx));
        end

        fprintf('  EMG: ch[%s]=%s | ROI: ch[%s]=%s\n', ...
            num2str(emg_idx), strjoin(all_lbl(emg_idx),' − '), ...
            num2str(roi_idx),  strjoin(all_lbl(roi_idx),' + '));

        % ── Señales continuas ───────────────────────────────────────────────
        if numel(emg_idx)>=2
            emg_cont = double(EEG.data(emg_idx(1),:)) - double(EEG.data(emg_idx(2),:));
        else
            emg_cont = double(EEG.data(emg_idx(1),:));
        end
        eeg_cont = double(mean(EEG.data(roi_idx,:), 1));

        emg_rms = rms(emg_cont);
        fprintf('  EMG RMS=%.2f µV\n', emg_rms);
        if emg_rms < 0.1
            fprintf('  [AVISO] EMG muy plano — verificar canales.\n');
        end

        % ── Filtrado continuo completo ──────────────────────────────────────
        fp_lo = max(0.5, pk - P.bw_pac);
        fp_hi = pk + P.bw_pac;
        [b_ph, a_ph] = butter(P.filt_order, [fp_lo fp_hi]/(fs/2), 'bandpass');
        phi_cont = angle(hilbert(filtfilt(b_ph, a_ph, emg_cont)));

        phi_var = var(phi_cont);
        fprintf('  Fase EMG [%.2f–%.2f Hz] | var=%.3f%s\n', fp_lo, fp_hi, phi_var, ...
            ternary(phi_var < 0.1, ' ← AVISO: fase plana', ''));

        amp_cont = nan(nF, EEG.pnts, 'single');
        for a = 1:nF
            fa = P.amp_freqs(a);
            [b_a,a_a] = butter(P.filt_order, [max(1,fa-2) fa+2]/(fs/2), 'bandpass');
            amp_cont(a,:) = single(abs(hilbert(filtfilt(b_a,a_a,eeg_cont))));
        end

        % ── Extraer latencias por condición ────────────────────────────────
        ev_types = pac_event_types(EEG);
        ev_lats  = round([EEG.event.latency]);
        lat_ch   = ev_lats(ev_types == P.trig_ch);
        lat_nc   = ev_lats(ev_types == P.trig_nc);
        fprintf('  Markers: %d×ch(40)  %d×nc(30)\n', numel(lat_ch), numel(lat_nc));

        if isempty(lat_ch) && isempty(lat_nc)
            fprintf('  [SKIP] Sin eventos.\n');
            clear EEG emg_cont eeg_cont phi_cont amp_cont; continue;
        end

        % ── PAC por condición ───────────────────────────────────────────────
        for c = 1:2
            cn   = CONDS{c};
            lats = eval(['lat_' lower(cn)]);
            n_tr = numel(lats);

            R.(grp).(cn).n_trials(j)  = n_tr;
            R.(grp).(cn).peak_used(j) = pk;

            if n_tr < 5
                fprintf('  [%s] %d trials → SKIP\n', cn, n_tr); continue;
            end
            fprintf('  [%s] %d trials → calculando PAC...\n', cn, n_tr);

            % Por ventana
            for w = 1:nW
                s0_samp = round(P.wins.ms(w,1)*fs/1000);
                s1_samp = round(P.wins.ms(w,2)*fs/1000);
                n_samp  = s1_samp - s0_samp + 1;

                phi_mat = nan(n_samp, n_tr, 'single');
                amp_mat = nan(nF, n_samp, n_tr, 'single');
                valid   = false(n_tr,1);

                for tr = 1:n_tr
                    i0 = lats(tr) + s0_samp;
                    i1 = lats(tr) + s1_samp;
                    if i0<1 || i1>EEG.pnts; continue; end
                    phi_mat(:,tr)    = phi_cont(i0:i1)';
                    amp_mat(:,:,tr)  = amp_cont(:,i0:i1);
                    valid(tr)        = true;
                end

                v_tr = find(valid);
                n_v  = numel(v_tr);
                R.(grp).(cn).n_valid(j,w) = n_v;

                if n_v < 5; continue; end

                phi_pool = reshape(phi_mat(:,v_tr), 1, []);

                % Temporales sliced por frecuencia (para parfor)
                pp_a = nan(nF,1); zmi_a = nan(nF,1); mi_a = nan(nF,1);
                mus_a = nan(nF,1); sds_a = nan(nF,1);
                n_bins = P.n_bins; n_surr = P.n_surr;   % broadcast escalares

                % parfor sobre frecuencias: cada una corre sus 200 surrogados en paralelo
                parfor a = 1:nF
                    amp_2d   = squeeze(amp_mat(a,:,v_tr));   % sliced input (dim1 = a)
                    amp_pool = reshape(amp_2d, 1, []);

                    v = ~isnan(phi_pool) & ~isnan(amp_pool) & ...
                        isfinite(phi_pool) & isfinite(amp_pool);
                    if sum(v) < 50; continue; end

                    phi_v = phi_pool(v);
                    amp_v = amp_pool(v);

                    MI_real = pac_mi(phi_v, amp_v, n_bins);
                    zv = mean(amp_v .* exp(1i*phi_v));
                    pp_a(a) = angle(zv);

                    % Surrogados: permutación de trials (RNG por worker → z-scores estables)
                    MI_surr = nan(n_surr,1);
                    for s = 1:n_surr
                        perm_tr  = v_tr(randperm(n_v));
                        phi_perm = reshape(phi_mat(:,perm_tr), 1, []);
                        vs = ~isnan(phi_perm) & ~isnan(amp_pool) & ...
                             isfinite(phi_perm) & isfinite(amp_pool);
                        if sum(vs)<50; continue; end
                        MI_surr(s) = pac_mi(phi_perm(vs), amp_pool(vs), n_bins);
                    end
                    mu_s = mean(MI_surr,'omitnan');
                    sd_s = std( MI_surr,'omitnan');

                    zmi_a(a) = (MI_real - mu_s)/(sd_s+eps);
                    mi_a(a)  = MI_real;
                    mus_a(a) = mu_s;
                    sds_a(a) = sd_s;
                end % frecuencias (parfor)

                % Ensamblar resultados de las frecuencias en R
                R.(grp).(cn).pref_phase(j,:,w)  = pp_a.';
                R.(grp).(cn).zMI(j,:,w)         = zmi_a.';
                R.(grp).(cn).MI(j,:,w)          = mi_a.';
                R.(grp).(cn).MI_surr_mu(j,:,w)  = mus_a.';
                R.(grp).(cn).MI_surr_sd(j,:,w)  = sds_a.';
            end % ventanas

            % Resumen consola
            th_m = P.amp_freqs>=P.bands_hz{1}(1) & P.amp_freqs<=P.bands_hz{1}(2);
            be_m = P.amp_freqs>=P.bands_hz{3}(1) & P.amp_freqs<=P.bands_hz{3}(2);
            fprintf('    [%s] zMI θ-Early=%+.2f | β-Early=%+.2f | β-Late=%+.2f\n', cn, ...
                mean(R.(grp).(cn).zMI(j,th_m,1),'omitnan'), ...
                mean(R.(grp).(cn).zMI(j,be_m,1),'omitnan'), ...
                mean(R.(grp).(cn).zMI(j,be_m,2),'omitnan'));
        end % condiciones

        clear EEG emg_cont eeg_cont phi_cont amp_cont phi_mat amp_mat;
        fprintf('  Tiempo: %.1f s\n', toc(t_sj));
    end % sujetos
end % grupos

fprintf('\n%s\n  LOOP FINALIZADO — %.1f min totales\n%s\n\n', ...
    repmat('-',1,82), toc/60, repmat('-',1,82));


%% ── 6. AGREGAR POR BANDA × VENTANA ─────────────────────────────────────────
get_band = @(zMI, b_hz) ...
    squeeze(mean(zMI(:, P.amp_freqs>=b_hz(1) & P.amp_freqs<=b_hz(2), :), 2, 'omitnan'));

B = struct();
for g = 1:2
    grp = GROUPS{g};
    for c = 1:2
        cn = CONDS{c};
        for bi = 1:3
            B.(grp).(cn).(P.bands{bi}) = get_band(R.(grp).(cn).zMI, P.bands_hz{bi});
            % B.(grp).(cn).(band): [nS × nW]
        end
    end
end


%% ── 7. GUARDAR WORKSPACE ────────────────────────────────────────────────────
out_mat = fullfile(P.dir_out, 'PAC_4Groups_Workspace.mat');
save(out_mat, 'R','B','P','casos','controles','peak_hz_cas','peak_med','SUBJ','NCAS', '-v7.3');
fprintf('>>> Workspace guardado:\n    %s\n\n', out_mat);


%% ── 8. CARGAR CONDUCTA (Casos) ──────────────────────────────────────────────
IES_Ch_cas = nan(nCas, 1);
if exist(P.f_behavior,'file')
    tmp_beh = load(P.f_behavior);
    BEH_Ch  = tmp_beh.tb_data_45.casos.chew;
    for i = 1:nCas
        idx = find(strcmp(string(BEH_Ch.Participantes), casos{i}));
        if ~isempty(idx); IES_Ch_cas(i) = BEH_Ch.ies_m(idx(1)); end
    end
    fprintf('>>> Conducta: %d/%d casos con IES_Ch\n\n', sum(~isnan(IES_Ch_cas)), nCas);
else
    fprintf('[AVISO] Archivo de conducta no encontrado. Correlaciones omitidas.\n\n');
end


%% ── 9. ESTADÍSTICA Y REPORTE ────────────────────────────────────────────────
f_rep = fullfile(P.dir_reports, 'Reporte_PAC_4Groups.txt');
fid   = fopen(f_rep, 'w');
W = @(s) (fprintf('%s\n',s) + fprintf(fid,'%s\n',s));

W(repmat('=',1,90));
W('  REPORTE PAC 4 GRUPOS — S6_PAC_4Groups.m');
W(sprintf('  Generado: %s', datestr(now)));
W(repmat('=',1,90));
W(sprintf('  Casos: N=%d | Controles: N=%d', nCas, nCtr));
W(sprintf('  Excluidos: %s', strjoin(P.excluded,', ')));
W(sprintf('  Peak EMG Casos: mediana=%.3f Hz | Controles (fijo): %.3f Hz', peak_med, peak_med));
W(sprintf('  Surrogados: permutación trials (N=%d)', P.n_surr));
W(sprintf('  ROI EEG: %s', strjoin(P.roi_labels,' + ')));
W('');

% Header común para tablas banda×ventana
hdr = sprintf('  %-22s', '');
for bi = 1:3
    for wi = 1:nW
        hdr = [hdr sprintf(' %-13s', [P.bands{bi}(1:3) '-' P.wins.labels{wi}(1:3)])]; %#ok<AGROW>
    end
end

% ────────────────────────────────────────────────────────────────────────────
W('Q1 — ¿PAC > surrogate? (signrank vs 0, one-tailed)');
W('     Mediana zMI  (* p<.05  ** p<.01  *** p<.001)');
W(hdr);
for g = 1:2; for c = 1:2
    grp = GROUPS{g}; cn = CONDS{c};
    row = sprintf('  %-22s', [grp '/' cn]);
    for bi = 1:3; for wi = 1:nW
        v = B.(grp).(cn).(P.bands{bi})(:,wi);
        v = v(~isnan(v));
        if numel(v)>=5
            p = signrank(v, 0, 'tail','right');
            row = [row sprintf(' %+.3f%-4s    ', median(v), stars(p))]; %#ok<AGROW>
        else
            row = [row sprintf(' NaN          ')]; %#ok<AGROW>
        end
    end; end
    W(row);
end; end
W('');

% ────────────────────────────────────────────────────────────────────────────
W('Q2 — Delta zMI (Ch − Nc), signrank pareado DENTRO DE CADA GRUPO');
W('     Interpretación: solo significativo en Casos indica especificidad masticatoria.');
W(hdr);
for g = 1:2
    grp = GROUPS{g};
    row_m = sprintf('  %-22s', ['Med Δ ' grp]);
    row_p = sprintf('  %-22s', ['p   Δ ' grp]);
    for bi = 1:3; for wi = 1:nW
        ch = B.(grp).Ch.(P.bands{bi})(:,wi);
        nc = B.(grp).Nc.(P.bands{bi})(:,wi);
        v  = ~isnan(ch) & ~isnan(nc);
        if sum(v)>=5
            d = ch(v)-nc(v);
            p = signrank(ch(v),nc(v));
            row_m = [row_m sprintf(' %+.3f        ', median(d))];      %#ok<AGROW>
            row_p = [row_p  sprintf(' %.4f%-4s    ', p, stars(p))];    %#ok<AGROW>
        else
            row_m = [row_m sprintf(' NaN          ')]; %#ok<AGROW>
            row_p = [row_p  sprintf(' NaN          ')]; %#ok<AGROW>
        end
    end; end
    W(row_m); W(row_p);
end
W('');

% ────────────────────────────────────────────────────────────────────────────
W('Q3 — Banda preferida: Theta vs Alpha vs Beta (Friedman, ventana Early)');
W(sprintf('  %-22s  %-6s  %-8s  %-8s  %-8s  %-14s', ...
    'Grupo/Cond', 'N', 'Med-θ', 'Med-α', 'Med-β', 'Friedman-p'));
for g = 1:2; for c = 1:2
    grp = GROUPS{g}; cn = CONDS{c};
    Th = B.(grp).(cn).Theta(:,1);
    Al = B.(grp).(cn).Alpha(:,1);
    Be = B.(grp).(cn).Beta(:,1);
    M  = [Th Al Be];
    v  = all(~isnan(M),2);
    if sum(v)>=5
        Mv = M(v,:);
        p_fr = friedman(Mv,1,'off');
        W(sprintf('  %-22s  %-6d  %-8.3f  %-8.3f  %-8.3f  %-14.4g%s', ...
            [grp '/' cn], sum(v), median(Mv(:,1)), median(Mv(:,2)), median(Mv(:,3)), ...
            p_fr, ternary(p_fr<0.05,' *','')));
    end
end; end
W('');

% ────────────────────────────────────────────────────────────────────────────
W('Q4 — Test de Rayleigh (fase preferida EMG): TODOS LOS GRUPOS × BANDAS × VENTANAS');
W('     R = mean vector length [0-1]. Sig solo en Casos/Ch indica especificidad masticatoria.');
W('     mu = ángulo preferido en grados (0°=cresta masticatoria, 180°=valle)');
W('');

band_masks = {P.amp_freqs>=P.bands_hz{1}(1) & P.amp_freqs<=P.bands_hz{1}(2), ...
              P.amp_freqs>=P.bands_hz{2}(1) & P.amp_freqs<=P.bands_hz{2}(2), ...
              P.amp_freqs>=P.bands_hz{3}(1) & P.amp_freqs<=P.bands_hz{3}(2)};

for bi = 1:3
    W(sprintf('  === BANDA: %s ===', P.bands{bi}));
    W(sprintf('  %-18s  %5s %7s %5s %8s    %5s %7s %5s %8s    %5s %7s %5s %8s    %5s %7s %5s %8s', ...
        'Ventana', ...
        'N','mu°','R','p [CasCh]', ...
        'N','mu°','R','p [CasNc]', ...
        'N','mu°','R','p [CtrCh]', ...
        'N','mu°','R','p [CtrNc]'));
    W(sprintf('  %s', repmat('-',1,96)));

    for wi = 1:nW
        row = sprintf('  %-18s', P.wins.labels{wi});
        for g = 1:2; for c = 1:2
            grp = GROUPS{g}; cn = CONDS{c};
            phi = squeeze(mean(R.(grp).(cn).pref_phase(:, band_masks{bi}, wi), 2, 'omitnan'));
            phi = phi(~isnan(phi));
            if numel(phi)>=5
                ray = pac_rayleigh(phi);
                row = [row sprintf('  %5d %7.1f %5.3f %7.4g%-3s', ...
                    numel(phi), rad2deg(ray.mu), ray.R, ray.p, stars(ray.p))]; %#ok<AGROW>
            else
                row = [row sprintf('  %5s %7s %5s %10s', 'NA','—','—','—')]; %#ok<AGROW>
            end
        end; end
        W(row);
    end
    W('');
end

% ────────────────────────────────────────────────────────────────────────────
W('Q5 — Selectividad temporal (Friedman Early/Late/Active, banda Theta)');
W(sprintf('  %-22s %-6s %-10s %-10s %-10s %-14s %s', ...
    'Grupo/Cond','N','Med-Early','Med-Late','Med-Active','Friedman-p','E<>L / E<>A / L<>A'));
for g = 1:2; for c = 1:2
    grp = GROUPS{g}; cn = CONDS{c};
    Th = B.(grp).(cn).Theta;
    v  = all(~isnan(Th),2);
    if sum(v)>=5
        p_fr = friedman(Th(v,:),1,'off');
        p_EL = signrank(Th(v,1),Th(v,2));
        p_EA = signrank(Th(v,1),Th(v,3));
        p_LA = signrank(Th(v,2),Th(v,3));
        W(sprintf('  %-22s %-6d %-10.3f %-10.3f %-10.3f %-14.4g %.3f/%.3f/%.3f%s', ...
            [grp '/' cn], sum(v), ...
            median(Th(v,1)), median(Th(v,2)), median(Th(v,3)), ...
            p_fr, p_EL, p_EA, p_LA, ternary(p_fr<0.05,' *','')));
    end
end; end
W('');

% ────────────────────────────────────────────────────────────────────────────
W('Q6 — Comparaciones ENTRE GRUPOS (ranksum, independiente)');
W('     Interpretación: Casos/Ch vs Controles/Ch = especificidad patológica del PAC');
W(hdr);

contrasts = {
    'CasCh vs CtrCh', 'Casos','Ch',  'Controles','Ch';
    'CasNc vs CtrNc', 'Casos','Nc',  'Controles','Nc';
    'CasCh vs CasNc', 'Casos','Ch',  'Casos','Nc';    % para referencia
    'CtrCh vs CtrNc', 'Controles','Ch','Controles','Nc';
};

for k = 1:size(contrasts,1)
    lbl  = contrasts{k,1};
    g1   = contrasts{k,2}; c1 = contrasts{k,3};
    g2   = contrasts{k,4}; c2 = contrasts{k,5};
    row  = sprintf('  %-22s', lbl);
    for bi = 1:3; for wi = 1:nW
        d1 = B.(g1).(c1).(P.bands{bi})(:,wi); d1 = d1(~isnan(d1));
        d2 = B.(g2).(c2).(P.bands{bi})(:,wi); d2 = d2(~isnan(d2));
        if numel(d1)>=5 && numel(d2)>=5
            if strcmp(g1,g2)   % mismo grupo → signrank pareado
                v = ~isnan(B.(g1).(c1).(P.bands{bi})(:,wi)) & ~isnan(B.(g2).(c2).(P.bands{bi})(:,wi));
                if sum(v)>=5
                    p = signrank(B.(g1).(c1).(P.bands{bi})(v,wi), B.(g2).(c2).(P.bands{bi})(v,wi));
                else; p = NaN;
                end
            else               % grupos distintos → ranksum
                p = ranksum(d1, d2);
            end
            row = [row sprintf(' %.4f%-4s    ', p, stars(p))]; %#ok<AGROW>
        else
            row = [row sprintf(' NaN          ')]; %#ok<AGROW>
        end
    end; end
    W(row);
end
W('');

% ────────────────────────────────────────────────────────────────────────────
W('Q7 — Correlaciones zMI_Ch vs IES_Ch (Spearman, FDR-BH, solo Casos)');
if any(~isnan(IES_Ch_cas))
    rho_beh = nan(3,nW); p_beh = nan(3,nW);
    for bi = 1:3; for wi = 1:nW
        z = B.Casos.Ch.(P.bands{bi})(:,wi);
        v = ~isnan(z) & ~isnan(IES_Ch_cas);
        if sum(v)>=5
            [rho_beh(bi,wi), p_beh(bi,wi)] = corr(z(v), IES_Ch_cas(v), 'Type','Spearman');
        end
    end; end
    % FDR-BH
    p_flat = p_beh(:);
    [p_s, idx_s] = sort(p_flat); m = numel(p_flat);
    p_fdr = nan(m,1);
    for k = 1:m
        if ~isnan(p_s(k)); p_fdr(idx_s(k)) = min(1, p_s(k)*m/k); end
    end
    p_fdr = reshape(p_fdr,3,nW);

    W(sprintf('  %-8s  %-8s  %7s  %8s  %8s', 'Banda','Ventana','rho','p_raw','p_FDR'));
    W(sprintf('  %s', repmat('-',1,48)));
    for bi = 1:3; for wi = 1:nW
        r=rho_beh(bi,wi); pr=p_beh(bi,wi); pf=p_fdr(bi,wi);
        if isnan(r); continue; end
        sym=''; if pf<0.05; sym='* FDR'; elseif pr<0.05; sym='* raw'; end
        W(sprintf('  %-8s  %-8s  %+7.3f  %8.4f  %8.4f   %s', ...
            P.bands{bi}, P.wins.labels{wi}, r, pr, pf, sym));
    end; end
else
    W('  [OMITIDO] Conducta no disponible.');
end
W('');

W(repmat('=',1,90));
W(sprintf('  Tiempo total: %.1f min', toc/60));
fclose(fid);
fprintf('>>> Reporte: %s\n\n', f_rep);


%% ── 10. FIGURAS RAYLEIGH (rose plots) ───────────────────────────────────────
% Una figura 4×3 (grupos×conds × bandas) para cada ventana temporal
colors4 = [0.20 0.55 0.40;   % Casos/Ch  verde
           0.65 0.80 0.70;   % Casos/Nc  verde pálido
           0.20 0.40 0.70;   % Controles/Ch  azul
           0.60 0.75 0.90];  % Controles/Nc  azul pálido

grp_cond_labels = {'Casos/Ch','Casos/Nc','Controles/Ch','Controles/Nc'};
GC_pairs = {{'Casos','Ch'},{'Casos','Nc'},{'Controles','Ch'},{'Controles','Nc'}};

for wi = 1:nW
    fig = figure('Color','w','Position',[50 50 1200 900],'Visible','off');
    tl  = tiledlayout(fig, 4, 3, 'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Rose Plots — Fase preferida EMG | Ventana: %s', P.wins.labels{wi}), ...
        'FontSize',13,'FontWeight','bold');

    for row = 1:4
        grp = GC_pairs{row}{1};
        cn  = GC_pairs{row}{2};
        for bi = 1:3
            ax = nexttile;
            phi = squeeze(mean(R.(grp).(cn).pref_phase(:, band_masks{bi}, wi), 2,'omitnan'));
            phi = phi(~isnan(phi));

            if numel(phi)>=5
                % Rose plot manual (compatible con todas las versiones MATLAB)
                n_bins_rose = 12;
                edges_rose  = linspace(-pi, pi, n_bins_rose+1);
                counts = histcounts(phi, edges_rose);
                bin_c  = (edges_rose(1:end-1)+edges_rose(2:end))/2;

                % Dibujar sectores usando patch
                hold(ax,'on');
                r_max = max(counts)*1.2 + 0.5;
                for b = 1:n_bins_rose
                    if counts(b)==0; continue; end
                    th1 = edges_rose(b); th2 = edges_rose(b+1);
                    r   = counts(b);
                    th_arc = linspace(th1,th2,20);
                    xp = [0, r*sin(th_arc), 0];
                    yp = [0, r*cos(th_arc), 0];
                    patch(ax, xp, yp, colors4(row,:), 'EdgeColor','w', ...
                        'FaceAlpha',0.75, 'LineWidth',0.5);
                end

                % Vector resultante
                ray = pac_rayleigh(phi);
                r_plot = r_max * 0.85 * ray.R;
                quiver(ax, 0, 0, r_plot*sin(ray.mu), r_plot*cos(ray.mu), 0, ...
                    'k', 'LineWidth',2.5, 'MaxHeadSize', 0.4);

                % Ejes
                plot(ax, [-r_max r_max], [0 0], 'k-', 'LineWidth',0.5, 'Color',[0.6 0.6 0.6]);
                plot(ax, [0 0], [-r_max r_max], 'k-', 'LineWidth',0.5, 'Color',[0.6 0.6 0.6]);
                viscircles(ax, [0 0], r_max*0.95, 'Color',[0.7 0.7 0.7], 'LineWidth',0.5, 'EnhanceVisibility',false);
                axis(ax, [-r_max r_max -r_max r_max]); axis(ax,'square','off');

                % Etiquetas de fase
                text(ax,  r_max*0.05,  r_max*0.95, '0',   'FontSize',8,'HorizontalAlignment','center');
                text(ax, -r_max*0.95,  r_max*0.05, 'π',   'FontSize',9,'HorizontalAlignment','center');
                text(ax,  r_max*0.95,  r_max*0.05, '0',   'FontSize',8,'HorizontalAlignment','center');

                p_sym = stars(ray.p);
                title(ax, sprintf('%s | %s\nR=%.3f  μ=%.0f°  %s', ...
                    grp_cond_labels{row}, P.bands{bi}, ...
                    ray.R, rad2deg(ray.mu), ...
                    ternary(isempty(p_sym),'n.s.',p_sym)), ...
                    'FontSize',8,'FontWeight','bold');
            else
                title(ax, sprintf('%s | %s\nN<5', grp_cond_labels{row}, P.bands{bi}), ...
                    'FontSize',8);
                axis(ax,'off');
            end
        end % bandas
    end % grupos×conds

    fn_fig = fullfile(P.dir_plots, sprintf('Fig_Rayleigh_%s.png', P.wins.labels{wi}));
    exportgraphics(fig, fn_fig, 'Resolution',300);
    close(fig);
    fprintf('>>> Rose plots (%s): %s\n', P.wins.labels{wi}, fn_fig);
end

fprintf('\n>>> S6_PAC_4Groups DONE — %.1f min totales.\n\n', toc/60);


%% ============================================================================
%  FUNCIONES LOCALES
%% ============================================================================

function mi = pac_mi(phi, amp, n_bins)
edges    = linspace(-pi, pi, n_bins+1);
mean_amp = zeros(1, n_bins);
for b = 1:n_bins
    in_b = phi >= edges(b) & phi < edges(b+1);
    if any(in_b); mean_amp(b) = mean(amp(in_b)); end
end
mean_amp(mean_amp<=0) = eps;
p  = mean_amp / sum(mean_amp);
H  = -sum(p .* log(p+eps));
mi = (log(n_bins) - H) / log(n_bins);
end

function ev_num = pac_event_types(EEG)
ev_num = nan(1, numel(EEG.event));
for k = 1:numel(EEG.event)
    t = EEG.event(k).type;
    if isnumeric(t);                ev_num(k) = t;
    elseif ischar(t)||isstring(t);  n=str2double(t); if ~isnan(n); ev_num(k)=n; end
    end
end
end

function ray = pac_rayleigh(phi)
phi = phi(:); n = numel(phi);
Rg  = abs(mean(exp(1i*phi)));
mu  = angle(mean(exp(1i*phi)));
Z   = n * Rg^2;
p   = exp(-Z) * (1 + (2*Z-Z^2)/(4*n) - (24*Z-132*Z^2+76*Z^3-9*Z^4)/(288*n^2));
p   = max(0, min(1,p));
ray = struct('R',Rg,'mu',mu,'Z',Z,'p',p,'n',n);
end

function s = stars(p)
if     isnan(p);  s = '';
elseif p < 0.001; s = '***';
elseif p < 0.01;  s = '**';
elseif p < 0.05;  s = '*';
else;             s = '';
end
end

function r = ternary(cond, a, b)
if cond; r = a; else; r = b; end
end