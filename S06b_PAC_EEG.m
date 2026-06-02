%% ============================================================================
%  S5c_PAC_EEG.m — PAC EEG-EEG: fase theta (F1/FC1) × amplitud (8-100 Hz)
%  ----------------------------------------------------------------------------
%  MECANISMO QUE SE TESTA:
%    La hipótesis de Lisman & Jensen (2005) / Canolty & Knight (2010):
%    el theta frontal organiza paquetes de actividad de alta frecuencia —
%    cada ciclo theta "abre una ventana" para codificar/mantener ítems en WM.
%    Si la masticación reorganiza la fase theta frontal, se espera:
%      1. Aumento del PAC theta-gamma en condición Ch vs Nc/surrogates.
%      2. zMI theta-gamma predice IES individualmente (Ch correlaciona con WM).
%
%  DIFERENCIA CLAVE con S5b_PAC_Continuous.m:
%    - Driver de fase: EEG theta (F1/FC1, 4-7 Hz), NO el EMG.
%    - Esto mide reorganización COGNITIVA de la fase theta frontal.
%    - S5b medía coordinación sensoriomotora (EMG→beta cortical).
%
%  PIPELINE:
%    1. Cargar archivo continuo S5_Final por sujeto.
%    2. Filtrar señal COMPLETA en theta (fase) y múltiples freqs (amplitud).
%    3. Extraer ventanas desde señal filtrada alrededor de marcas 40/30.
%    4. Calcular MI (Tort 2010) + surrogates permutación de trials.
%    5. PARFOR sobre sujetos → Ryzen 5900X 12 cores.
%
%  OUTPUT:
%    EEG/PAC/EEG_EEG/PAC_EEG_Workspace.mat
%    EEG/PAC/EEG_EEG/Reports/Reporte_PAC_EEG.txt
%    EEG/PAC/EEG_EEG/Plots/*.png
%
%  Sebastián, 2026-05-19
%  ============================================================================
clear; clc; close all; rng(2026); tic;

%% ── 1. CONFIGURACIÓN ────────────────────────────────────────────────────────
P = struct();

S = S0_paths();                    % rutas + parámetros UNIFICADOS
P.dir_cont    = S.dir_cont;        % *_final.set (referenciado, solo lectura)
P.dir_lists   = S.dir_data;        % incluidos45.mat (copiado al proyecto)
P.dir_out     = S.fig04;           % PAC_EEG_Workspace.mat → proyecto/Figure04_PAC
P.dir_reports = S.fig04;
P.dir_plots   = S.fig04;
P.eeglab_path = S.eeglab_path;
P.f_behavior  = S.file_beh;
P.f_lists     = S.file_incluidos;

% ── Canales EEG ROI (por label; fallback por índice) ────────────────────────
P.phase_labels  = {'F1','FC1'};      % driver de FASE theta frontal
P.phase_idx_fb  = [12, 16];
P.amp_labels    = {'F1','FC1'};      % fuente de AMPLITUD (mismo ROI, intra-areal)
P.amp_idx_fb    = [12, 16];

% ── Parámetros PAC ──────────────────────────────────────────────────────────
P.phase_band    = S.bands_hz{1};     % theta desde config (4-7 Hz)
P.amp_freqs     = 8:4:100;          % espectro de amplitud (24 frecuencias)
P.n_bins        = 18;               % bins MI (Tort 2010)
P.n_surr        = 200;             % surrogates por permutación de trials
P.filt_order    = 4;               % orden Butterworth

% ── Ventanas y marcadores ────────────────────────────────────────────────────
P.wins.labels   = {'Early','Late','Active'};
P.wins.ms       = [0 300; 300 900; 200 700];
P.trig_ch       = 40;
P.trig_nc       = 30;
P.run_nc_verification = true;

% ── Paralelo ────────────────────────────────────────────────────────────────
% Ryzen 5900X: 12 cores físicos / 24 threads
% Se usan los 12 cores físicos (1 worker por core = máximo throughput sin
% hyperthreading overhead en operaciones FP intensivas como filtfilt/hilbert)
P.n_workers = 12;

% ── Sujetos excluidos (desde config) ─────────────────────────────────────────
P.excluded = S.excluded;

% ── Crear directorios ────────────────────────────────────────────────────────
for d = {P.dir_out, P.dir_reports, P.dir_plots}
    if ~exist(d{1},'dir'); mkdir(d{1}); end
end

% ── EEGLAB al path (una vez, antes del parfor) ───────────────────────────────
if exist(P.eeglab_path,'dir')
    addpath(P.eeglab_path);
    addpath(genpath(fullfile(P.eeglab_path,'functions')));
    eeglab nogui;
end

% ── Encabezado consola ───────────────────────────────────────────────────────
sep = repmat('=',1,82);
fprintf('\n%s\n', sep);
fprintf('  S5c_PAC_EEG.m — PAC EEG-EEG (theta-phase × broadband-amplitude)\n');
fprintf('  %s\n', datestr(now));
fprintf('%s\n', sep);
fprintf('  Fase        : theta [%d–%d Hz] en F1/FC1 — señal continua filtrada\n', P.phase_band);
fprintf('  Amplitud    : %d–%d Hz en pasos de 4 Hz (%d frecuencias)\n', ...
    P.amp_freqs(1), P.amp_freqs(end), numel(P.amp_freqs));
fprintf('  Surrogates  : permutación de trials (N=%d)\n', P.n_surr);
fprintf('  Paralelo    : %d workers (Ryzen 5900X)\n', P.n_workers);
fprintf('%s\n\n', repmat('-',1,82));


%% ── 2. SUJETOS (CASOS) ──────────────────────────────────────────────────────
tmp      = load(P.f_lists);
fn       = fieldnames(tmp); incluidos = tmp.(fn{1})(:);
casos    = sort(incluidos(startsWith(incluidos, 'E3S')));
casos    = setdiff(casos, P.excluded, 'stable');
nS       = numel(casos);
fprintf('  Casos a procesar : %d\n', nS);
fprintf('  Excluidos        : %s\n\n', strjoin(P.excluded, ', '));


%% ── 3. CARGAR CONDUCTA ──────────────────────────────────────────────────────
IES_Ch = nan(nS,1); IES_Nc = nan(nS,1);
if exist(P.f_behavior,'file')
    tmp_beh  = load(P.f_behavior);
    BEH_Ch   = tmp_beh.tb_data_45.casos.chew;
    BEH_Nc   = tmp_beh.tb_data_45.casos.nochew;
    for i = 1:nS
        suj = casos{i};
        idx = find(strcmp(string(BEH_Ch.Participantes), suj));
        if ~isempty(idx); IES_Ch(i) = BEH_Ch.ies_m(idx(1)); end
        idx = find(strcmp(string(BEH_Nc.Participantes), suj));
        if ~isempty(idx); IES_Nc(i) = BEH_Nc.ies_m(idx(1)); end
    end
    fprintf('  Conducta cargada : %d/%d sujetos con IES_Ch (mediana=%.0f ms)\n\n', ...
        sum(~isnan(IES_Ch)), nS, median(IES_Ch,'omitnan'));
else
    fprintf('  [AVISO] data_beh_tb_45.mat no encontrado.\n\n');
end
Delta_IES = IES_Ch - IES_Nc;


%% ── 4. PRE-ALLOCACIÓN (arrays para parfor) ───────────────────────────────────
nF = numel(P.amp_freqs);
nW = size(P.wins.ms,1);
CONDS = {'Ch','Nc'};

% Arrays de salida: indexados por [sujeto × freq × ventana]
% parfor requiere arrays indexables, no structs con campos dinámicos
zMI_ch      = nan(nS, nF, nW, 'single');
zMI_nc      = nan(nS, nF, nW, 'single');
MI_ch       = nan(nS, nF, nW, 'single');
MI_nc       = nan(nS, nF, nW, 'single');
pref_ch     = nan(nS, nF, nW, 'single');  % fase preferida (rad)
pref_nc     = nan(nS, nF, nW, 'single');
nv_ch       = zeros(nS, nW, 'uint8');
nv_nc       = zeros(nS, nW, 'uint8');
qc_phi_var  = nan(nS, 1);   % varianza fase theta (QC)
qc_amp_rms  = nan(nS, nF);  % RMS amplitud por frecuencia (QC)
qc_ch_label = repmat({'?'}, nS, 1);  % canales detectados
qc_fs       = zeros(nS, 1);
qc_n_ch     = zeros(nS, 1);
qc_n_nc     = zeros(nS, 1);
failed      = false(nS, 1);

% Copias locales de parámetros para broadcast en parfor
p_phase_band   = P.phase_band;
p_amp_freqs    = P.amp_freqs;
p_wins_ms      = P.wins.ms;
p_n_bins       = P.n_bins;
p_n_surr       = P.n_surr;
p_filt_order   = P.filt_order;
p_trig_ch      = P.trig_ch;
p_trig_nc      = P.trig_nc;
p_run_nc       = P.run_nc_verification;
p_phase_labels = P.phase_labels;
p_phase_fb     = P.phase_idx_fb;
p_amp_labels   = P.amp_labels;
p_amp_fb       = P.amp_idx_fb;
p_dir_cont     = P.dir_cont;
p_wins_labels  = P.wins.labels;


%% ── 5. PARFOR PRINCIPAL ─────────────────────────────────────────────────────
fprintf('%s\n  INICIANDO PARFOR — %d sujetos × %d workers\n%s\n\n', ...
    repmat('-',1,82), nS, P.n_workers, repmat('-',1,82));

% Inicializar pool (reutiliza si ya existe)
pool = gcp('nocreate');
if isempty(pool) || pool.NumWorkers ~= P.n_workers
    if ~isempty(pool); delete(pool); end
    parpool('local', P.n_workers);
end

parfor j = 1:nS  %#ok<*PFBNS>
    suj = casos{j};

    % ── Cargar archivo continuo ────────────────────────────────────────────
    f_set = fullfile(p_dir_cont, [suj '_final.set']);
    if ~exist(f_set,'file')
        fprintf('  [%s] SKIP — archivo no encontrado\n', suj);
        failed(j) = true; %#ok<PFOUS>
        continue
    end

    EEG = [];
    try
        EEG = pop_loadset('filename',[suj '_final.set'], 'filepath',p_dir_cont);
    catch ME
        fprintf('  [%s] ERROR carga: %s\n', suj, ME.message);
        failed(j) = true;
        continue
    end

    fs      = EEG.srate;
    nPnts   = EEG.pnts;

    % ── Detectar canales por label ─────────────────────────────────────────
    if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs,'labels')
        all_lbl = {EEG.chanlocs.labels};
    else
        all_lbl = arrayfun(@(k) sprintf('ch%d',k), 1:EEG.nbchan, 'Uni',false);
    end

    % Fase: canales theta — índice numérico (requerido en parfor)
    phi_idx = zeros(1, numel(p_phase_labels));
    phi_count = 0;
    for ki = 1:numel(p_phase_labels)
        h = find(strcmpi(all_lbl, p_phase_labels{ki}));
        if ~isempty(h)
            phi_count = phi_count + 1;
            phi_idx(phi_count) = h(1);
        end
        if phi_count == 2; break; end
    end
    phi_idx = phi_idx(1:phi_count);
    if isempty(phi_idx)
        fb = p_phase_fb(p_phase_fb <= EEG.nbchan);
        phi_idx = fb(1:min(2,numel(fb)));
    end

    % Amplitud: mismos canales (índice numérico)
    amp_idx = zeros(1, numel(p_amp_labels));
    amp_count = 0;
    for ki = 1:numel(p_amp_labels)
        h = find(strcmpi(all_lbl, p_amp_labels{ki}));
        if ~isempty(h)
            amp_count = amp_count + 1;
            amp_idx(amp_count) = h(1);
        end
        if amp_count == 2; break; end
    end
    amp_idx = amp_idx(1:amp_count);
    if isempty(amp_idx)
        fb = p_amp_fb(p_amp_fb <= EEG.nbchan);
        amp_idx = fb(1:min(2,numel(fb)));
    end

    ch_str = strjoin(all_lbl(phi_idx), '+');
    fprintf('  [%s] fs=%d Hz | %d pts | Canales: %s\n', suj, fs, nPnts, ch_str);

    % ── Señal media sobre canales ROI ─────────────────────────────────────
    sig_phase = double(mean(EEG.data(phi_idx,:), 1));   % [1 × nPnts]
    sig_amp   = double(mean(EEG.data(amp_idx,:),  1));   % [1 × nPnts]

    % Extraer eventos y liberar EEG asignando a variable vacía (clear no válido en parfor)
    ev_types_raw = {EEG.event.type};
    ev_lats_raw  = round([EEG.event.latency]);
    EEG = [];   % liberar memoria sin usar clear

    % ── Convertir tipos de evento a numérico ─────────────────────────────
    ev_num = nan(1, numel(ev_types_raw));
    for k = 1:numel(ev_types_raw)
        t = ev_types_raw{k};
        if isnumeric(t)
            ev_num(k) = t;
        elseif ischar(t) || isstring(t)
            n = str2double(t);
            if ~isnan(n); ev_num(k) = n; end
        end
    end
    lat_ch = ev_lats_raw(ev_num == p_trig_ch);
    lat_nc = ev_lats_raw(ev_num == p_trig_nc);
    fprintf('    Markers: %d × Ch40  |  %d × Nc30\n', numel(lat_ch), numel(lat_nc));

    % ── Filtrar señal completa — FASE theta ──────────────────────────────
    fp_lo = p_phase_band(1); fp_hi = p_phase_band(2);
    [b_ph, a_ph] = butter(p_filt_order, [fp_lo fp_hi]/(fs/2), 'bandpass');
    phi_filt   = filtfilt(b_ph, a_ph, sig_phase);
    phi_cont   = angle(hilbert(phi_filt));          % [1 × nPnts]

    phi_var = var(phi_cont);
    fprintf('    Varianza fase theta: %.4f\n', phi_var);

    % ── Filtrar señal completa — AMPLITUD (todas las frecuencias) ────────
    nFreqs    = numel(p_amp_freqs);
    amp_cont  = nan(nFreqs, nPnts, 'single');

    for a = 1:nFreqs
        fa   = p_amp_freqs(a);
        f_lo = max(1, fa - 2);
        f_hi = min(fs/2 - 1, fa + 2);
        [b_a, a_a] = butter(p_filt_order, [f_lo f_hi]/(fs/2), 'bandpass');
        amp_cont(a,:) = single(abs(hilbert(filtfilt(b_a, a_a, sig_amp))));
    end

    % ── Extracción de ventanas y PAC ─────────────────────────────────────
    nW_loc = size(p_wins_ms, 1);

    % Variables locales de salida para este sujeto
    loc_zMI_ch   = nan(nFreqs, nW_loc, 'single');
    loc_zMI_nc   = nan(nFreqs, nW_loc, 'single');
    loc_MI_ch    = nan(nFreqs, nW_loc, 'single');
    loc_MI_nc    = nan(nFreqs, nW_loc, 'single');
    loc_pref_ch  = nan(nFreqs, nW_loc, 'single');
    loc_pref_nc  = nan(nFreqs, nW_loc, 'single');
    loc_nv_ch    = zeros(1, nW_loc, 'uint8');
    loc_nv_nc    = zeros(1, nW_loc, 'uint8');

    for c = 1:2
        if c==2 && ~p_run_nc
            continue
        end
        if c==1; lats=lat_ch; else; lats=lat_nc; end
        n_tr = numel(lats);
        if n_tr < 5; continue; end

        for w = 1:nW_loc
            win_s = round(p_wins_ms(w,1)*fs/1000);
            win_e = round(p_wins_ms(w,2)*fs/1000);
            n_smp = win_e - win_s + 1;

            phi_mat = nan(n_smp, n_tr, 'single');
            amp_mat = nan(nFreqs, n_smp, n_tr, 'single');
            valid   = false(n_tr, 1);

            for tr = 1:n_tr
                s0 = lats(tr) + win_s;
                s1 = lats(tr) + win_e;
                if s0 < 1 || s1 > nPnts; continue; end
                phi_mat(:,tr)     = phi_cont(s0:s1)';
                amp_mat(:,:,tr)   = amp_cont(:,s0:s1);
                valid(tr)         = true;
            end

            v_tr = find(valid);
            n_v  = numel(v_tr);

            if c==1; loc_nv_ch(w)=uint8(min(n_v,255));
            else;    loc_nv_nc(w)=uint8(min(n_v,255)); end

            if n_v < 5; continue; end

            phi_pool = reshape(phi_mat(:,v_tr), 1, []);

            for a = 1:nFreqs
                amp_2d   = squeeze(amp_mat(a,:,v_tr));   % [n_smp × n_v]
                amp_pool = reshape(amp_2d, 1, []);

                ok = ~isnan(phi_pool) & ~isnan(amp_pool) & ...
                     isfinite(phi_pool) & isfinite(amp_pool);
                if sum(ok) < 50; continue; end

                phi_v = phi_pool(ok);
                amp_v = amp_pool(ok);

                MI_real = pac_mi_local(phi_v, amp_v, p_n_bins);

                % Fase preferida
                pref_angle = single(angle(mean(amp_v .* exp(1i*phi_v))));

                % Surrogates: permutación de trials
                MI_surr = nan(p_n_surr, 1);
                for s = 1:p_n_surr
                    perm_tr   = v_tr(randperm(n_v));
                    phi_perm  = reshape(phi_mat(:,perm_tr), 1, []);
                    vp = ~isnan(phi_perm) & ~isnan(amp_pool) & ...
                         isfinite(phi_perm) & isfinite(amp_pool);
                    if sum(vp) < 50; continue; end
                    MI_surr(s) = pac_mi_local(phi_perm(vp), amp_pool(vp), p_n_bins);
                end

                mu_s = mean(MI_surr,'omitnan');
                sd_s = std(MI_surr,0,'omitnan');
                zmi  = single((MI_real - mu_s) / (sd_s + eps));

                if c==1
                    loc_zMI_ch(a,w)  = zmi;
                    loc_MI_ch(a,w)   = single(MI_real);
                    loc_pref_ch(a,w) = pref_angle;
                else
                    loc_zMI_nc(a,w)  = zmi;
                    loc_MI_nc(a,w)   = single(MI_real);
                    loc_pref_nc(a,w) = pref_angle;
                end
            end % frecuencias
        end % ventanas
    end % condiciones

    % Resumen en consola por sujeto
    th_m = p_amp_freqs >= 4  & p_amp_freqs <= 7;
    be_m = p_amp_freqs >= 14 & p_amp_freqs <= 30;
    ga_m = p_amp_freqs >= 30 & p_amp_freqs <= 80;
    fprintf('    [Ch] θ-Late=%+.2f | β-Late=%+.2f | γ-Late=%+.2f\n', ...
        mean(loc_zMI_ch(th_m,2),'omitnan'), ...
        mean(loc_zMI_ch(be_m,2),'omitnan'), ...
        mean(loc_zMI_ch(ga_m,2),'omitnan'));

    % ── Asignación sliced para parfor ────────────────────────────────────
    zMI_ch(j,:,:)     = loc_zMI_ch;
    zMI_nc(j,:,:)     = loc_zMI_nc;
    MI_ch(j,:,:)      = loc_MI_ch;
    MI_nc(j,:,:)      = loc_MI_nc;
    pref_ch(j,:,:)    = loc_pref_ch;
    pref_nc(j,:,:)    = loc_pref_nc;
    nv_ch(j,:)        = loc_nv_ch;
    nv_nc(j,:)        = loc_nv_nc;
    qc_phi_var(j)     = single(phi_var);
    qc_ch_label{j}    = ch_str;
    qc_fs(j)          = fs;
    qc_n_ch(j)        = numel(lat_ch);
    qc_n_nc(j)        = numel(lat_nc);

end % parfor

fprintf('\n%s\n  PARFOR FINALIZADO — %.1f min totales\n%s\n\n', ...
    repmat('-',1,82), toc/60, repmat('-',1,82));


%% ── 6. AGREGAR POR BANDA ────────────────────────────────────────────────────
band_defs = {
    'Theta',   [4  7];
    'Alpha',   [8  13];
    'Beta',    [14 30];
    'LowGamma',[30 50];
    'HiGamma', [50 80];
    'Gamma',   [30 80];
};
nBands = size(band_defs,1);

get_band = @(z3d, hz) squeeze(mean(z3d(:, P.amp_freqs>=hz(1) & P.amp_freqs<=hz(2), :), 2, 'omitnan'));

B = struct();
for bi = 1:nBands
    bn = band_defs{bi,1};
    B.Ch.(bn) = get_band(zMI_ch, band_defs{bi,2});   % [nS × nW]
    B.Nc.(bn) = get_band(zMI_nc, band_defs{bi,2});
end


%% ── 7. GUARDAR WORKSPACE ────────────────────────────────────────────────────
out_mat = fullfile(P.dir_out, 'PAC_EEG_Workspace.mat');
save(out_mat, 'zMI_ch','zMI_nc','MI_ch','MI_nc','pref_ch','pref_nc', ...
    'nv_ch','nv_nc','B','P','casos','IES_Ch','IES_Nc','Delta_IES', ...
    'qc_phi_var','qc_ch_label','qc_fs','qc_n_ch','qc_n_nc','band_defs', '-v7.3');
fprintf('>>> Workspace guardado:\n    %s\n\n', out_mat);


%% ── 8. ESTADÍSTICA Y REPORTE ────────────────────────────────────────────────
f_rep = fullfile(P.dir_reports, 'Reporte_PAC_EEG.txt');
fid   = fopen(f_rep,'w');
W     = @(s) (fprintf('%s\n',s) + fprintf(fid,'%s\n',s));

n_ok  = sum(~failed);

W(repmat('=',1,90));
W('  REPORTE PAC EEG-EEG — S5c_PAC_EEG.m');
W(sprintf('  Generado: %s', datestr(now)));
W(repmat('=',1,90));
W(sprintf('  Casos procesados  : %d / %d', n_ok, nS));
W(sprintf('  Excluidos         : %s', strjoin(P.excluded,', ')));
W(sprintf('  Driver fase       : theta [%d–%d Hz] — F1+FC1, señal continua', P.phase_band));
W(sprintf('  Amplitud          : %d–%d Hz en pasos de 4 Hz (%d frecuencias)', ...
    P.amp_freqs(1), P.amp_freqs(end), nF));
W(sprintf('  Surrogates        : permutación de trials (N=%d)', P.n_surr));
W(sprintf('  Workers           : %d (Ryzen 5900X)', P.n_workers));
W('');

% ── Tabla QC ─────────────────────────────────────────────────────────────────
W('CONTROL DE CALIDAD POR SUJETO:');
W(sprintf('  %-10s  %-8s  %-8s  %-8s  %-10s  %-18s', ...
    'Sujeto','fs','N_Ch','N_Nc','phi_var','Canales'));
W(sprintf('  %s', repmat('-',1,72)));
for j = 1:nS
    W(sprintf('  %-10s  %-8d  %-8d  %-8d  %-10.4f  %s', ...
        casos{j}, qc_fs(j), qc_n_ch(j), qc_n_nc(j), qc_phi_var(j), qc_ch_label{j}));
end
W('');

% ── Q1: PAC > 0 por banda × ventana ─────────────────────────────────────────
W('Q1 — ¿PAC EEG-EEG theta-X > surrogate? (signrank vs 0, one-tailed)');
W('');
hdr_bands = {'Theta','Alpha','Beta','Gamma'};
hdr = sprintf('  %-14s', '');
for bn = hdr_bands
    for wi = 1:nW
        hdr = [hdr sprintf(' %-11s', [bn{1}(1:min(3,end)) '-' P.wins.labels{wi}(1:3)])]; %#ok<AGROW>
    end
end
W(hdr);
for c = 1:2
    cn  = CONDS{c};
    row = sprintf('  Casos/%-8s', cn);
    for bn = hdr_bands
        for wi = 1:nW
            v = B.(cn).(bn{1})(:,wi);
            v = v(~isnan(v));
            if numel(v)>=5
                p   = signrank(v, 0, 'tail','right');
                sym = ''; if p<0.05;sym='*';end;if p<0.01;sym='**';end;if p<0.001;sym='***';end
                row = [row sprintf(' %+.2f%-3s    ', median(v), sym)]; %#ok<AGROW>
            else
                row = [row sprintf(' NaN         ')]; %#ok<AGROW>
            end
        end
    end
    W(row);
end
W('  (* p<0.05  ** p<0.01  *** p<0.001, sin corrección múltiple)');
W('');

% ── Q2: Delta Ch-Nc ──────────────────────────────────────────────────────────
W('Q2 — Delta zMI (Ch − Nc), signrank pareado');
W('');
W(hdr);
row_m = sprintf('  %-14s','Med(Ch-Nc)');
row_p = sprintf('  %-14s','p-value');
for bn = hdr_bands
    for wi = 1:nW
        ch = B.Ch.(bn{1})(:,wi); nc = B.Nc.(bn{1})(:,wi);
        v  = ~isnan(ch) & ~isnan(nc);
        if sum(v)>=5
            d   = ch(v)-nc(v);
            p   = signrank(ch(v),nc(v));
            sym = ''; if p<0.05;sym='*';end;if p<0.01;sym='**';end;if p<0.001;sym='***';end
            row_m = [row_m sprintf(' %+.3f         ', median(d))]; %#ok<AGROW>
            row_p = [row_p sprintf(' %.4f%-3s      ', p, sym)]; %#ok<AGROW>
        else
            row_m = [row_m sprintf(' NaN           ')]; %#ok<AGROW>
            row_p = [row_p sprintf(' NaN           ')]; %#ok<AGROW>
        end
    end
end
W(row_m); W(row_p);
W('');

% ── Q3: Mapa zMI media ───────────────────────────────────────────────────────
W('Q3 — Mapa zMI (media sobre sujetos válidos)');
W('');
W(hdr);
for c = 1:2
    cn  = CONDS{c};
    row = sprintf('  Casos/%-8s', cn);
    for bn = hdr_bands
        for wi = 1:nW
            v = B.(cn).(bn{1})(:,wi);
            row = [row sprintf(' %+.3f         ', mean(v,'omitnan'))]; %#ok<AGROW>
        end
    end
    W(row);
end
W('');

% ── Q4: Rayleigh por banda × ventana ─────────────────────────────────────────
W('Q4 — Test de Rayleigh (concentración de fase preferida theta)');
W('    R_Ch >> R_Nc → acoplamiento cognitivo específico de masticación');
W('');
W(sprintf('  %-18s  %5s %7s %5s %9s  %5s %7s %5s %9s', ...
    'Banda-Ventana','N_Ch','mu°_Ch','R_Ch','p_Ch','N_Nc','mu°_Nc','R_Nc','p_Nc'));
W(sprintf('  %s',repmat('-',1,82)));
for bni = 1:numel(hdr_bands)
    bn_name = hdr_bands{bni};
    % Buscar índice de esta banda en band_defs
    bd_idx = find(strcmp(band_defs(:,1), bn_name), 1);
    hz_ray = band_defs{bd_idx, 2};
    f_mask_ray = P.amp_freqs >= hz_ray(1) & P.amp_freqs <= hz_ray(2);
    for wi = 1:nW
        row = sprintf('  %-18s', [bn_name '-' P.wins.labels{wi}]);
        for c = 1:2
            if c==1
                phi_sub = squeeze(mean(pref_ch(:, f_mask_ray, wi), 2, 'omitnan'));
            else
                phi_sub = squeeze(mean(pref_nc(:, f_mask_ray, wi), 2, 'omitnan'));
            end
            phi_sub = phi_sub(~isnan(phi_sub));
            if numel(phi_sub)>=5
                Rg = abs(mean(exp(1i*phi_sub)));
                mu = angle(mean(exp(1i*phi_sub)));
                Z  = numel(phi_sub)*Rg^2;
                p  = exp(-Z)*(1+(2*Z-Z^2)/(4*numel(phi_sub)));
                p  = max(0,min(1,p));
                sym='';if p<0.001;sym='***';elseif p<0.01;sym='**';elseif p<0.05;sym='*';end
                row = [row sprintf('  %5d %7.1f %5.3f %6.4g%-3s', ...
                    numel(phi_sub),rad2deg(mu),Rg,p,sym)]; %#ok<AGROW>
            else
                row = [row sprintf('  %5s %7s %5s %9s','NA','—','—','—')]; %#ok<AGROW>
            end
        end
        W(row);
    end
    W('');
end

% ── Q5: Correlaciones PAC-EEG vs conducta ────────────────────────────────────
W('Q5 — Correlaciones zMI_Ch vs IES_Ch (Spearman, FDR-BH sobre 4 bandas × 3 ventanas = 12)');
W('');
if any(~isnan(IES_Ch))
    rho_beh_eeg = nan(nBands, nW);
    p_beh_eeg   = nan(nBands, nW);
    for bi_beh = 1:nBands
        bn = band_defs{bi_beh,1};
        for wi = 1:nW
            z = B.Ch.(bn)(:,wi);
            v = ~isnan(z) & ~isnan(IES_Ch);
            if sum(v)>=5
                [rho_beh_eeg(bi_beh,wi), p_beh_eeg(bi_beh,wi)] = ...
                    corr(z(v), IES_Ch(v), 'Type','Spearman');
            end
        end
    end

    % FDR Benjamini-Hochberg
    p_flat = p_beh_eeg(:);
    [p_s, idx_s] = sort(p_flat);
    m = numel(p_flat);
    p_fdr_eeg = nan(m,1);
    for k = 1:m
        if ~isnan(p_s(k))
            p_fdr_eeg(idx_s(k)) = min(1, p_s(k)*m/k);
        end
    end
    p_fdr_eeg = reshape(p_fdr_eeg, nBands, nW);

    W(sprintf('  %-12s  %-8s  %+7s  %8s  %8s', 'Banda','Ventana','rho','p_raw','p_FDR'));
    W(sprintf('  %s',repmat('-',1,50)));
    for bi_beh = 1:nBands
        bn = band_defs{bi_beh,1};
        for wi = 1:nW
            r = rho_beh_eeg(bi_beh,wi); pr = p_beh_eeg(bi_beh,wi); pf = p_fdr_eeg(bi_beh,wi);
            if isnan(r); continue; end
            sym='';
            if ~isnan(pf) && pf<0.05; sym='* FDR'; end
            if pr<0.05 && (isnan(pf)||pf>=0.05); sym='* raw'; end
            W(sprintf('  %-12s  %-8s  %+7.3f  %8.4f  %8.4f   %s', bn, P.wins.labels{wi}, r, pr, pf, sym));
        end
    end
    W('');

    % También Delta
    W('  Delta zMI (Ch-Nc) vs Delta IES:');
    W(sprintf('  %-12s  %-8s  %+7s  %8s','Banda','Ventana','rho','p_raw'));
    W(sprintf('  %s',repmat('-',1,38)));
    for bi_beh = 1:nBands
        bn = band_defs{bi_beh,1};
        for wi = 1:nW
            dz = B.Ch.(bn)(:,wi) - B.Nc.(bn)(:,wi);
            v  = ~isnan(dz) & ~isnan(Delta_IES);
            if sum(v)>=5
                [rd,pd] = corr(dz(v),Delta_IES(v),'Type','Spearman');
                sym='';if pd<0.05;sym='*';end;if pd<0.01;sym='**';end
                W(sprintf('  %-12s  %-8s  %+7.3f  %8.4f   %s', bn, P.wins.labels{wi}, rd, pd, sym));
            end
        end
    end
    W('');
end

W(repmat('=',1,90));
W(sprintf('  Tiempo total: %.1f min', toc/60));
fclose(fid);
fprintf('>>> Reporte guardado:\n    %s\n\n', f_rep);


%% ── 9. FIGURAS ──────────────────────────────────────────────────────────────

% ── Fig A: Comodulogram zMI medio (Ch vs Nc) — 2 paneles × 3 ventanas ───────
c_maps = {[1 0.95 0.8; 1 0.6 0; 0.7 0 0], ...  % naranja→rojo para Ch
          [0.9 0.9 1;  0.4 0.4 0.8; 0.1 0.1 0.5]};  % azul para Nc

for c = 1:2
    cn = CONDS{c};
    f  = figure('Color','w','Position',[50 50 1100 380],'Visible','off');
    tl = tiledlayout(f, 1, nW, 'TileSpacing','compact','Padding','compact');
    for wi = 1:nW
        nexttile;
        zmap = squeeze(mean(zMI_ch .* (c==1) + zMI_nc .* (c==2), 1, 'omitnan'))'; %#ok
        if c==1
            zm = squeeze(mean(zMI_ch(:,:,wi),1,'omitnan'));
        else
            zm = squeeze(mean(zMI_nc(:,:,wi),1,'omitnan'));
        end
        bar(P.amp_freqs, zm, 'FaceColor', [0.3 0.5+0.4*(c==1) 0.8*(c==2)], ...
            'EdgeColor','none','FaceAlpha',0.85);
        hold on;
        yline(0,'k--','LineWidth',1);
        yline(2,'r:','LineWidth',0.8);
        xlabel('Frecuencia amplitud (Hz)','FontSize',9);
        ylabel('zMI (media)','FontSize',9);
        title(sprintf('%s — %s',cn, P.wins.labels{wi}),'FontSize',10,'FontWeight','bold');
        xlim([P.amp_freqs(1)-2 P.amp_freqs(end)+2]);
        grid on; box off;
    end
    title(tl, sprintf('Comodulogram PAC EEG-EEG (θ-phase → amp) — Condición %s', cn), ...
        'FontSize',12,'FontWeight','bold');
    out_f = fullfile(P.dir_plots, sprintf('Fig_Comodo_%s.png', cn));
    exportgraphics(f, out_f,'Resolution',300);
    close(f);
    fprintf('>>> Fig Comodo_%s guardada\n', cn);
end

% ── Fig B: Ch - Nc difference spectrum por ventana ───────────────────────────
f = figure('Color','w','Position',[50 50 1100 380],'Visible','off');
tl = tiledlayout(f, 1, nW,'TileSpacing','compact','Padding','compact');
for wi = 1:nW
    nexttile; hold on;
    zm_ch = squeeze(mean(zMI_ch(:,:,wi),1,'omitnan'));
    zm_nc = squeeze(mean(zMI_nc(:,:,wi),1,'omitnan'));
    diff_z = zm_ch - zm_nc;
    bar(P.amp_freqs, diff_z, 'FaceColor',[0.2 0.6 0.3],'EdgeColor','none','FaceAlpha',0.85);
    yline(0,'k--'); yline(1,'r:','LineWidth',0.8);
    xlabel('Frecuencia amplitud (Hz)','FontSize',9);
    ylabel('\DeltazMI (Ch-Nc)','FontSize',9);
    title(P.wins.labels{wi},'FontSize',10,'FontWeight','bold');
    xlim([P.amp_freqs(1)-2 P.amp_freqs(end)+2]);
    grid on; box off;
end
title(tl,'PAC EEG-EEG: Ch − Nc (\DeltazMI) — Especificidad de masticación', ...
    'FontSize',12,'FontWeight','bold');
exportgraphics(f, fullfile(P.dir_plots,'Fig_Comodo_Diff_Ch_Nc.png'),'Resolution',300);
close(f);
fprintf('>>> Fig Comodo_Diff guardada\n');

% ── Fig C: Gamma-Late zMI Ch vs IES_Ch ───────────────────────────────────────
if any(~isnan(IES_Ch))
    % Encontrar la banda+ventana con menor p_raw
    [~, best_flat] = min(p_beh_eeg(:));
    [best_bi_e, best_wi_e] = ind2sub([nBands nW], best_flat);
    bn_best = band_defs{best_bi_e,1};

    z_best = B.Ch.(bn_best)(:,best_wi_e);
    v_best = ~isnan(z_best) & ~isnan(IES_Ch);
    if sum(v_best)>=5
        f = figure('Color','w','Position',[100 100 480 420],'Visible','off');
        hold on;
        scatter(z_best(v_best), IES_Ch(v_best), 65, [0.2 0.55 0.8],'filled', ...
            'MarkerFaceAlpha',0.8,'MarkerEdgeColor',[0.1 0.4 0.65]);
        lf = polyfit(z_best(v_best), IES_Ch(v_best),1);
        xl = linspace(min(z_best(v_best)), max(z_best(v_best)),60);
        plot(xl, polyval(lf,xl),'--','Color',[0.2 0.55 0.8 0.7],'LineWidth',2);
        xlabel(sprintf('zMI EEG-EEG — %s-%s (θ-phase)',bn_best,P.wins.labels{best_wi_e}),'FontSize',11);
        ylabel('IES\_m Chew (ms)','FontSize',11);
        title(sprintf('PAC EEG %s-%s vs IES\n\\rho=%+.3f  p=%.4f', ...
            bn_best, P.wins.labels{best_wi_e}, ...
            rho_beh_eeg(best_bi_e,best_wi_e), p_beh_eeg(best_bi_e,best_wi_e)), ...
            'FontSize',11,'FontWeight','bold');
        grid on; box off;
        exportgraphics(f, fullfile(P.dir_plots,'Fig_PAC_EEG_Best_vs_IES.png'),'Resolution',300);
        close(f);
        fprintf('>>> Fig PAC_EEG_Best_vs_IES guardada\n');
    end
end

fprintf('\n>>> S5c_PAC_EEG DONE — %.1f min totales.\n\n', toc/60);


%% ============================================================================
%  FUNCIONES LOCALES
%% ============================================================================

function mi = pac_mi_local(phi, amp, n_bins)
% Modulation Index (Tort et al. 2010). phi en radianes, amp envolvente positiva.
edges    = linspace(-pi, pi, n_bins+1);
mean_amp = zeros(1, n_bins);
for b = 1:n_bins
    in_b = phi >= edges(b) & phi < edges(b+1);
    if any(in_b); mean_amp(b) = mean(amp(in_b)); end
end
mean_amp(mean_amp <= 0) = eps;
p  = mean_amp / sum(mean_amp);
H  = -sum(p .* log(p + eps));
mi = (log(n_bins) - H) / log(n_bins);
end
