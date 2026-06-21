% S4b_PAC_ROI.m  —  PAC con ROI theta-interacción (18 electrodos CBPT)
%
% Corrige el análisis S4 original que usó solo 4 electrodos (F1,F2,FC1,FC2).
% Usa el ROI de la interacción significativa theta (sig_int_fdr, 18 electrodos
% frontales/frontopolares identificados por el CBPT en S2c).
%
% BLOQUE A — PAC CONTINUO (señal completa)
%   Fuente: Data_PAC/{subj}_{Ch|Nc}_clean_emg.set
%   Fase:   EMG → filtro f_chew±0.5Hz → Hilbert
%   Amp:    EEG ROI-18 → filtro theta/alpha/beta → Hilbert → MI Tort
%   Null:   500 circular-shift por banda → zMI por banda
%   Grupos: Cases_Ch, Cases_Nc (comparación intra-grupo × condición)
%
% BLOQUE B — PAC POR VENTANA (sliding window, 500ms, step 100ms)
%   Ventanas: Base [-500,0]ms | Early [300,900]ms | Late [900,1300]ms
%   Fuente:   bloque continuo + markers de trials (evento 40/30)
%
% Preguntas:
%   1. ¿Existe PAC theta (vs null)? → zMI_theta > 0
%   2. ¿Es específico de banda?    → Friedman + post-hoc θ/α/β (con null)
%   3. ¿Es específico de ventana?  → Friedman Base/Early/Late × banda
%   4. ¿Difiere Ch vs Nc?          → Wilcoxon zMI_theta Ch vs Nc
%
% Salida: Analysis_V1_Final/outputs/stats/v1_S4b_PAC_ROI.mat

% S0_config — ruta directa (Analysis_v1 no tiene subcarpeta code)
s0_path = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'Analysis_v1', 'S0_config.m');
if ~exist(s0_path,'file')
    % fallback: buscar en el mismo directorio
    s0_path = fullfile(fileparts(mfilename('fullpath')), 'S0_config.m');
end
if exist(s0_path,'file')
    run(s0_path);
else
    error(['S0_config.m no encontrado en: ' s0_path ...
           '\nEjecuta >> S0_config manualmente y vuelve a correr este script.']);
end

DIR_OUT = fullfile(ROOT, 'Analysis_V1_Final', 'outputs', 'stats');
if ~exist(DIR_OUT,'dir'), mkdir(DIR_OUT); end

%% ── ROI: 18 electrodos theta-interacción FDR (canónico — igual que S4g) ─
% Definición autoritativa en S0_config.m → ROI_CBPT
% Estos 18 electrodos son los visibles en topo_theta_interaction_cases_minus_controls.png
ROI_THETA_LABELS = {'Fp1','AF7','AF3','F1','F3','F5','F7','FC1','Fpz',...
                    'Fp2','AF8','AF4','AFz','Fz','F2','F4','FC2','FCz'};
fprintf('ROI theta (%d electrodos): %s\n', numel(ROI_THETA_LABELS), strjoin(ROI_THETA_LABELS,', '));

%% ── Parámetros ───────────────────────────────────────────────────────────
EMG_BAND   = [0.8  2.2];     % Hz búsqueda f_chew (solo descriptivo)
PHASE_BAND = [0.5  2.0];     % Hz filtro EMG broadband (R1.M8: sin frecuencia fija)
BANDS_HZ   = {[4 7],[8 13],[13 30]};   % theta / alpha / beta
BAND_NAMES = {'theta','alpha','beta'};
nBands     = 3;
N_SURR     = 500;
N_BINS     = MI_BINS;        % 18 bins Tort (desde S0_config)

% Ventanas (segundos desde onset trial)
WIN_EP    = [-2.0  2.5];   % ventana de época completa
WIN_LATE  = [ 0.9  1.3];   % Late
WIN_EARLY = [ 0.3  0.9];   % Early
WIN_BASE  = [-0.5  0.0];   % Baseline

EVENT_CH  = 40;   % marker chewing
EVENT_NC  = 30;   % marker no-chew

CONDITIONS = {'Ch','Nc'};
N_COND     = 2;
SUFFIXES   = {EMG_SUFFIX, EMG_SUFFIX_NC};   % desde S0_config

%% ── f_chew individual ────────────────────────────────────────────────────
tmp = load(FILE_CHEW,'T_freq');
T_freq = tmp.T_freq;
F_chew_v2 = nan(N_CASES,1);
for s_ = 1:N_CASES
    idx_ = strcmp(T_freq.Sujeto, CASES{s_});
    if any(idx_)
        F_chew_v2(s_) = mean([T_freq.Freq_Left(idx_), T_freq.Freq_Right(idx_)],'omitnan');
    end
end
fprintf('f_chew: M=%.2f±%.2f Hz\n', mean(F_chew_v2,'omitnan'), std(F_chew_v2,'omitnan'));

%% ── Pre-allocate ─────────────────────────────────────────────────────────
% BLOQUE A — continuo (por condición × banda)
MI_cont_ch  = nan(nBands, N_CASES);   % Cases-Chewing
MI_cont_nc  = nan(nBands, N_CASES);   % Cases-NoChew
MI_null_m   = nan(N_CASES, nBands);   % media null por sujeto × banda
MI_null_s   = nan(N_CASES, nBands);   % std null
zMI_ch      = nan(N_CASES, nBands);   % z-score Ch
zMI_nc      = nan(N_CASES, nBands);   % z-score Nc

f_chew_hz   = nan(N_CASES,1);
snr_db      = nan(N_CASES,1);
pref_ch     = nan(N_CASES, nBands);
pref_nc     = nan(N_CASES, nBands);
emg_R_s     = nan(N_CASES, 1);      % per-subject mean resultant length EMG phase
emg_pref_s  = nan(N_CASES, 1);      % per-subject preferred EMG phase angle

% BLOQUE B — ventanas
MI_late_ch  = nan(nBands, N_CASES);
MI_early_ch = nan(nBands, N_CASES);
MI_base_ch  = nan(nBands, N_CASES);
MI_late_nc  = nan(nBands, N_CASES);
MI_early_nc = nan(nBands, N_CASES);
MI_base_nc  = nan(nBands, N_CASES);

n_trials_ch = nan(N_CASES,1);
n_trials_nc = nan(N_CASES,1);

%% ── Pool paralelo ────────────────────────────────────────────────────────
if isempty(gcp('nocreate'))
    try parpool('local',min(8, feature('numcores'))); catch; end
end
eeg_root = fileparts(which('eeglab'));
if ~isempty(eeg_root)
    cur_path = path();
    try
        spmd, addpath(cur_path); addpath(genpath(eeg_root)); end
    catch; end
end

%% ── Variables broadcast para parfor ─────────────────────────────────────
cases_l      = CASES;
f_chew_l     = F_chew_v2;
path_pac_l   = DATA_PAC;
suf_ch_l     = SUFFIXES{1};
suf_nc_l     = SUFFIXES{2};
roi_labs_l   = ROI_THETA_LABELS;
emg_chs_l    = EMG_USE_PAC;   % [65 66]
emg_band_l    = EMG_BAND;
phase_band_l  = PHASE_BAND;
bands_l       = BANDS_HZ;
n_surr_l     = N_SURR;
n_bins_l     = N_BINS;
win_ep_l     = WIN_EP;
win_late_l   = WIN_LATE;
win_early_l  = WIN_EARLY;
win_base_l   = WIN_BASE;
ev_ch_l      = EVENT_CH;
ev_nc_l      = EVENT_NC;

%% ── Loop principal ───────────────────────────────────────────────────────
fprintf('Procesando %d sujetos × 2 condiciones con ROI=%d electrodos...\n',...
        N_CASES, numel(ROI_THETA_LABELS));

parfor s = 1:N_CASES
    fname_ch = fullfile(path_pac_l, [cases_l{s} suf_ch_l]);
    fname_nc = fullfile(path_pac_l, [cases_l{s} suf_nc_l]);

    if ~exist(fname_ch,'file'), fprintf('  Ch no encontrado: %s\n',cases_l{s}); continue; end

    % ── Cargar archivo Chewing (continuo) ─────────────────────────────────
    EEG = pop_loadset(fname_ch); %#ok<PFBNS>
    fs  = EEG.srate;
    T   = EEG.pnts;

    % ── EMG: mejor canal por SNR ──────────────────────────────────────────
    ch_emg = emg_chs_l(emg_chs_l <= EEG.nbchan);
    if isempty(ch_emg), continue; end
    [b_hf,a_hf] = butter(4,[20 min(400,fs/2-1)]/(fs/2),'bandpass');
    snr_v = nan(1,numel(ch_emg));
    for ci = 1:numel(ch_emg)
        sig_ci  = detrend(double(EEG.data(ch_emg(ci),:)));
        env_ci  = detrend(abs(filtfilt(b_hf,a_hf,sig_ci)));
        [px,fx] = pwelch(env_ci,round(fs*2),round(fs),[],fs);
        in_b  = fx>=1.0 & fx<=2.5;
        out_b = (fx>=0.1 & fx<1.0) | (fx>2.5 & fx<=5);
        if any(in_b) && any(out_b)
            snr_v(ci) = 10*log10(mean(px(in_b))/mean(px(out_b)));
        end
    end
    [best_snr, best_ci] = max(snr_v);
    snr_db(s) = best_snr; %#ok<PFOUS>
    emg_raw   = detrend(double(EEG.data(ch_emg(best_ci),:)));

    % f_chew individual
    if ~isnan(f_chew_l(s))
        f_c = f_chew_l(s);
    else
        [b_m,a_m] = butter(4,[1.0 2.5]/(fs/2),'bandpass');
        env_m = detrend(abs(filtfilt(b_m,a_m,emg_raw)));
        [pxx,fx] = pwelch(env_m,round(fs*8),round(fs*4),[],fs);
        fi = fx>=emg_band_l(1) & fx<=emg_band_l(2);
        if any(fi), [~,pk]=max(pxx(fi)); fxs=fx(fi); f_c=fxs(pk);
        else, f_c=1.5; end
    end
    f_chew_hz(s) = f_c; %#ok<PFOUS>

    % ── Fase EMG: f_chew_individual ± 0.5 Hz (R1.M8) ────────────────────────
    lo_ph = max(0.3, f_c - 0.5);
    hi_ph = min(fs/2 - 0.5, f_c + 0.5);
    [b_ph,a_ph] = butter(4, [lo_ph hi_ph]/(fs/2), 'bandpass');
    phase_ch = angle(hilbert(filtfilt(b_ph,a_ph,emg_raw)));  % [1 x T]

    % Uniformidad fase EMG por sujeto (debe ser ~0 para oscilador uniforme)
    ep_s = mean(exp(1i*phase_ch));
    emg_R_s(s)    = abs(ep_s);     %#ok<PFOUS>
    emg_pref_s(s) = angle(ep_s);   %#ok<PFOUS>

    % ── ROI EEG: average de 18 electrodos ────────────────────────────────
    n_eeg = min(64, EEG.nbchan);
    [~, roi_idx] = ismember(roi_labs_l, {EEG.chanlocs(1:n_eeg).labels});
    roi_idx(roi_idx==0) = [];
    if isempty(roi_idx), continue; end
    eeg_roi_ch = mean(double(EEG.data(roi_idx,:)), 1);   % [1 x T]

    % ── BLOQUE A: MI continuo + null (3 bandas) ───────────────────────────
    be  = linspace(-pi,pi,n_bins_l+1);
    bc  = (be(1:end-1)+be(2:end))/2;
    loc_mi_ch   = nan(nBands,1);
    loc_null_m  = nan(nBands,1);
    loc_null_s  = nan(nBands,1);
    loc_zmi_ch  = nan(nBands,1);
    loc_pref_ch = nan(nBands,1);
    min_sh = round(1.0*fs);   % ≥1 s shift → >0.5 ciclos a 0.5 Hz (banda más lenta)

    for b = 1:nBands
        [bb,ab] = butter(4, bands_l{b}/(fs/2), 'bandpass');
        amp_b   = abs(hilbert(filtfilt(bb,ab,eeg_roi_ch)));
        mi_b    = tort_mi(phase_ch', amp_b', n_bins_l);
        loc_mi_ch(b) = mi_b;

        % Fase preferida
        ad = zeros(n_bins_l,1);
        for k = 1:n_bins_l
            ik = phase_ch >= be(k) & phase_ch < be(k+1);
            if any(ik), ad(k) = mean(amp_b(ik)); end
        end
        ad = ad/(sum(ad)+eps);
        loc_pref_ch(b) = angle(sum(ad .* exp(1i*bc')));

        % Null (surrogados)
        mi_surr = nan(n_surr_l,1);
        for k = 1:n_surr_l
            sh = randi([min_sh, T-min_sh]);
            mi_surr(k) = tort_mi(phase_ch', circshift(amp_b,sh)', n_bins_l);
        end
        loc_null_m(b) = mean(mi_surr,'omitnan');
        loc_null_s(b) = std(mi_surr,'omitnan');
        loc_zmi_ch(b) = (mi_b - loc_null_m(b)) / (loc_null_s(b)+eps);
    end

    MI_cont_ch(:,s)  = loc_mi_ch;  %#ok<PFOUS>
    MI_null_m(s,:)   = loc_null_m'; %#ok<PFOUS>
    MI_null_s(s,:)   = loc_null_s'; %#ok<PFOUS>
    zMI_ch(s,:)      = loc_zmi_ch'; %#ok<PFOUS>
    pref_ch(s,:)     = loc_pref_ch'; %#ok<PFOUS>

    % ── BLOQUE B: Ventanas (sliding window desde señal Ch) ────────────────
    wS   = round(win_ep_l * fs);
    Lwin = diff(wS)+1;

    % Encontrar markers evento 40
    ev_b = [];
    for e = 1:numel(EEG.event)
        et  = EEG.event(e).type;
        if isequal(et,ev_ch_l) || strcmp(et,num2str(ev_ch_l))
            lat = round(EEG.event(e).latency);
            ix  = lat+wS(1):lat+wS(2);
            if ix(1)>=1 && ix(end)<=T
                ev_b(end+1) = lat; %#ok<AGROW>
            end
        end
    end
    nev = numel(ev_b);
    n_trials_ch(s) = nev; %#ok<PFOUS>

    loc_mi_late  = nan(nBands,1);
    loc_mi_early = nan(nBands,1);
    loc_mi_base  = nan(nBands,1);

    if nev >= 5
        % Señal inducida (ERP subtract)
        erp = zeros(1,Lwin);
        for k = 1:nev
            erp = erp + eeg_roi_ch(ev_b(k)+wS(1):ev_b(k)+wS(2));
        end
        erp = erp/nev;
        eeg_ind = eeg_roi_ch;
        for k = 1:nev
            ix_k = ev_b(k)+wS(1):ev_b(k)+wS(2);
            eeg_ind(ix_k) = eeg_roi_ch(ix_k) - erp;
        end

        % Amplitud por banda (señal inducida)
        amp_bands = nan(nBands, T);
        for b = 1:nBands
            [bb,ab] = butter(4, bands_l{b}/(fs/2), 'bandpass');
            amp_bands(b,:) = abs(hilbert(filtfilt(bb,ab,eeg_ind)));
        end

        % Epocar fase + amplitud
        phase_ep = nan(Lwin, nev);
        amp_ep   = nan(nBands, Lwin, nev);
        for k = 1:nev
            ix_k = ev_b(k)+wS(1):ev_b(k)+wS(2);
            phase_ep(:,k) = phase_ch(ix_k)';
            for b = 1:nBands
                amp_ep(b,:,k) = amp_bands(b,ix_k);
            end
        end

        % Sliding window MI
        slide_s = round(0.5*fs); step_s = round(0.1*fs);
        t_cen   = 1:step_s:Lwin;
        hw      = floor(slide_s/2);
        MI_tc   = nan(nBands, numel(t_cen));
        for ti = 1:numel(t_cen)
            t_s = max(1, t_cen(ti)-hw);  t_e = min(Lwin, t_cen(ti)+hw);
            idx_w = t_s:t_e;
            if numel(idx_w) < 10, continue; end
            phi_p = reshape(phase_ep(idx_w,:),[],1);
            for b = 1:nBands
                amp_p = reshape(amp_ep(b,idx_w,:),[],1);
                MI_tc(b,ti) = tort_mi(phi_p, amp_p, n_bins_l);
            end
        end

        % Extraer ventanas
        tc_ms = (t_cen-1)/fs*1000 + win_ep_l(1)*1000;
        for b = 1:nBands
            il  = tc_ms >= win_late_l(1)*1000  & tc_ms <= win_late_l(2)*1000;
            ie  = tc_ms >= win_early_l(1)*1000 & tc_ms <= win_early_l(2)*1000;
            ib_ = tc_ms >= win_base_l(1)*1000   & tc_ms <= win_base_l(2)*1000;
            if any(il),  loc_mi_late(b)  = mean(MI_tc(b,il), 'omitnan');  end
            if any(ie),  loc_mi_early(b) = mean(MI_tc(b,ie), 'omitnan');  end
            if any(ib_), loc_mi_base(b)  = mean(MI_tc(b,ib_),'omitnan'); end
        end
    end

    MI_late_ch(:,s)  = loc_mi_late;  %#ok<PFOUS>
    MI_early_ch(:,s) = loc_mi_early; %#ok<PFOUS>
    MI_base_ch(:,s)  = loc_mi_base;  %#ok<PFOUS>

    % ── NoChew: solo Bloque A (MI continuo + zMI) ────────────────────────
    if exist(fname_nc,'file')
        EEG_nc = pop_loadset(fname_nc);
        fs_nc  = EEG_nc.srate;
        T_nc   = EEG_nc.pnts;
        ch_emg_nc = emg_chs_l(emg_chs_l <= EEG_nc.nbchan);
        if ~isempty(ch_emg_nc)
            % EMG NoChew — mismo f_chew individual ± 0.5 Hz, promedio de canales
            emg_nc = detrend(mean(double(EEG_nc.data(ch_emg_nc,:)), 1));
            [b_pn,a_pn] = butter(4, [lo_ph hi_ph]/(fs_nc/2), 'bandpass');
            phase_nc = angle(hilbert(filtfilt(b_pn,a_pn,emg_nc)));

            n_eeg_nc = min(64, EEG_nc.nbchan);
            [~, roi_nc] = ismember(roi_labs_l, {EEG_nc.chanlocs(1:n_eeg_nc).labels});
            roi_nc(roi_nc==0) = [];
            if ~isempty(roi_nc)
                eeg_roi_nc = mean(double(EEG_nc.data(roi_nc,:)), 1);
                loc_mi_nc  = nan(nBands,1);
                loc_zmi_nc = nan(nBands,1);
                loc_pref_nc = nan(nBands,1);
                min_sh_nc  = round(1.0*fs_nc);
                for b = 1:nBands
                    [bb,ab] = butter(4, bands_l{b}/(fs_nc/2), 'bandpass');
                    amp_nc  = abs(hilbert(filtfilt(bb,ab,eeg_roi_nc)));
                    mi_nc_b = tort_mi(phase_nc', amp_nc', n_bins_l);
                    loc_mi_nc(b) = mi_nc_b;

                    % Null (usa null calculado en Ch para comparabilidad)
                    null_m_b = loc_null_m(b);
                    null_s_b = loc_null_s(b);
                    loc_zmi_nc(b) = (mi_nc_b - null_m_b)/(null_s_b+eps);

                    % Fase preferida Nc
                    ad_nc = zeros(n_bins_l,1);
                    for k = 1:n_bins_l
                        ik = phase_nc >= be(k) & phase_nc < be(k+1);
                        if any(ik), ad_nc(k)=mean(amp_nc(ik)); end
                    end
                    ad_nc = ad_nc/(sum(ad_nc)+eps);
                    loc_pref_nc(b) = angle(sum(ad_nc.*exp(1i*bc')));
                end
                MI_cont_nc(:,s) = loc_mi_nc;  %#ok<PFOUS>
                zMI_nc(s,:)     = loc_zmi_nc'; %#ok<PFOUS>
                pref_nc(s,:)    = loc_pref_nc'; %#ok<PFOUS>

                % Bloque B NoChew (ventanas)
                ev_nc_b = [];
                for e = 1:numel(EEG_nc.event)
                    et_n = EEG_nc.event(e).type;
                    if isequal(et_n,ev_nc_l)||strcmp(et_n,num2str(ev_nc_l))
                        lat_n = round(EEG_nc.event(e).latency);
                        ix_n  = lat_n+wS(1):lat_n+wS(2);
                        if ix_n(1)>=1 && ix_n(end)<=T_nc
                            ev_nc_b(end+1) = lat_n; %#ok<AGROW>
                        end
                    end
                end
                n_trials_nc(s) = numel(ev_nc_b); %#ok<PFOUS>

                if numel(ev_nc_b) >= 5
                    erp_nc = zeros(1,Lwin);
                    for k = 1:numel(ev_nc_b)
                        erp_nc=erp_nc+eeg_roi_nc(ev_nc_b(k)+wS(1):ev_nc_b(k)+wS(2));
                    end
                    erp_nc = erp_nc/numel(ev_nc_b);
                    eeg_ind_nc = eeg_roi_nc;
                    for k = 1:numel(ev_nc_b)
                        ix_kn = ev_nc_b(k)+wS(1):ev_nc_b(k)+wS(2);
                        eeg_ind_nc(ix_kn) = eeg_roi_nc(ix_kn)-erp_nc;
                    end
                    amp_bnds_nc = nan(nBands, T_nc);
                    for b = 1:nBands
                        [bb,ab]=butter(4,bands_l{b}/(fs_nc/2),'bandpass');
                        amp_bnds_nc(b,:)=abs(hilbert(filtfilt(bb,ab,eeg_ind_nc)));
                    end
                    ph_nc_ep = nan(Lwin,numel(ev_nc_b));
                    amp_nc_ep= nan(nBands,Lwin,numel(ev_nc_b));
                    for k=1:numel(ev_nc_b)
                        ix_kn=ev_nc_b(k)+wS(1):ev_nc_b(k)+wS(2);
                        ph_nc_ep(:,k)=phase_nc(ix_kn)';
                        for b=1:nBands, amp_nc_ep(b,:,k)=amp_bnds_nc(b,ix_kn); end
                    end
                    MI_tc_nc = nan(nBands,numel(t_cen));
                    for ti=1:numel(t_cen)
                        t_s=max(1,t_cen(ti)-hw); t_e=min(Lwin,t_cen(ti)+hw);
                        idx_w=t_s:t_e; if numel(idx_w)<10, continue; end
                        phi_pn=reshape(ph_nc_ep(idx_w,:),[],1);
                        for b=1:nBands
                            amp_pn=reshape(amp_nc_ep(b,idx_w,:),[],1);
                            MI_tc_nc(b,ti)=tort_mi(phi_pn,amp_pn,n_bins_l);
                        end
                    end
                    for b=1:nBands
                        il_n =tc_ms>=win_late_l(1)*1000 &tc_ms<=win_late_l(2)*1000;
                        ie_n =tc_ms>=win_early_l(1)*1000&tc_ms<=win_early_l(2)*1000;
                        ib_n =tc_ms>=win_base_l(1)*1000 &tc_ms<=win_base_l(2)*1000;
                        if any(il_n),  MI_late_nc(b,s) =mean(MI_tc_nc(b,il_n), 'omitnan');  end %#ok<PFOUS>
                        if any(ie_n),  MI_early_nc(b,s)=mean(MI_tc_nc(b,ie_n), 'omitnan');  end %#ok<PFOUS>
                        if any(ib_n),  MI_base_nc(b,s) =mean(MI_tc_nc(b,ib_n), 'omitnan');  end %#ok<PFOUS>
                    end
                end
            end
        end
    end

    if mod(s,5)==0, fprintf('  %d/%d\n',s,N_CASES); end
end

fprintf('\n✓ Procesamiento completo.\n');

%% ── ESTADÍSTICOS ─────────────────────────────────────────────────────────
fprintf('\n======================================================\n');
fprintf('PREGUNTA 1: ¿Existe PAC theta en Cases-Chewing?\n');
fprintf('======================================================\n');
k_sig_ch = sum(zMI_ch(:,1)>1.96,'omitnan');
nVal     = sum(~isnan(zMI_ch(:,1)));
p_binom  = 1 - binocdf(k_sig_ch-1, nVal, 0.05);
[~,p_w1] = signrank(zMI_ch(:,1));
fprintf('zMI theta Ch: M=%.2f  Mdn=%.2f  | %d/%d sig | binom p=%.2e | Wilcoxon p=%.2e\n',...
        mean(zMI_ch(:,1),'omitnan'), median(zMI_ch(:,1),'omitnan'),...
        k_sig_ch, nVal, p_binom, p_w1);

fprintf('\n======================================================\n');
fprintf('PREGUNTA 2: ¿Es específico de banda? (zMI por banda)\n');
fprintf('======================================================\n');
for b = 1:nBands
    k_b = sum(zMI_ch(:,b)>1.96,'omitnan');
    [~,pw] = signrank(zMI_ch(:,b));
    fprintf('%s: M=%.2f  Mdn=%.2f  %d/%d sig  p=%.4f\n',...
            BAND_NAMES{b}, mean(zMI_ch(:,b),'omitnan'),...
            median(zMI_ch(:,b),'omitnan'), k_b, nVal, pw);
end
% Friedman sobre zMI × banda
[p_fr_z,~] = friedman(zMI_ch, 1, 'off');
fprintf('Friedman zMI Ch × banda: p=%.4f\n', p_fr_z);

fprintf('\n======================================================\n');
fprintf('PREGUNTA 3: ¿Es específico de ventana? (theta)\n');
fprintf('======================================================\n');
for b = 1:nBands
    [~,pb] = signrank(MI_late_ch(b,:)', MI_base_ch(b,:)');
    [~,pe] = signrank(MI_early_ch(b,:)', MI_base_ch(b,:)');
    [~,pl] = signrank(MI_late_ch(b,:)', MI_early_ch(b,:)');
    fprintf('%s — Base→Early p=%.3f | Base→Late p=%.3f | Early→Late p=%.3f\n',...
            BAND_NAMES{b}, pe, pb, pl);
end
data_win = [MI_base_ch(1,:)' MI_early_ch(1,:)' MI_late_ch(1,:)'];
[p_fw,~] = friedman(data_win,1,'off');
fprintf('Friedman theta × ventana (Base/Early/Late): p=%.4f\n', p_fw);

fprintf('\n======================================================\n');
fprintf('PREGUNTA 4: ¿Difiere Ch vs Nc (theta)?\n');
fprintf('======================================================\n');
[~,p_cond] = signrank(zMI_ch(:,1), zMI_nc(:,1));
fprintf('zMI theta Ch: M=%.2f  vs  Nc: M=%.2f  |  paired Wilcoxon p=%.4f\n',...
        mean(zMI_ch(:,1),'omitnan'), mean(zMI_nc(:,1),'omitnan'), p_cond);
[~,p_mi_cond] = signrank(MI_cont_ch(1,:)', MI_cont_nc(1,:)');
fprintf('MI theta Ch vs Nc (abs): p=%.4f\n', p_mi_cond);

%% ── RAYLEIGH: fase preferida por banda (responde R1.M9) ─────────────────
fprintf('\n======================================================\n');
fprintf('RAYLEIGH: Fase preferida (grupo, n=%d)\n', N_CASES);
fprintf('======================================================\n');

rayl_cont_R = nan(nBands,1);
rayl_cont_Z = nan(nBands,1);
rayl_cont_p = nan(nBands,1);

for b = 1:nBands
    pv_c = pref_ch(:,b);  pv_c = pv_c(~isnan(pv_c));
    nv   = numel(pv_c);
    Rc   = abs(mean(exp(1i*pv_c)));
    Zc   = nv * Rc^2;
    pc   = exp(-Zc);   % Rayleigh p-value, exacto para n>=10
    rayl_cont_R(b) = Rc;
    rayl_cont_Z(b) = Zc;
    rayl_cont_p(b) = pc;
    fprintf('Continuo %s: R=%.3f  Z=%.2f  p=%.4f  (n=%d)\n',...
            BAND_NAMES{b}, Rc, Zc, pc, nv);
end

% Rayleigh EMG uniformidad (debe ser NS → confirma que la fase no está sesgada)
fprintf('\n-- Uniformidad fase EMG (R≈0 esperado para oscilador uniforme) --\n');
fprintf('R individual: M=%.4f  SD=%.4f  [rango: %.4f – %.4f]\n',...
        mean(emg_R_s,'omitnan'), std(emg_R_s,'omitnan'),...
        min(emg_R_s,[],'omitnan'), max(emg_R_s,[],'omitnan'));
% Rayleigh sobre ángulo preferido EMG entre sujetos
pref_emg_v = emg_pref_s(~isnan(emg_pref_s));
ne = numel(pref_emg_v);
Re = abs(mean(exp(1i*pref_emg_v)));
Ze = ne * Re^2;
pe = max(0, exp(-Ze) * (1 + (2*Ze - Ze^2)/(4*ne)));
fprintf('Rayleigh EMG preferred angle (N=%d): R=%.3f  Z=%.2f  p=%.4f\n', ne, Re, Ze, pe);
fprintf('(NS esperado: los sujetos no comparten un ángulo EMG preferido en común)\n');

%% ── CORRELACIÓN PAC × RT (R1.M2) ────────────────────────────────────────
fprintf('\n======================================================\n');
fprintf('CORRELACIÓN PAC × Conducta (R1.M2)\n');
fprintf('======================================================\n');

FILE_BEHAV = fullfile(ROOT, 'Analysis_V1_Final', 'data', 'computed', 'v1_S1_behavior_stats.mat');
Beh = load(FILE_BEHAV, 'rt_cas_b2', 'ies_cas_b2');
RT_ch  = Beh.rt_cas_b2(:);    % RT mediana bloque 2 (Chewing), n=31
IES_ch = Beh.ies_cas_b2(:);   % IES bloque 2 (Chewing), n=31

[rho_zrt,  p_zrt]  = corr(zMI_ch(:,1),     RT_ch,  'Type','Spearman','Rows','complete');
[rho_mrt,  p_mrt]  = corr(MI_late_ch(1,:)', RT_ch,  'Type','Spearman','Rows','complete');
[rho_zies, p_zies] = corr(zMI_ch(:,1),     IES_ch, 'Type','Spearman','Rows','complete');
[rho_mies, p_mies] = corr(MI_late_ch(1,:)', IES_ch, 'Type','Spearman','Rows','complete');

fprintf('rho(zMI_theta_cont × RT_b2)  = %+.3f  p = %.4f\n', rho_zrt,  p_zrt);
fprintf('rho(MI_theta_late  × RT_b2)  = %+.3f  p = %.4f\n', rho_mrt,  p_mrt);
fprintf('rho(zMI_theta_cont × IES_b2) = %+.3f  p = %.4f\n', rho_zies, p_zies);
fprintf('rho(MI_theta_late  × IES_b2) = %+.3f  p = %.4f\n', rho_mies, p_mies);

%% ── Guardar ──────────────────────────────────────────────────────────────
ROI_PAC = ROI_THETA_LABELS;
save(fullfile(DIR_OUT,'v1_S4b_PAC_ROI.mat'),...
     'MI_cont_ch','MI_cont_nc','MI_null_m','MI_null_s',...
     'zMI_ch','zMI_nc','pref_ch','pref_nc',...
     'MI_late_ch','MI_early_ch','MI_base_ch',...
     'MI_late_nc','MI_early_nc','MI_base_nc',...
     'f_chew_hz','snr_db','n_trials_ch','n_trials_nc',...
     'emg_R_s','emg_pref_s',...
     'rayl_cont_R','rayl_cont_Z','rayl_cont_p',...
     'ROI_PAC','BAND_NAMES','BANDS_HZ','PHASE_BAND',...
     'WIN_LATE','WIN_EARLY','WIN_BASE','N_SURR','-v7.3');
fprintf('\nGuardado: v1_S4b_PAC_ROI.mat\n');

% tort_mi está en tort_mi.m (mismo directorio) — no función local para
% evitar que MATLAB ejecute este script en workspace propio.
