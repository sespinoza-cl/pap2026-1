function res = pac_subject(s, P)
% pac_subject  Todo el cómputo PAC de UN sujeto (Ch + Nc), reproducible.
%   s : índice de sujeto
%   P : struct con parámetros broadcast (ver S4b_PAC_ROI.m)
%   res : struct con todos los resultados del sujeto (NaN si falta dato)

    nB = numel(P.bands);
    res = init_res(nB);

    stream = RandStream('Threefry','Seed',P.seed); stream.Substream = s;

    fch = fullfile(P.pac, [P.cases{s} P.suf_ch]);
    fnc = fullfile(P.pac, [P.cases{s} P.suf_nc]);
    if ~exist(fch,'file'), fprintf('  [SKIP] %s sin Ch\n',P.cases{s}); return; end

    % ===== CHEW =====
    EEG = pop_loadset(fch); fs = EEG.srate; T = EEG.pnts;
    ch_emg = P.emg_chs(P.emg_chs <= EEG.nbchan); if isempty(ch_emg), return; end
    [emg, used_e, snr_e] = emg_bilateral(EEG.data, ch_emg, fs);
    res.emg_used_ch = used_e; res.emg_snr_ch = snr_e;

    if ~isnan(P.fchew(s)), f_c = P.fchew(s);
    else
        [bm,am]=butter(4,[1.0 2.5]/(fs/2),'bandpass');
        env=detrend(abs(filtfilt(bm,am,emg)));
        [pxx,fx]=pwelch(env,round(fs*8),round(fs*4),[],fs);
        fi=fx>=0.8&fx<=2.2; if any(fi),[~,pk]=max(pxx(fi));fxs=fx(fi);f_c=fxs(pk);else,f_c=1.5;end
    end
    res.f_chew = f_c;

    lo=max(0.3,f_c-P.ph_half); hi=min(fs/2-0.5,f_c+P.ph_half);
    [bp,ap]=butter(4,[lo hi]/(fs/2),'bandpass');
    phase_ch = angle(hilbert(filtfilt(bp,ap,emg)));

    n_eeg=min(64,EEG.nbchan);
    [~,ri]=ismember(P.roi,{EEG.chanlocs(1:n_eeg).labels}); ri(ri==0)=[];
    if isempty(ri), return; end
    roi_ch = mean(double(EEG.data(ri,:)),1);

    Rc = compute_pac_cont(phase_ch, roi_ch, P.bands, fs, P.nb, P.ns, P.minsh, stream);
    res.MI_ch=Rc.mi; res.zC_ch=Rc.zc'; res.zA_ch=Rc.za'; res.pref_ch=Rc.pref';
    res.nullC_m=Rc.ncm'; res.nullC_s=Rc.ncs'; res.nullA_m=Rc.nam'; res.nullA_s=Rc.nas';

    ev_ch = grab_events(EEG, P.ev_ch, P.wS, T);
    res.n_trials_ch = numel(ev_ch);
    res.MIw_ch = compute_pac_windows(roi_ch, phase_ch, ev_ch, P.wS, P.bands, fs, P.nb, P.wins);

    % ===== NOCHEW (mismo f_chew, su PROPIO null) =====
    if exist(fnc,'file')
        EEGn=pop_loadset(fnc); fsn=EEGn.srate; Tn=EEGn.pnts;
        ch_emgn=P.emg_chs(P.emg_chs<=EEGn.nbchan);
        if ~isempty(ch_emgn)
            [emgn,used_n,snr_n]=emg_bilateral(EEGn.data,ch_emgn,fsn);
            res.emg_used_nc=used_n; res.emg_snr_nc=snr_n;
            lon=max(0.3,f_c-P.ph_half); hin=min(fsn/2-0.5,f_c+P.ph_half);
            [bpn,apn]=butter(4,[lon hin]/(fsn/2),'bandpass');
            phase_nc=angle(hilbert(filtfilt(bpn,apn,emgn)));
            ne=min(64,EEGn.nbchan);
            [~,rin]=ismember(P.roi,{EEGn.chanlocs(1:ne).labels}); rin(rin==0)=[];
            if ~isempty(rin)
                roi_nc=mean(double(EEGn.data(rin,:)),1);
                Rn=compute_pac_cont(phase_nc, roi_nc, P.bands, fsn, P.nb, P.ns, P.minsh, stream);
                res.MI_nc=Rn.mi; res.zC_nc=Rn.zc'; res.zA_nc=Rn.za'; res.pref_nc=Rn.pref';
                ev_nc=grab_events(EEGn, P.ev_nc, P.wS, Tn);
                res.n_trials_nc=numel(ev_nc);
                res.MIw_nc=compute_pac_windows(roi_nc, phase_nc, ev_nc, P.wS, P.bands, fsn, P.nb, P.wins);
            end
        end
    end
    fprintf('  [%d] %s OK (f_chew=%.2f Hz, EMG used=%d)\n', s, P.cases{s}, f_c, used_e);
end

function res = init_res(nB)
    z=@()nan(1,nB);
    res=struct('MI_ch',z(),'zC_ch',z(),'zA_ch',z(),'pref_ch',z(),...
        'nullC_m',z(),'nullC_s',z(),'nullA_m',z(),'nullA_s',z(),...
        'MI_nc',z(),'zC_nc',z(),'zA_nc',z(),'pref_nc',z(),...
        'MIw_ch',nan(nB,3),'MIw_nc',nan(nB,3),'f_chew',NaN,...
        'emg_used_ch',NaN,'emg_used_nc',NaN,'emg_snr_ch',NaN,'emg_snr_nc',NaN,...
        'n_trials_ch',NaN,'n_trials_nc',NaN);
    res.MI_ch=res.MI_ch(:); res.MI_nc=res.MI_nc(:);  % columnas
end

function ev = grab_events(EEG, code, wS, T)
    ev = [];
    for e=1:numel(EEG.event)
        et=EEG.event(e).type;
        if isequal(et,code)||strcmp(num2str(et),num2str(code))
            lat=round(EEG.event(e).latency);
            if lat+wS(1)>=1 && lat+wS(2)<=T, ev(end+1)=lat; end %#ok<AGROW>
        end
    end
end
