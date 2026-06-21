%% S6b_A3_muscle.m — A3 corregido: proxies musculares LIMPIOS
% El A3 de S6 usó 30-40 Hz en el ROI frontal (mezcla gamma neural + frontalis) → ambiguo.
% Aquí el confound muscular se prueba con proxies apropiados:
%   (a) RMS del EMG masetero real (canales 65-66, 20-200 Hz) durante chew
%   (b) Δ(30-40 Hz) en electrodos TEMPORALES (T7/T8/FT7/FT8/TP7/TP8) — donde proyecta masetero
% Hipótesis nula deseada: Δtheta(ROI) NO correlaciona con ninguno → efecto no muscular.
%
% Salida: outputs/stats/S6b_A3_muscle.mat

THIS_DIR=fileparts(mfilename('fullpath')); ROOT_V1F=fileparts(THIS_DIR);
run(fullfile(ROOT_V1F,'S0_config.m')); rng(RNG_SEED);

EEG_ref=pop_loadset(fullfile(DIR_EPOCHS,[CASES{1} '_Nc_ep.set']));
labels={EEG_ref.chanlocs(1:EEG_N).labels};
TEMPORAL={'T7','T8','FT7','FT8','TP7','TP8'};
[~,tp_idx]=ismember(TEMPORAL,labels); tp_idx(tp_idx==0)=[];

% dtheta canónico (Δtheta ROI Ch-Nc Late dB) ya calculado en S6
S6=load(fullfile(OUT_STATS,'S6_artifact_controls.mat'),'dtheta','dmusc');
dtheta=S6.dtheta; dmusc_frontal=S6.dmusc;

emg_rms=nan(N_CASES,1); dmusc_temp=nan(N_CASES,1);
cases_l=CASES; ep_l=DIR_EPOCHS; pac_l=DATA_PAC; suf_ch=EMG_SUFFIX;
emg_l=EMG_CHANS; tp_l=tp_idx; wl=WIN_LATE; wb=WIN_BASE; bmu=BAND_MUSC;

cur=path(); try, spmd, addpath(cur); end, catch; end
parfor s=1:N_CASES
    % (b) Δ(30-40) temporal desde épocas
    fch=fullfile(ep_l,[cases_l{s} '_Ch_ep.set']); fnc=fullfile(ep_l,[cases_l{s} '_Nc_ep.set']);
    if exist(fch,'file')&&exist(fnc,'file')
        Ec=pop_loadset(fch); En=pop_loadset(fnc); fse=Ec.srate; %#ok<PFBNS>
        twc=Ec.times>=wl(1)&Ec.times<=wl(2); twn=En.times>=wl(1)&En.times<=wl(2);
        tbc=Ec.times>=wb(1)&Ec.times<=wb(2); tbn=En.times>=wb(1)&En.times<=wb(2);
        mc=bandpow_db_roi(Ec,bmu(1),bmu(2),fse,EEG_N,twc,tbc,tp_l);
        mn=bandpow_db_roi(En,bmu(1),bmu(2),fse,EEG_N,twn,tbn,tp_l);
        dmusc_temp(s)=mc-mn; %#ok<PFOUS>
    end
    % (a) RMS EMG masetero (20-200 Hz) durante chew continuo
    fc=fullfile(pac_l,[cases_l{s} suf_ch]);
    if exist(fc,'file')
        E=pop_loadset(fc); fsm=E.srate; ch=emg_l(emg_l<=E.nbchan);
        if ~isempty(ch)
            emg=emg_bilateral(E.data,ch,fsm);
            [bh,ah]=butter(4,[20 min(200,fsm/2-1)]/(fsm/2),'bandpass');
            ef=filtfilt(bh,ah,emg); emg_rms(s)=sqrt(mean(ef.^2)); %#ok<PFOUS>
        end
    end
    fprintf('  %d/%d %s\n',s,N_CASES,cases_l{s});
end

%% Correlaciones (Spearman)
[r_rms,p_rms]=corr(dtheta,emg_rms,'Type','Spearman','Rows','complete');
[r_tmp,p_tmp]=corr(dtheta,dmusc_temp,'Type','Spearman','Rows','complete');
[r_fro,p_fro]=corr(dtheta,dmusc_frontal,'Type','Spearman','Rows','complete');
fprintf('\n===== A3 corregido — Δtheta(ROI) × proxies musculares =====\n');
fprintf('(a) EMG masetero RMS (real)      : rho=%+.3f p=%.4f  %s\n',r_rms,p_rms,verdict(p_rms));
fprintf('(b) Δ(30-40Hz) TEMPORAL          : rho=%+.3f p=%.4f  %s\n',r_tmp,p_tmp,verdict(p_tmp));
fprintf('(ref) Δ(30-40Hz) FRONTAL (ROI)   : rho=%+.3f p=%.4f  %s\n',r_fro,p_fro,verdict(p_fro));
fprintf('\nInterpretación: si (a) y (b) n.s. pero (ref) sig, el ρ frontal refleja\n');
fprintf('co-activación gamma neural en el mismo ROI, NO spillover del masetero.\n');

prov=struct('script',mfilename('fullpath'),'date',datestr(now,'yyyy-mm-dd HH:MM:SS'),...
    'proxies','EMG masetero RMS + Δ(30-40) temporal vs frontal','temporal',{TEMPORAL});
save(fullfile(OUT_STATS,'S6b_A3_muscle.mat'),'emg_rms','dmusc_temp','dmusc_frontal','dtheta',...
    'r_rms','p_rms','r_tmp','p_tmp','r_fro','p_fro','prov','-v7.3');
fprintf('\nGuardado: S6b_A3_muscle.mat\n');

function s=verdict(p), if p<0.05, s='SIG (preocupa)'; else, s='n.s. (ok)'; end, end
function v=bandpow_db_roi(EEG,flo,fhi,fs,n_ch,tw,tb,roi)
    pe=bandpow_db(EEG,flo,fhi,fs,n_ch,tw,tb); v=mean(pe(roi),'omitnan');
end
function pow_db=bandpow_db(EEG,flo,fhi,fs,n_ch,tw,tb)
    nt=EEG.trials; np=EEG.pnts; dc=reshape(EEG.data(1:n_ch,:,:),n_ch,nt*np);
    try, df=eegfilt(dc,fs,flo,fhi);
    catch, [bc,ac]=butter(4,[flo fhi]/(fs/2),'bandpass'); df=filtfilt(bc,ac,dc')'; end
    env=abs(hilbert(df'))'; e3=reshape(env,n_ch,np,nt); p2=e3.^2;
    pw=squeeze(mean(mean(p2(:,tw,:),2),3)); pb=squeeze(mean(mean(p2(:,tb,:),2),3));
    pow_db=(10*log10((pw+eps)./(pb+eps)))';
end
