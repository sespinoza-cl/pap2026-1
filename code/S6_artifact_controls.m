%% S6_artifact_controls.m — Controles de artefacto motor A1-A6 (converge el blindaje)
% A1 disociación espacial  : theta frontal-medial sí / temporal no (theta_topo, t+FDR)
% A3 proxy muscular         : Δtheta(ROI) × Δ(30-40Hz)(ROI), Spearman → debe ser n.s.
% A4(ii) sharpness          : razón armónica envolvente EMG (2f/f) × zMI theta → n.s. debilita artefacto
% (A2 FOOOF, A4(i) AAFT, A5 Ch>Nc, A6 conducta ya en sus .mat — se consolidan en el .md)
%
% Salida: outputs/stats/S6_artifact_controls.mat (+ provenance)

THIS_DIR=fileparts(mfilename('fullpath')); ROOT_V1F=fileparts(THIS_DIR);
run(fullfile(ROOT_V1F,'S0_config.m')); rng(RNG_SEED);

%% chanlocs + índices ROI / temporal
EEG_ref=pop_loadset(fullfile(DIR_EPOCHS,[CASES{1} '_Nc_ep.set']));
labels={EEG_ref.chanlocs(1:EEG_N).labels};
[~,roi_idx]=ismember(ROI_CBPT,labels); roi_idx(roi_idx==0)=[];
FRONTAL_MED={'Fpz','AFz','Fz','FCz','AF3','AF4','F1','F2'};
TEMPORAL   ={'T7','T8','FT7','FT8','TP7','TP8'};
[~,fm_idx]=ismember(FRONTAL_MED,labels); fm_idx(fm_idx==0)=[];
[~,tp_idx]=ismember(TEMPORAL,labels);    tp_idx(tp_idx==0)=[];

%% ===== A1: disociación espacial (theta_topo Morlet) =====
S=load(FILE_TF,'theta_topo_ch','theta_topo_nc');
dth=S.theta_topo_ch - S.theta_topo_nc;           % [31×64] Δtheta Ch-Nc
[~,p_el,~,st]=ttest(dth); t_el=st.tstat;
p_fdr=fdr_bh(p_el(:));
fprintf('\n===== A1: DISOCIACIÓN ESPACIAL (Δtheta Ch-Nc por electrodo) =====\n');
fprintf('Frontal-medial:\n');
for i=fm_idx, fprintf('  %-4s t=%+.2f p=%.4f pFDR=%.4f %s\n',labels{i},t_el(i),p_el(i),p_fdr(i),tern(p_fdr(i)<0.05,'*','')); end
fprintf('Temporal (proxy muscular):\n');
for i=tp_idx, fprintf('  %-4s t=%+.2f p=%.4f pFDR=%.4f %s\n',labels{i},t_el(i),p_el(i),p_fdr(i),tern(p_fdr(i)<0.05,'*','')); end
n_fm_sig=sum(p_fdr(fm_idx)<0.05); n_tp_sig=sum(p_fdr(tp_idx)<0.05);
fprintf('Veredicto A1: frontal-medial sig=%d/%d | temporal sig=%d/%d\n',...
    n_fm_sig,numel(fm_idx),n_tp_sig,numel(tp_idx));

%% ===== A3 + A4(ii): loop por caso =====
dtheta=nan(N_CASES,1); dmusc=nan(N_CASES,1); HR=nan(N_CASES,1);
cases_l=CASES; ep_l=DIR_EPOCHS; pac_l=DATA_PAC; suf_ch=EMG_SUFFIX; suf_nc=EMG_SUFFIX_NC;
emg_l=EMG_CHANS; roi_l=roi_idx; wl=WIN_LATE; wb=WIN_BASE; bth=BAND_THETA; bmu=BAND_MUSC; fs_cfg=FS;

cur=path(); try, spmd, addpath(cur); end, catch; end
parfor s=1:N_CASES
    % A3: ROI Δtheta y Δmusc desde épocas (Late dB, Ch-Nc)
    fch=fullfile(ep_l,[cases_l{s} '_Ch_ep.set']); fnc=fullfile(ep_l,[cases_l{s} '_Nc_ep.set']);
    if exist(fch,'file')&&exist(fnc,'file')
        Ec=pop_loadset(fch); En=pop_loadset(fnc); fse=Ec.srate; %#ok<PFBNS>
        twc=Ec.times>=wl(1)&Ec.times<=wl(2); twn=En.times>=wl(1)&En.times<=wl(2);
        tbc=Ec.times>=wb(1)&Ec.times<=wb(2); tbn=En.times>=wb(1)&En.times<=wb(2);
        th_c=bandpow_db_roi(Ec,bth(1),bth(2),fse,EEG_N,twc,tbc,roi_l);
        th_n=bandpow_db_roi(En,bth(1),bth(2),fse,EEG_N,twn,tbn,roi_l);
        mu_c=bandpow_db_roi(Ec,bmu(1),bmu(2),fse,EEG_N,twc,tbc,roi_l);
        mu_n=bandpow_db_roi(En,bmu(1),bmu(2),fse,EEG_N,twn,tbn,roi_l);
        dtheta(s)=th_c-th_n; dmusc(s)=mu_c-mu_n; %#ok<PFOUS>
    end
    % A4(ii): razón armónica de la envolvente EMG (2f/f) — bloque chewing continuo
    fc=fullfile(pac_l,[cases_l{s} suf_ch]);
    if exist(fc,'file')
        E=pop_loadset(fc); fsm=E.srate; ch=emg_l(emg_l<=E.nbchan);
        if ~isempty(ch)
            emg=emg_bilateral(E.data,ch,fsm);
            % envolvente del ritmo masticatorio: rectificar EMG de alta freq y suavizar
            [bh,ah]=butter(4,[20 min(200,fsm/2-1)]/(fsm/2),'bandpass');
            env=abs(filtfilt(bh,ah,emg)); env=detrend(env);
            [pxx,fx]=pwelch(env,round(fsm*8),round(fsm*4),[],fsm);
            % f_chew como pico en 0.8-2.2
            fi=fx>=0.8&fx<=2.2; fxs=fx(fi); pxs=pxx(fi); [~,pk]=max(pxs); f0=fxs(pk);
            pf =band_pk(pxx,fx,f0,0.25); p2f=band_pk(pxx,fx,2*f0,0.35);
            HR(s)=p2f/(pf+eps); %#ok<PFOUS>
        end
    end
    fprintf('  A3/A4 %d/%d %s\n',s,N_CASES,cases_l{s});
end

%% Correlaciones
[rho_A3,p_A3]=corr(dtheta,dmusc,'Type','Spearman','Rows','complete');
fprintf('\n===== A3: PROXY MUSCULAR =====\n');
fprintf('Δtheta(ROI) × Δ(30-40Hz)(ROI): rho=%+.3f p=%.4f (n=%d) → %s\n',...
    rho_A3,p_A3,sum(~isnan(dtheta+dmusc)),tern(p_A3<0.05,'SIG (preocupa)','n.s. (ok)'));

% A4(ii): correlación HR × zMI theta (de S4b)
D=load(fullfile(OUT_STATS,'v1_S4b_PAC_ROI.mat'),'zC_ch');
zth=D.zC_ch(:,1);
[rho_A4,p_A4]=corr(HR,zth,'Type','Spearman','Rows','complete');
fprintf('\n===== A4(ii): SHARPNESS =====\n');
fprintf('HR(envolvente EMG 2f/f): M=%.3f±%.3f\n',mean(HR,'omitnan'),std(HR,'omitnan'));
fprintf('HR × zMI theta: rho=%+.3f p=%.4f (n=%d) → %s\n',...
    rho_A4,p_A4,sum(~isnan(HR+zth)),tern(p_A4<0.05,'SIG (MI escala con forma)','n.s. (MI NO lo explica la forma)'));

%% Guardar
prov=struct('script',mfilename('fullpath'),'date',datestr(now,'yyyy-mm-dd HH:MM:SS'),...
    'matlab',version,'roi',{ROI_CBPT},'frontal_med',{FRONTAL_MED},'temporal',{TEMPORAL},...
    'band_theta',BAND_THETA,'band_musc',BAND_MUSC,'win_late',WIN_LATE);
save(fullfile(OUT_STATS,'S6_artifact_controls.mat'),'t_el','p_el','p_fdr','fm_idx','tp_idx',...
    'labels','dtheta','dmusc','HR','rho_A3','p_A3','rho_A4','p_A4',...
    'n_fm_sig','n_tp_sig','prov','-v7.3');
fprintf('\nGuardado: S6_artifact_controls.mat\n');

%% ===== locales =====
function s=tern(c,a,b), if c, s=a; else, s=b; end, end

function v=bandpow_db_roi(EEG,flo,fhi,fs,n_ch,tw,tb,roi)
    pe=bandpow_db(EEG,flo,fhi,fs,n_ch,tw,tb);  % [1×n_ch]
    v=mean(pe(roi),'omitnan');
end

function pow_db=bandpow_db(EEG,flo,fhi,fs,n_ch,tw,tb)
    nt=EEG.trials; np=EEG.pnts; dc=reshape(EEG.data(1:n_ch,:,:),n_ch,nt*np);
    try, df=eegfilt(dc,fs,flo,fhi);
    catch, [bc,ac]=butter(4,[flo fhi]/(fs/2),'bandpass'); df=filtfilt(bc,ac,dc')'; end
    env=abs(hilbert(df'))'; e3=reshape(env,n_ch,np,nt); p2=e3.^2;
    pw=squeeze(mean(mean(p2(:,tw,:),2),3)); pb=squeeze(mean(mean(p2(:,tb,:),2),3));
    pow_db=(10*log10((pw+eps)./(pb+eps)))';
end

function pk=band_pk(pxx,fx,f0,bw)
    ix=fx>=f0-bw & fx<=f0+bw; if any(ix), pk=max(pxx(ix)); else, pk=eps; end
end

function p_adj=fdr_bh(p)
    p=p(:); n=numel(p); [ps,si]=sort(p); pa=ps.*n./(1:n)';
    for k=n-1:-1:1, pa(k)=min(pa(k),pa(k+1)); end
    p_adj=zeros(n,1); p_adj(si)=pa;
end
