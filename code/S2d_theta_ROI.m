%% S2d_theta_ROI.m — Derivación CANÓNICA y reproducible del ROI (corrige I4)
% PRIMARIO (canónico): ROI = electrodos con interacción Cases×Cond FDR<0.05 sobre
%   theta_topo (potencia theta NPL/Morlet en WIN_LATE, dB) precomputada en FILE_TF
%   (v1_S2_TF_data.mat) — EXACTAMENTE el método de S2c_TF_GroupFigure.m que generó
%   los 18 canónicos. Esta es la fuente documentada del ROI_CBPT de S0_config.
% SECUNDARIO (robustez): ROI por método de épocas (eegfilt+Hilbert², bandpow_db).
%   Sirve para mostrar que el ROI frontal es robusto al método (no para el canónico).
%
% Salida: outputs/stats/ROI_canonical.mat (.LABELS .IDX prov + robustez) + .txt

THIS_DIR = fileparts(mfilename('fullpath'));
ROOT_V1F = fileparts(THIS_DIR);
run(fullfile(ROOT_V1F, 'S0_config.m'));
fprintf('Banda theta canónica: [%g %g] Hz\n', BAND_THETA);

%% chanlocs
EEG_ref  = pop_loadset(fullfile(DIR_EPOCHS,[CASES{1} '_Nc_ep.set']));
labels   = {EEG_ref.chanlocs(1:EEG_N).labels};

%% ===== PRIMARIO: theta_topo Morlet (FILE_TF) → interacción FDR =====
S = load(FILE_TF,'theta_topo_ch','theta_topo_nc','theta_topo_ch_ctr','theta_topo_nc_ctr');
dcas = S.theta_topo_ch     - S.theta_topo_nc;       % [31×64]
dctr = S.theta_topo_ch_ctr - S.theta_topo_nc_ctr;   % [15×64]
[~,p_int_M] = ttest2(dcas,dctr);
[~,p_cas_M] = ttest(dcas);
sig_int_M = fdr_bh(p_int_M(:))<0.05;
sig_cas_M = fdr_bh(p_cas_M(:))<0.05;
ROI_IDX    = find(sig_int_M);
ROI_LABELS = labels(ROI_IDX);
theta_topo_int_M = median(dcas,1)-median(dctr,1);

fprintf('\n===== PRIMARIO (Morlet theta_topo, método S2c) =====\n');
fprintf('ROI interacción θ FDR: N=%d → %s\n', numel(ROI_LABELS), strjoin(ROI_LABELS,', '));
fprintf('Casos Ch>Nc θ FDR: %s\n', strjoin(labels(sig_cas_M),', '));

CANON18 = {'Fp1','AF7','AF3','F1','F3','F5','F7','FC1','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','FC2','FCz'};
fprintf('vs ROI-18 S0_config: +[%s]  -[%s]\n', ...
    strjoin(setdiff(ROI_LABELS,CANON18),','), strjoin(setdiff(CANON18,ROI_LABELS),','));

%% ===== SECUNDARIO: método épocas (eegfilt+Hilbert²) — robustez =====
bands = struct('name',{'theta','alpha','beta'},'range',{BAND_THETA,BAND_ALPHA,BAND_BETA});
NB=numel(bands);
bp_cc=z3(N_CASES,EEG_N,NB); bp_cn=z3(N_CASES,EEG_N,NB);
bp_tc=z3(N_CONTROLS,EEG_N,NB); bp_tn=z3(N_CONTROLS,EEG_N,NB);
fprintf('\nSecundario (épocas)... casos\n');
for s=1:N_CASES
    [ch,nc]=epoch_multibandpow(DIR_EPOCHS,CASES{s},bands,WIN_LATE,WIN_BASE,EEG_N);
    for b=1:NB, bp_cc(s,:,b)=ch(b,:); bp_cn(s,:,b)=nc(b,:); end
end
for s=1:N_CONTROLS
    [ch,nc]=epoch_multibandpow(DIR_EPOCHS,CONTROLS{s},bands,WIN_LATE,WIN_BASE,EEG_N);
    for b=1:NB, bp_tc(s,:,b)=ch(b,:); bp_tn(s,:,b)=nc(b,:); end
end
sig_int_ep=cell(NB,1);
for b=1:NB
    dc=bp_cc(:,:,b)-bp_cn(:,:,b); dt=bp_tc(:,:,b)-bp_tn(:,:,b);
    [~,pe]=ttest2(dc,dt); sig_int_ep{b}=fdr_bh(pe(:))<0.05;
    fprintf('  %-6s interacción FDR (épocas): %d\n', bands(b).name, sum(sig_int_ep{b}));
end
ROI_ep = labels(sig_int_ep{1});
fprintf('ROI θ épocas: N=%d → %s\n', numel(ROI_ep), strjoin(ROI_ep,', '));
fprintf('Solapamiento Morlet∩épocas: %d electrodos\n', numel(intersect(ROI_LABELS,ROI_ep)));

%% Guardar
prov=struct('script',mfilename('fullpath'),'date',datestr(now,'yyyy-mm-dd HH:MM:SS'),...
    'matlab',version,'band_theta',BAND_THETA,'win_late',WIN_LATE,'win_base',WIN_BASE,...
    'primary_source','theta_topo (Morlet NPL dB, WIN_LATE) en v1_S2_TF_data.mat; interacción ttest2 + FDR-BH<0.05 (método S2c)',...
    'secondary_source','épocas eegfilt+Hilbert² bandpow_db (robustez)',...
    'roi_ep_labels',{ROI_ep},'overlap_ep',numel(intersect(ROI_LABELS,ROI_ep)));
save(fullfile(OUT_STATS,'ROI_canonical.mat'),'ROI_IDX','ROI_LABELS','labels',...
    'sig_int_M','sig_cas_M','theta_topo_int_M','sig_int_ep','ROI_ep','prov','-v7.3');
fid=fopen(fullfile(OUT_STATS,'ROI_canonical.txt'),'w');
fprintf(fid,'ROI canónico (PRIMARIO, Morlet theta_topo, interacción θ %g-%g FDR<0.05)\n',BAND_THETA);
fprintf(fid,'%s\n\nN=%d | %s\n',strjoin(ROI_LABELS,' '),numel(ROI_LABELS),prov.primary_source);
fprintf(fid,'Robustez (épocas eegfilt+Hilbert): N=%d → %s | solapamiento=%d\n',...
    numel(ROI_ep),strjoin(ROI_ep,' '),prov.overlap_ep);
fclose(fid);
fprintf('\nGuardado: ROI_canonical.mat + .txt\n');

%% ===== locales =====
function a=z3(n,m,k), a=zeros(n,m,k); end

function [pw_ch,pw_nc]=epoch_multibandpow(dir_ep,subj,bands,win_ms,base_ms,n_ch)
    nb=numel(bands); pw_ch=nan(nb,n_ch); pw_nc=nan(nb,n_ch);
    f_ch=fullfile(dir_ep,[subj '_Ch_ep.set']); f_nc=fullfile(dir_ep,[subj '_Nc_ep.set']);
    if ~exist(f_ch,'file')||~exist(f_nc,'file'), warning('faltan épocas %s',subj); return; end
    E_ch=pop_loadset(f_ch); E_nc=pop_loadset(f_nc); fs=E_ch.srate;
    tw_c=E_ch.times>=win_ms(1)&E_ch.times<=win_ms(2); tw_n=E_nc.times>=win_ms(1)&E_nc.times<=win_ms(2);
    tb_c=E_ch.times>=base_ms(1)&E_ch.times<=base_ms(2); tb_n=E_nc.times>=base_ms(1)&E_nc.times<=base_ms(2);
    for b=1:nb
        pw_ch(b,:)=bandpow_db(E_ch,bands(b).range(1),bands(b).range(2),fs,n_ch,tw_c,tb_c);
        pw_nc(b,:)=bandpow_db(E_nc,bands(b).range(1),bands(b).range(2),fs,n_ch,tw_n,tb_n);
    end
end

function pow_db=bandpow_db(EEG,flo,fhi,fs,n_ch,tw,tb)
    nt=EEG.trials; np=EEG.pnts; dc=reshape(EEG.data(1:n_ch,:,:),n_ch,nt*np);
    try, df=eegfilt(dc,fs,flo,fhi);
    catch, [bc,ac]=butter(4,[flo fhi]/(fs/2),'bandpass'); df=filtfilt(bc,ac,dc')'; end
    env=abs(hilbert(df'))'; e3=reshape(env,n_ch,np,nt); p2=e3.^2;
    pw=squeeze(mean(mean(p2(:,tw,:),2),3)); pb=squeeze(mean(mean(p2(:,tb,:),2),3));
    pow_db=(10*log10((pw+eps)./(pb+eps)))';
end

function p_adj=fdr_bh(p)
    p=p(:); n=numel(p); [ps,si]=sort(p); pa=ps.*n./(1:n)';
    for k=n-1:-1:1, pa(k)=min(pa(k),pa(k+1)); end
    p_adj=zeros(n,1); p_adj(si)=pa;
end
