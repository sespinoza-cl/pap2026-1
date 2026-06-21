%% S8_controls_PAC.m — Control negativo: PAC en CONTROLES (no mastican) bloque Ch
% Paralelo a S4b. Como los controles no mastican, su fase EMG en f_chew±0.5 es ruido →
% zMI theta debe estar en el null. Usa f_chew grupal de casos (1.547 Hz) como referencia.
THIS=fileparts(mfilename('fullpath')); run(fullfile(fileparts(THIS),'S0_config.m')); rng(RNG_SEED);
FCHEW_GRP=1.547;                  % media de casos (controles no tienen f_chew individual)
BANDS={BAND_THETA,BAND_ALPHA,BAND_BETA}; MINSH=round(5*FS);
ctrl=CONTROLS; pac=DATA_PAC; suf=EMG_SUFFIX; roi=ROI_CBPT; emgc=EMG_CHANS; nb=MI_BINS; ns=N_SURR; ph=0.5;
cur=path(); if isempty(gcp('nocreate')), try parpool('local',min(12,feature('numcores'))); catch; end; end
try, spmd, addpath(cur); end, catch; end
nC=numel(ctrl); MI=nan(3,nC); zC=nan(nC,3);
parfor s=1:nC
  st=RandStream('Threefry','Seed',42); st.Substream=s;
  f=fullfile(pac,[ctrl{s} suf]); if ~exist(f,'file'), continue; end
  E=pop_loadset(f); fs=E.srate; ch=emgc(emgc<=E.nbchan); if isempty(ch), continue; end %#ok<PFBNS>
  emg=emg_bilateral(E.data,ch,fs);
  lo=max(0.3,FCHEW_GRP-ph); hi=min(fs/2-0.5,FCHEW_GRP+ph);
  [b,a]=butter(4,[lo hi]/(fs/2),'bandpass'); phase=angle(hilbert(filtfilt(b,a,emg)));
  n_eeg=min(64,E.nbchan); [~,ri]=ismember(roi,{E.chanlocs(1:n_eeg).labels}); ri(ri==0)=[];
  if isempty(ri), continue; end
  rs=mean(double(E.data(ri,:)),1);
  R=compute_pac_cont(phase,rs,BANDS,fs,nb,ns,MINSH,st);
  MI(:,s)=R.mi; zC(s,:)=R.zc'; %#ok<PFOUS>
  fprintf('  %d/%d %s\n',s,nC,ctrl{s});
end
k=sum(zC(:,1)>1.96,'omitnan'); nv=sum(~isnan(zC(:,1)));
pw=signrank(zC(:,1),0,'tail','right'); pb=1-binocdf(k-1,nv,0.05);
fprintf('\n===== CONTROL NEGATIVO: PAC theta en CONTROLES (Ch, no mastican) =====\n');
fprintf('zMI theta: M=%.2f Mdn=%.2f | %d/%d>1.96 | binom p=%.3f | Wilcoxon p=%.3f\n',...
  mean(zC(:,1),'omitnan'),median(zC(:,1),'omitnan'),k,nv,pb,pw);
fprintf('(Comparar con CASOS: M_z=13.32, 22/31. Esperado controles ~null.)\n');
save(fullfile(OUT_STATS,'S8_controls_PAC.mat'),'MI','zC','ctrl','FCHEW_GRP','k','nv','pw','pb');
fprintf('Guardado: S8_controls_PAC.mat\n');
