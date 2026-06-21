addpath('D:\EEGLAB'); evalc('eeglab nogui'); close all force;
THIS=fileparts(mfilename('fullpath')); run(fullfile(fileparts(THIS),'S0_config.m'));
cur=path(); try, spmd, addpath(cur); end, catch; end
ALL=[CASES, CONTROLS]; isc=[false(1,N_CASES) true(1,N_CONTROLS)]; n=numel(ALL);
snrCh=nan(n,1); snrNc=nan(n,1);
emg_l=EMG_CHANS; pac_l=DATA_PAC;
parfor i=1:n
  snrCh(i)=chewsnr(fullfile(pac_l,[ALL{i} '_Ch_clean_emg.set']),emg_l); %#ok<PFBNS>
  snrNc(i)=chewsnr(fullfile(pac_l,[ALL{i} '_Nc_clean_emg.set']),emg_l);
end
casd=snrCh(~isc)-snrNc(~isc); ctrd=snrCh(isc)-snrNc(isc);
fprintf('\n===== Masseter chew-band SNR (dB): Ch vs Nc =====\n');
fprintf('CASOS (n=%d):     Ch=%+.2f±%.2f  Nc=%+.2f±%.2f  Ch-Nc=%+.2f±%.2f\n',...
  N_CASES,mean(snrCh(~isc),'omitnan'),std(snrCh(~isc),'omitnan'),mean(snrNc(~isc),'omitnan'),std(snrNc(~isc),'omitnan'),mean(casd,'omitnan'),std(casd,'omitnan'));
fprintf('CONTROLES (n=%d): Ch=%+.2f±%.2f  Nc=%+.2f±%.2f  Ch-Nc=%+.2f±%.2f\n',...
  N_CONTROLS,mean(snrCh(isc),'omitnan'),std(snrCh(isc),'omitnan'),mean(snrNc(isc),'omitnan'),std(snrNc(isc),'omitnan'),mean(ctrd,'omitnan'),std(ctrd,'omitnan'));
p_cas=signrank(casd); p_ctr=signrank(ctrd);
fprintf('Ch>Nc dentro de grupo (Wilcoxon): casos p=%.4g | controles p=%.4g\n',p_cas,p_ctr);
fprintf('Casos con Ch-Nc>0: %d/%d | Controles con Ch-Nc>0: %d/%d\n',sum(casd>0),N_CASES,sum(ctrd>0),N_CONTROLS);
[p_bt]=ranksum(casd,ctrd);
fprintf('Ch-Nc casos vs controles (ranksum): p=%.4g\n',p_bt);
save(fullfile(OUT_STATS,'chew_engagement_check.mat'),'ALL','isc','snrCh','snrNc','casd','ctrd','p_cas','p_ctr','p_bt');
function s=chewsnr(f,emg)
  s=NaN; if ~exist(f,'file'), return; end
  E=pop_loadset(f); fs=E.srate; ch=emg(emg<=E.nbchan); if isempty(ch), return; end
  x=mean(double(E.data(ch,:)),1);
  [b,a]=butter(4,[20 min(200,fs/2-1)]/(fs/2),'bandpass'); env=detrend(abs(filtfilt(b,a,x)));
  [px,fx]=pwelch(env,round(fs*4),round(fs*2),[],fs);
  inb=fx>=1&fx<=2.5; outb=(fx>=0.1&fx<1)|(fx>2.5&fx<=5);
  s=10*log10(mean(px(inb))/mean(px(outb)));
end
