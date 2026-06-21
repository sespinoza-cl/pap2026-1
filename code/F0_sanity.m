% F0_sanity.m — verificación de config canónica (no carga datos pesados)
cd(fileparts(fileparts(mfilename('fullpath'))));   % Analysis_V1_Final/
run('S0_config.m');
fprintf('\n==================== F0 SANITY ====================\n');
fprintf('N_CASES=%d  N_CONTROLS=%d\n', N_CASES, N_CONTROLS);
fprintf('ROI_CBPT (%d): %s\n', numel(ROI_CBPT), strjoin(ROI_CBPT,', '));
fprintf('BAND_THETA=[%g %g] ALPHA=[%g %g] BETA=[%g %g] EMG=[%g %g] MUSC=[%g %g]\n',...
        BAND_THETA, BAND_ALPHA, BAND_BETA, BAND_EMG, BAND_MUSC);
fprintf('WIN_EARLY=[%g %g] WIN_LATE=[%g %g] WIN_BASE=[%g %g] WIN_EPOCH=[%g %g] WIN_ANALYSIS=[%g %g]\n',...
        WIN_EARLY, WIN_LATE, WIN_BASE, WIN_EPOCH, WIN_ANALYSIS);
fprintf('FS=%d N_PERM=%d MI_BINS=%d N_SURR=%d RNG_SEED=%d\n', FS, N_PERM, MI_BINS, N_SURR, RNG_SEED);
fprintf('EMG_CHANS=[%s]\n', num2str(EMG_CHANS));
fprintf('--- input files exist? ---\n');
ff = {FILE_CHEW,FILE_TF,FILE_FOOOF,FILE_PAC,FILE_BEH,FILE_CSV};
fn = {'FILE_CHEW','FILE_TF','FILE_FOOOF','FILE_PAC','FILE_BEH','FILE_CSV'};
for i=1:numel(ff), fprintf('  %-10s exist=%d  %s\n', fn{i}, exist(ff{i},'file')>0, ff{i}); end
fprintf('  DATA_PAC dir exists: %d  (%s)\n', exist(DATA_PAC,'dir')>0, DATA_PAC);
d = dir(fullfile(DATA_PAC,'*_clean_emg.set'));
fprintf('  *_clean_emg.set count: %d\n', numel(d));
% chew_metrics sanity
tmp = load(FILE_CHEW,'T_freq'); T_freq = tmp.T_freq;
fc = nan(N_CASES,1);
for s=1:N_CASES, ix=strcmp(T_freq.Sujeto,CASES{s}); if any(ix), fc(s)=mean([T_freq.Freq_Left(ix),T_freq.Freq_Right(ix)],'omitnan'); end; end
fprintf('  f_chew: M=%.3f SD=%.3f  N_with_value=%d/%d\n', mean(fc,'omitnan'), std(fc,'omitnan'), sum(~isnan(fc)), N_CASES);
fprintf('==================================================\n');
