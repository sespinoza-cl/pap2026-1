% rev_beta_cmc.m — Coherencia corticomuscular (CMC) clásica EEG-EMG (masetero), banda β.
% Responde al panel D1/DA-2: convierte el "β-PAC análogo a CMC" en evidencia directa de CMC.
% CMC = coherencia fase-fase (magnitude-squared coherence) EEG×EMG; clásica en β (13-30 Hz).
clear; clc;
HERE = fileparts(mfilename('fullpath')); ROOT = fileparts(HERE); cd(ROOT);
run(fullfile(ROOT,'S0_config.m'));
if exist('pop_loadset','file')~=2, addpath('D:\EEGLAB'); eeglab nogui; end

% ROI frontocentral (igual que el usado para β-PAC): ROI_CBPT
roiLabels = ROI_CBPT;
betaBand = [13 30]; thetaBand=[4 7]; alphaBand=[8 13];

zC = load(fullfile(ROOT,'outputs','stats','v1_S4b_PAC_ROI.mat'),'zC_ch','BAND_NAMES');
zbeta_pac = zC.zC_ch(:,3);   % zMI beta (PAC)
ztheta_pac= zC.zC_ch(:,1);

N = numel(CASES);
bcmc = nan(N,1); tcmc=nan(N,1); acmc=nan(N,1); cl95=nan(N,1); pkf=nan(N,1);
fprintf('== rev_beta_cmc: CMC EEG-EMG (masetero) banda beta ==\n');
for s=1:N
    f = fullfile(DATA_PAC, [CASES{s} EMG_SUFFIX]);
    if ~exist(f,'file'), fprintf('  [skip] %s sin set\n',CASES{s}); continue; end
    EEG = pop_loadset('filename',[CASES{s} EMG_SUFFIX],'filepath',DATA_PAC);
    fs = EEG.srate;
    X = double(EEG.data);
    if ndims(X)==3, X = X(:,:); end           % concatenar epochs si las hay
    labels = {EEG.chanlocs.labels};
    % EMG = promedio de canales 65/66 (fallback al vivo)
    emgc = EMG_CHANS(EMG_CHANS<=size(X,1));
    emg = mean(X(emgc,:),1);
    if all(emg==0) || any(isnan(emg)), emg = X(emgc(1),:); end
    % ROI EEG indices
    [~,roi] = ismember(upper(roiLabels), upper(labels)); roi = roi(roi>0);
    % coherencia por canal ROI, promedio
    win = hann(round(fs)); nov = round(fs/2); nfft = max(256,2^nextpow2(round(fs)));
    Csum = 0; nseg = 0; ff=[];
    for c = roi(:)'
        [Cxy,ff] = mscohere(detrend(double(X(c,:))), detrend(emg), win, nov, nfft, fs);
        Csum = Csum + Cxy;
    end
    C = Csum / numel(roi);
    % nº segmentos para CL: L = floor((Nsamp-nov)/(win-nov))
    L = floor((size(X,2)-nov)/(numel(win)-nov));
    cl95(s) = 1-(0.05)^(1/max(L-1,1));
    inB = ff>=betaBand(1) & ff<=betaBand(2);
    inT = ff>=thetaBand(1)& ff<=thetaBand(2);
    inA = ff>=alphaBand(1)& ff<=alphaBand(2);
    bcmc(s)=mean(C(inB)); tcmc(s)=mean(C(inT)); acmc(s)=mean(C(inA));
    [~,ip]=max(C(ff>=5 & ff<=45)); fsub=ff(ff>=5&ff<=45); pkf(s)=fsub(ip);
    fprintf('  %-7s fs=%d  CMCβ=%.3f (CL95=%.3f %s) θ=%.3f α=%.3f peak=%.1fHz\n', ...
        CASES{s}, fs, bcmc(s), cl95(s), ternary(bcmc(s)>cl95(s),'SIG','ns'), tcmc(s),acmc(s),pkf(s));
end

ok = ~isnan(bcmc);
nsig = sum(bcmc(ok) > cl95(ok));
fprintf('\n=== RESUMEN CMC β (Cases-Chew, n=%d) ===\n', sum(ok));
fprintf('CMCβ media=%.3f  mediana=%.3f  | sig>CL95: %d/%d\n', ...
    mean(bcmc(ok)), median(bcmc(ok)), nsig, sum(ok));
fprintf('CMC por banda (media): θ=%.3f α=%.3f β=%.3f  (β>θ esperado si motor)\n', ...
    mean(tcmc(ok)),mean(acmc(ok)),mean(bcmc(ok)));
[rb,pb]=corr(bcmc(ok), zbeta_pac(ok),'type','Spearman','rows','complete');
[rt,pt]=corr(bcmc(ok), ztheta_pac(ok),'type','Spearman','rows','complete');
fprintf('CMCβ × zMI_β(PAC): rho=%.3f p=%.3f  | CMCβ × zMI_θ(PAC): rho=%.3f p=%.3f\n', rb,pb,rt,pt);
save(fullfile(ROOT,'outputs','stats','rev_beta_cmc.mat'),'bcmc','tcmc','acmc','cl95','pkf','CASES','nsig');
fprintf('Guardado outputs/stats/rev_beta_cmc.mat\n');

function o=ternary(c,a,b), if c, o=a; else, o=b; end, end
