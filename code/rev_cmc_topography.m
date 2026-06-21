% rev_cmc_topography.m — Topografía de la CMC β (EEG-EMG masetero), Cases-Chew.
% Distingue CMC genuina (frontocentral/central) de leakage EMG puro (temporal).
% Si CMC β es mayor central > temporal -> apoya origen corticomuscular, no artefacto temporal.
clear; clc;
HERE=fileparts(mfilename('fullpath')); ROOT=fileparts(HERE); cd(ROOT);
run(fullfile(ROOT,'S0_config.m'));
if exist('pop_loadset','file')~=2, addpath('D:\EEGLAB'); eeglab nogui; end

CENTRAL = {'Fz','FCz','Cz','FC1','FC2','C1','C2','F1','F2','AFz','FCC1h','FCC2h'};
TEMPORAL= {'T7','T8','FT7','FT8','TP7','TP8','C5','C6','FC5','FC6'};
betaBand=[13 30];
N=numel(CASES);
cmc_topo=nan(64,N); chl=[];
fprintf('== rev_cmc_topography: CMC beta por canal ==\n');
for s=1:N
    f=fullfile(DATA_PAC,[CASES{s} EMG_SUFFIX]);
    if ~exist(f,'file'), continue; end
    EEG=pop_loadset('filename',[CASES{s} EMG_SUFFIX],'filepath',DATA_PAC);
    if isempty(chl), chl=EEG.chanlocs(1:EEG_N); end
    fs=EEG.srate; X=double(EEG.data); if ndims(X)==3, X=X(:,:); end
    emgc=EMG_CHANS(EMG_CHANS<=size(X,1)); emg=mean(X(emgc,:),1);
    if all(emg==0)||any(isnan(emg)), emg=X(emgc(1),:); end
    win=hann(round(fs)); nov=round(fs/2); nfft=max(256,2^nextpow2(round(fs)));
    for c=1:EEG_N
        [Cxy,ff]=mscohere(detrend(double(X(c,:))),detrend(emg),win,nov,nfft,fs);
        cmc_topo(c,s)=mean(Cxy(ff>=betaBand(1)&ff<=betaBand(2)));
    end
    fprintf('  %s ok\n',CASES{s});
end
ok=~all(isnan(cmc_topo),1);
topo=mean(cmc_topo(:,ok),2,'omitnan');         % 64x1 group mean
labels={chl.labels};
ci=find(ismember(upper(labels),upper(CENTRAL)));
ti=find(ismember(upper(labels),upper(TEMPORAL)));
cen=mean(cmc_topo(ci,ok),1,'omitnan')';  tem=mean(cmc_topo(ti,ok),1,'omitnan')';
[p_ct,~,st]=signrank(cen,tem);
fprintf('\n=== CMC beta topografia (n=%d) ===\n',sum(ok));
fprintf('central (%d el) media=%.4f | temporal (%d el) media=%.4f | Wilcoxon central vs temporal p=%.4f\n',...
    numel(ci),mean(cen),numel(ti),mean(tem),p_ct);
[~,pk]=max(topo); fprintf('Canal pico CMC beta: %s (%.4f)\n',labels{pk},topo(pk));

% topoplot
OUT=fullfile(ROOT,'outputs','figures');
fig=figure('Visible','off','Color','w','Units','inches','Position',[1 1 3.5 3.5]);
topoplot(topo,chl,'electrodes','off','whitebk','on');
colormap(hot); cb=colorbar; cb.Label.String='\beta CMC (EEG-EMG)'; cb.FontName='Arial'; cb.FontSize=9;
title(sprintf('\\beta corticomuscular coherence\n(central %.3f vs temporal %.3f, p=%.3f)',mean(cen),mean(tem),p_ct),...
    'FontName','Arial','FontSize',11,'FontWeight','normal');
set(gca,'FontName','Arial');
set(fig,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 3.5 3.5]);
print(fig,fullfile(OUT,'artifactS_betaCMC_topography.png'),'-dpng','-r300'); close(fig);
fprintf('saved artifactS_betaCMC_topography.png\n');
save(fullfile(ROOT,'outputs','stats','rev_cmc_topography.mat'),'cmc_topo','topo','cen','tem','p_ct','labels','ci','ti');
