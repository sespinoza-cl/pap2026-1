%% S4b_recompute_stats.m — recomputa P1-P4/Rayleigh/conducta desde el .mat guardado
% (sin re-correr el parfor). Corrige el orden de salida de signrank ([p,h,stats]).
THIS_DIR=fileparts(mfilename('fullpath')); ROOT_V1F=fileparts(THIS_DIR);
run(fullfile(ROOT_V1F,'S0_config.m'));
D=load(fullfile(OUT_STATS,'v1_S4b_PAC_ROI.mat'));
zC=D.zC_ch; zA=D.zA_ch; zCn=D.zC_nc; zAn=D.zA_nc;
MIw=D.MIw_ch; MIc=D.MI_ch; MIn=D.MI_nc; pref=D.pref_ch; BN=D.BAND_NAMES;
nB=numel(BN); nVal=sum(~isnan(zC(:,1)));

fprintf('\n===== P1: PAC theta > null (N=%d) =====\n',nVal);
kC=sum(zC(:,1)>1.96); kA=sum(zA(:,1)>1.96); kBoth=sum(zC(:,1)>1.96 & zA(:,1)>1.96);
pbC=1-binocdf(kC-1,nVal,0.05); pbA=1-binocdf(kA-1,nVal,0.05); pbB=1-binocdf(kBoth-1,nVal,0.05);
pwC=signrank(zC(:,1),0,'tail','right'); pwA=signrank(zA(:,1),0,'tail','right');
fprintf('theta vs CIRC: M_z=%.2f Mdn=%.2f | %d/%d | binom p=%.2e | Wilcoxon(1-sided) p=%.2e\n',...
        mean(zC(:,1)),median(zC(:,1)),kC,nVal,pbC,pwC);
fprintf('theta vs AAFT: M_z=%.2f Mdn=%.2f | %d/%d | binom p=%.2e | Wilcoxon(1-sided) p=%.2e\n',...
        mean(zA(:,1)),median(zA(:,1)),kA,nVal,pbA,pwA);
fprintf('theta supera AMBOS: %d/%d (binom p=%.2e)\n',kBoth,nVal,pbB);

fprintf('\n===== P2: especificidad de banda =====\n');
for b=1:nB
    kcb=sum(zC(:,b)>1.96); kab=sum(zA(:,b)>1.96);
    pcb=signrank(zC(:,b),0,'tail','right'); pab=signrank(zA(:,b),0,'tail','right');
    fprintf('%-6s | circ M_z=%.2f %d/%d p=%.2e || aaft M_z=%.2f %d/%d p=%.2e\n',...
        BN{b},mean(zC(:,b)),kcb,nVal,pcb,mean(zA(:,b)),kab,nVal,pab);
end
fprintf('Friedman zMI x banda: circ p=%.4f | aaft p=%.4f\n',friedman(zC,1,'off'),friedman(zA,1,'off'));

fprintf('\n===== P3: especificidad de ventana (MI Base/Early/Late) =====\n');
for b=1:nB
    ba=squeeze(MIw(b,1,:)); ea=squeeze(MIw(b,2,:)); la=squeeze(MIw(b,3,:));
    fprintf('%-6s | B->E p=%.3f | B->L p=%.3f | E->L p=%.3f | MI: B=%.2e E=%.2e L=%.2e\n',...
        BN{b},signrank(ea,ba),signrank(la,ba),signrank(la,ea),...
        median(ba,'omitnan'),median(ea,'omitnan'),median(la,'omitnan'));
end
dataW=[squeeze(MIw(1,1,:)) squeeze(MIw(1,2,:)) squeeze(MIw(1,3,:))];
fprintf('Friedman theta x ventana: p=%.4f\n',friedman(dataW,1,'off'));

fprintf('\n===== P4: Ch vs Nc (null por condicion) =====\n');
kCn=sum(zCn(:,1)>1.96);
fprintf('zMI theta circ: Ch M=%.2f vs Nc M=%.2f | Wilcoxon p=%.4f (Nc sig %d/%d)\n',...
    mean(zC(:,1)),mean(zCn(:,1),'omitnan'),signrank(zC(:,1),zCn(:,1)),kCn,nVal);
fprintf('zMI theta aaft: Ch M=%.2f vs Nc M=%.2f | Wilcoxon p=%.4f\n',...
    mean(zA(:,1)),mean(zAn(:,1),'omitnan'),signrank(zA(:,1),zAn(:,1)));
fprintf('MI theta abs Ch vs Nc | Wilcoxon p=%.4f\n',signrank(MIc(1,:)',MIn(1,:)'));

fprintf('\n===== Rayleigh fase preferida =====\n');
for b=1:nB
    pv=pref(:,b); pv=pv(~isnan(pv)); nv=numel(pv);
    Rr=abs(mean(exp(1i*pv))); Zr=nv*Rr^2; pr=exp(-Zr);
    fprintf('%-6s R=%.3f Z=%.2f p=%.4f\n',BN{b},Rr,Zr,pr);
end

fprintf('\n===== PAC x conducta (Spearman) =====\n');
Beh=load(FILE_BEH);
fn=fieldnames(Beh); fprintf('campos beh: %s\n',strjoin(fn,', '));

% Actualizar p-valores corregidos en provenance del .mat
prov=D.prov; prov.p1_corrected=struct('k_circ',kC,'k_aaft',kA,'k_both',kBoth,...
    'binom_circ',pbC,'binom_aaft',pbA,'wilcoxon_circ',pwC,'wilcoxon_aaft',pwA);
prov.p4_corrected=struct('p_zc',signrank(zC(:,1),zCn(:,1)),'nc_sig',kCn);
prov.signrank_fix='2026-06-19: corregido orden de salida signrank [p,h]';
save(fullfile(OUT_STATS,'v1_S4b_PAC_ROI.mat'),'-struct','D');  % conserva data
m=matfile(fullfile(OUT_STATS,'v1_S4b_PAC_ROI.mat'),'Writable',true); m.prov=prov;
fprintf('\nProvenance actualizado con p-valores corregidos.\n');
