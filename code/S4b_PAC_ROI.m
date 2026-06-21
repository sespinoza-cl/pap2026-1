%% S4b_PAC_ROI.m  —  PAC principal (ROI-18 θ-interacción), versión a prueba de balas
%
% CORRIGE respecto a S4b_PAC_ROI_OLD_noncompliant.m:
%   [I1] cargaba Analysis_v1\S0_config.m (VIEJA) → ahora SOLO el canónico (R1).
%   [I3] zMI_nc usaba el null de Ch → cada condición tiene su PROPIO null.
%   [R6] EMG inconsistente → criterio único: promedio bilateral 65+66 + fallback (emg_bilateral.m).
%   [R7] sólo circ-shift → DOBLE NULL: circular-shift (5 s) + AAFT (compute_pac_cont.m).
%   [F5] ERP sobre continuo solapado → epoca-primero, resta media-por-época (compute_pac_windows.m).
%   [R3] rng(RNG_SEED) + RandStream substream por sujeto (parfor reproducible).
%   [R4] sidecar provenance JSON.
%
% P1 ¿PAC θ>null? | P2 ¿banda? | P3 ¿ventana? | P4 ¿Ch vs Nc? (null por condición)
% Salida: outputs/stats/v1_S4b_PAC_ROI.mat (+ provenance.json)

%% ── FUENTE ÚNICA DE VERDAD (R1) ──────────────────────────────────────────
THIS_DIR = fileparts(mfilename('fullpath'));
ROOT_V1F = fileparts(THIS_DIR);
run(fullfile(ROOT_V1F, 'S0_config.m'));
rng(RNG_SEED);

if ~exist('SMOKE','var'), SMOKE = false; end
N_RUN = N_CASES; smoke_tag = '';
if SMOKE, N_RUN = min(2, N_CASES); smoke_tag = '_SMOKE'; end
% N_LIMIT: prueba de parfor en pocos sujetos sin marcar como SMOKE (no sobrescribe el .mat full)
if exist('N_LIMIT','var') && ~isempty(N_LIMIT)
    N_RUN = min(N_LIMIT, N_CASES); smoke_tag = sprintf('_PARTEST%d', N_RUN);
end

%% ── Parámetros del PAC (derivados de S0_config) ──────────────────────────
PHASE_HALF  = 0.5;
BANDS_HZ    = {BAND_THETA, BAND_ALPHA, BAND_BETA};
BAND_NAMES  = {'theta','alpha','beta'};
nBands      = numel(BANDS_HZ);
MIN_SHIFT_S = 5.0;
MIN_SHIFT   = round(MIN_SHIFT_S * FS);
WINS_MS     = [WIN_BASE; WIN_EARLY; WIN_LATE];
wS          = round([WIN_EPOCH(1) WIN_EPOCH(2)]/1000 * FS);

fprintf('ROI-18: %s\n', strjoin(ROI_CBPT,', '));
fprintf('θ[%g %g] α[%g %g] β[%g %g] | min_shift=%.1fs(%d) | N_SURR=%d×2 nulls (circ+AAFT)\n',...
        BAND_THETA,BAND_ALPHA,BAND_BETA,MIN_SHIFT_S,MIN_SHIFT,N_SURR);

%% ── f_chew individual ────────────────────────────────────────────────────
tmp = load(FILE_CHEW,'T_freq'); T_freq = tmp.T_freq;
F_chew = nan(N_CASES,1);
for s_ = 1:N_CASES
    ix_ = strcmp(T_freq.Sujeto, CASES{s_});
    if any(ix_), F_chew(s_) = mean([T_freq.Freq_Left(ix_), T_freq.Freq_Right(ix_)],'omitnan'); end
end
fprintf('f_chew: M=%.3f±%.3f Hz (%d/%d con valor)\n',...
        mean(F_chew,'omitnan'),std(F_chew,'omitnan'),sum(~isnan(F_chew)),N_CASES);

%% ── Parámetros broadcast ─────────────────────────────────────────────────
P = struct('cases',{CASES},'fchew',F_chew,'pac',DATA_PAC,...
    'suf_ch',EMG_SUFFIX,'suf_nc',EMG_SUFFIX_NC,'roi',{ROI_CBPT},'emg_chs',EMG_CHANS,...
    'bands',{BANDS_HZ},'nb',MI_BINS,'ns',N_SURR,'minsh',MIN_SHIFT,'ph_half',PHASE_HALF,...
    'wS',wS,'wins',WINS_MS,'seed',RNG_SEED,'ev_ch',40,'ev_nc',30);

%% ── Pool paralelo + EEGLAB en workers ────────────────────────────────────
% El cliente ya tiene el path correcto de EEGLAB (run_S4b → eeglab nogui).
% Replicamos ESE path a los workers; NO usar genpath(eeglab) porque añade
% Fieldtrip/compat/octave que shadowa builtins (fullfile, contains) y atasca parfor.
if ~SMOKE && isempty(gcp('nocreate'))
    try parpool('local', min(12, feature('numcores'))); catch; end
end
cur_path = path();
if ~SMOKE
    try, spmd, addpath(cur_path); end, catch; end
end

%% ── Loop principal (pac_subject por sujeto) ──────────────────────────────
fprintf('Procesando %d sujetos %s...\n', N_RUN, smoke_tag);
res = repmat(pac_subject_empty(nBands), N_RUN, 1);
if SMOKE
    for s = 1:N_RUN, res(s) = pac_subject(s, P); end
else
    parfor s = 1:N_RUN, res(s) = pac_subject(s, P); end
end
fprintf('\n✓ Procesamiento completo.\n');

%% ── Unpack ───────────────────────────────────────────────────────────────
MI_ch=nan(nBands,N_RUN); MI_nc=nan(nBands,N_RUN);
zC_ch=nan(N_RUN,nBands); zC_nc=nan(N_RUN,nBands);
zA_ch=nan(N_RUN,nBands); zA_nc=nan(N_RUN,nBands);
pref_ch=nan(N_RUN,nBands); pref_nc=nan(N_RUN,nBands);
nullC_m=nan(N_RUN,nBands); nullC_s=nan(N_RUN,nBands);
nullA_m=nan(N_RUN,nBands); nullA_s=nan(N_RUN,nBands);
MIw_ch=nan(nBands,3,N_RUN); MIw_nc=nan(nBands,3,N_RUN);
f_chew_hz=nan(N_RUN,1);
emg_used_ch=nan(N_RUN,1); emg_used_nc=nan(N_RUN,1);
emg_snr_ch=nan(N_RUN,1);  emg_snr_nc=nan(N_RUN,1);
n_trials_ch=nan(N_RUN,1); n_trials_nc=nan(N_RUN,1);
for s=1:N_RUN
    r=res(s);
    MI_ch(:,s)=r.MI_ch; MI_nc(:,s)=r.MI_nc;
    zC_ch(s,:)=r.zC_ch; zC_nc(s,:)=r.zC_nc;
    zA_ch(s,:)=r.zA_ch; zA_nc(s,:)=r.zA_nc;
    pref_ch(s,:)=r.pref_ch; pref_nc(s,:)=r.pref_nc;
    nullC_m(s,:)=r.nullC_m; nullC_s(s,:)=r.nullC_s;
    nullA_m(s,:)=r.nullA_m; nullA_s(s,:)=r.nullA_s;
    MIw_ch(:,:,s)=r.MIw_ch; MIw_nc(:,:,s)=r.MIw_nc;
    f_chew_hz(s)=r.f_chew;
    emg_used_ch(s)=r.emg_used_ch; emg_used_nc(s)=r.emg_used_nc;
    emg_snr_ch(s)=r.emg_snr_ch; emg_snr_nc(s)=r.emg_snr_nc;
    n_trials_ch(s)=r.n_trials_ch; n_trials_nc(s)=r.n_trials_nc;
end

%% ── ESTADÍSTICA ──────────────────────────────────────────────────────────
nVal = sum(~isnan(zC_ch(:,1)));
fprintf('\n===== P1: ¿PAC θ > null? (Cases-Chew, N=%d) =====\n', nVal);
kC=sum(zC_ch(:,1)>1.96,'omitnan'); kA=sum(zA_ch(:,1)>1.96,'omitnan');
kBoth=sum(zC_ch(:,1)>1.96 & zA_ch(:,1)>1.96,'omitnan');
pbC=1-binocdf(kC-1,nVal,0.05); pbA=1-binocdf(kA-1,nVal,0.05);
% signrank devuelve [p,h,stats] (p primero, distinto a ttest) — one-sided z>0
pwC=signrank(zC_ch(:,1),0,'tail','right'); pwA=signrank(zA_ch(:,1),0,'tail','right');
fprintf('θ vs CIRC: M_z=%.2f Mdn=%.2f | %d/%d sig | binom p=%.2e | Wilcoxon p=%.2e\n',...
    mean(zC_ch(:,1),'omitnan'),median(zC_ch(:,1),'omitnan'),kC,nVal,pbC,pwC);
fprintf('θ vs AAFT: M_z=%.2f Mdn=%.2f | %d/%d sig | binom p=%.2e | Wilcoxon p=%.2e\n',...
    mean(zA_ch(:,1),'omitnan'),median(zA_ch(:,1),'omitnan'),kA,nVal,pbA,pwA);
fprintf('θ supera AMBOS nulls: %d/%d\n', kBoth,nVal);

fprintf('\n===== P2: ¿Específico de banda? =====\n');
for b=1:nBands
    kcb=sum(zC_ch(:,b)>1.96,'omitnan'); kab=sum(zA_ch(:,b)>1.96,'omitnan');
    pcb=signrank(zC_ch(:,b),0,'tail','right'); pab=signrank(zA_ch(:,b),0,'tail','right');
    fprintf('%-6s | circ M_z=%.2f %d/%d p=%.4f || aaft M_z=%.2f %d/%d p=%.4f\n',...
        BAND_NAMES{b},mean(zC_ch(:,b),'omitnan'),kcb,nVal,pcb,...
        mean(zA_ch(:,b),'omitnan'),kab,nVal,pab);
end
[p_fr_c,~]=friedman(zC_ch,1,'off'); [p_fr_a,~]=friedman(zA_ch,1,'off');
fprintf('Friedman zMI×banda: circ p=%.4f | aaft p=%.4f\n', p_fr_c,p_fr_a);

fprintf('\n===== P3: ¿Específico de ventana? (MI Base/Early/Late) =====\n');
for b=1:nBands
    base=squeeze(MIw_ch(b,1,:)); early=squeeze(MIw_ch(b,2,:)); late=squeeze(MIw_ch(b,3,:));
    pbe=signrank(early,base); pbl=signrank(late,base); pel=signrank(late,early);
    fprintf('%-6s | B→E p=%.3f | B→L p=%.3f | E→L p=%.3f\n',BAND_NAMES{b},pbe,pbl,pel);
end
dataW=[squeeze(MIw_ch(1,1,:)) squeeze(MIw_ch(1,2,:)) squeeze(MIw_ch(1,3,:))];
[p_fw,~]=friedman(dataW,1,'off');
fprintf('Friedman θ × ventana: p=%.4f\n', p_fw);

fprintf('\n===== P4: ¿Difiere Ch vs Nc? (θ, null por condición) =====\n');
p_zc_cond=signrank(zC_ch(:,1),zC_nc(:,1));
p_za_cond=signrank(zA_ch(:,1),zA_nc(:,1));
p_mi_cond=signrank(MI_ch(1,:)',MI_nc(1,:)');
kC_nc=sum(zC_nc(:,1)>1.96,'omitnan');
fprintf('zMI θ circ: Ch M=%.2f vs Nc M=%.2f | Wilcoxon p=%.4f (Nc sig %d/%d)\n',...
    mean(zC_ch(:,1),'omitnan'),mean(zC_nc(:,1),'omitnan'),p_zc_cond,kC_nc,nVal);
fprintf('zMI θ aaft: Ch M=%.2f vs Nc M=%.2f | Wilcoxon p=%.4f\n',...
    mean(zA_ch(:,1),'omitnan'),mean(zA_nc(:,1),'omitnan'),p_za_cond);
fprintf('MI θ abs Ch vs Nc | Wilcoxon p=%.4f\n', p_mi_cond);

%% ── Rayleigh fase preferida ──────────────────────────────────────────────
fprintf('\n===== Rayleigh fase preferida (Cases-Chew) =====\n');
rayl_R=nan(nBands,1); rayl_Z=nan(nBands,1); rayl_p=nan(nBands,1);
for b=1:nBands
    pv=pref_ch(:,b); pv=pv(~isnan(pv)); nv=numel(pv);
    Rr=abs(mean(exp(1i*pv))); Zr=nv*Rr^2; pr=exp(-Zr);
    rayl_R(b)=Rr; rayl_Z(b)=Zr; rayl_p(b)=pr;
    fprintf('%-6s R=%.3f Z=%.2f p=%.4f (n=%d)\n',BAND_NAMES{b},Rr,Zr,pr,nv);
end

%% ── Correlaciones PAC × conducta (A6) ────────────────────────────────────
fprintf('\n===== PAC × conducta (Spearman, Cases-Chew) =====\n');
Beh=load(FILE_BEH);
RT  = pick_field(Beh,{'rt_cas_b2','medrt_cas_b2','rt_cas_ch'});
IES = pick_field(Beh,{'ies_cas_b2','ies_cas_ch'});
corr_tbl={};
if numel(RT)==N_RUN
    [r1,p1]=corr(zC_ch(:,1),RT,'Type','Spearman','Rows','complete');
    [r2,p2]=corr(squeeze(MIw_ch(1,3,:)),RT,'Type','Spearman','Rows','complete');
    fprintf('zMI θ circ × RT  ρ=%+.3f p=%.4f | MI θ Late × RT  ρ=%+.3f p=%.4f\n',r1,p1,r2,p2);
    corr_tbl=[corr_tbl;{'zMIcirc_x_RT',r1,p1};{'MIlate_x_RT',r2,p2}];
end
if numel(IES)==N_RUN
    [r3,p3]=corr(zC_ch(:,1),IES,'Type','Spearman','Rows','complete');
    [r4,p4]=corr(squeeze(MIw_ch(1,3,:)),IES,'Type','Spearman','Rows','complete');
    fprintf('zMI θ circ × IES ρ=%+.3f p=%.4f | MI θ Late × IES ρ=%+.3f p=%.4f\n',r3,p3,r4,p4);
    corr_tbl=[corr_tbl;{'zMIcirc_x_IES',r3,p3};{'MIlate_x_IES',r4,p4}];
end

%% ── Provenance + guardar (R3/R4) ─────────────────────────────────────────
prov=struct('script',mfilename('fullpath'),'date',datestr(now,'yyyy-mm-dd HH:MM:SS'),...
    'matlab',version,'rng_seed',RNG_SEED,'subjects',{CASES(1:N_RUN)},'roi',{ROI_CBPT},...
    'bands',{BANDS_HZ},'band_names',{BAND_NAMES},'phase_half_hz',PHASE_HALF,...
    'min_shift_s',MIN_SHIFT_S,'n_surr',N_SURR,'nulls',{{'circular-shift','AAFT'}},...
    'wins_ms',WINS_MS,'win_epoch_ms',WIN_EPOCH,'N_valid',nVal,...
    'p1',struct('k_circ',kC,'k_aaft',kA,'k_both',kBoth,'binom_circ',pbC,'binom_aaft',pbA,...
                'wilcoxon_circ',pwC,'wilcoxon_aaft',pwA),...
    'p4',struct('p_zc_cond',p_zc_cond,'p_za_cond',p_za_cond,'p_mi_cond',p_mi_cond),'smoke',SMOKE);

ROI_PAC=ROI_CBPT;
out_mat = fullfile(OUT_STATS, ['v1_S4b_PAC_ROI' smoke_tag '.mat']);
save(out_mat,'MI_ch','MI_nc','zC_ch','zC_nc','zA_ch','zA_nc','pref_ch','pref_nc',...
   'nullC_m','nullC_s','nullA_m','nullA_s','MIw_ch','MIw_nc','f_chew_hz',...
   'emg_used_ch','emg_used_nc','emg_snr_ch','emg_snr_nc','n_trials_ch','n_trials_nc',...
   'rayl_R','rayl_Z','rayl_p','corr_tbl','ROI_PAC','BAND_NAMES','BANDS_HZ','prov','-v7.3');
fid=fopen([out_mat(1:end-4) '.provenance.json'],'w'); fprintf(fid,'%s',jsonencode(prov)); fclose(fid);
fprintf('\nGuardado: %s (+ provenance.json)\n', out_mat);

%% ── Helpers locales (sólo fuera de parfor) ───────────────────────────────
function v = pick_field(S,names)
    v=[];
    for i=1:numel(names)
        if isfield(S,names{i}), v=S.(names{i})(:); return; end
    end
end

function r = pac_subject_empty(nB)
    z=@()nan(1,nB);
    r=struct('MI_ch',nan(nB,1),'zC_ch',z(),'zA_ch',z(),'pref_ch',z(),...
        'nullC_m',z(),'nullC_s',z(),'nullA_m',z(),'nullA_s',z(),...
        'MI_nc',nan(nB,1),'zC_nc',z(),'zA_nc',z(),'pref_nc',z(),...
        'MIw_ch',nan(nB,3),'MIw_nc',nan(nB,3),'f_chew',NaN,...
        'emg_used_ch',NaN,'emg_used_nc',NaN,'emg_snr_ch',NaN,'emg_snr_nc',NaN,...
        'n_trials_ch',NaN,'n_trials_nc',NaN);
end
