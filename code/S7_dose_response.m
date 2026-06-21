%% S7_dose_response.m — Dosis-respuesta: intensidad de masticación × efectos (solo casos)
% Aprovecha el diseño: ¿los casos que mastican MÁS muestran MÁS theta/PAC/mejora?
% engagement = masseter Ch-Nc SNR (rhythmicity, de S0c). Todo en orden CASES (1:31).
THIS=fileparts(mfilename('fullpath')); run(fullfile(fileparts(THIS),'S0_config.m'));
E=load(fullfile(OUT_STATS,'chew_engagement_check.mat'),'casd');      eng=E.casd(:);          % Ch-Nc SNR
S6=load(fullfile(OUT_STATS,'S6_artifact_controls.mat'),'dtheta');    dth=S6.dtheta(:);       % Δtheta ROI
P=load(fullfile(OUT_STATS,'v1_S4b_PAC_ROI.mat'),'zC_ch','MIw_ch','f_chew_hz'); 
zmi=P.zC_ch(:,1); miL=squeeze(P.MIw_ch(1,3,:)); fch=P.f_chew_hz(:);
B=load(FILE_BEH,'rt_delta_cas','ies_delta_cas'); drt=B.rt_delta_cas(:); dies=B.ies_delta_cas(:);
sp=@(x,y) deal_corr(x,y);
fprintf('\n===== DOSIS-RESPUESTA (Spearman, casos n=31) =====\n');
pairs={ 'engagement × Δtheta(ROI)',eng,dth;
        'engagement × zMI theta',eng,zmi;
        'engagement × MI theta Late',eng,miL;
        'engagement × ΔRT',eng,drt;
        'engagement × ΔIES',eng,dies;
        'f_chew × Δtheta',fch,dth;
        'f_chew × zMI theta',fch,zmi;
        'f_chew × MI theta Late',fch,miL;
        'Δtheta × zMI theta',dth,zmi;
        'Δtheta × ΔIES',dth,dies};
R=nan(size(pairs,1),2);
for i=1:size(pairs,1)
  [r,p]=corr(pairs{i,2},pairs{i,3},'Type','Spearman','Rows','complete'); R(i,:)=[r p];
  fprintf('%-28s rho=%+.3f p=%.4f %s\n',pairs{i,1},r,p,tern(p<0.05,'*',''));
end
save(fullfile(OUT_STATS,'S7_dose_response.mat'),'pairs','R','eng','dth','zmi','miL','fch','drt','dies');
fprintf('\nGuardado: S7_dose_response.mat\n');
function s=tern(c,a,b), if c,s=a;else,s=b;end,end
function [r,p]=deal_corr(x,y), [r,p]=corr(x,y,'Type','Spearman','Rows','complete'); end
