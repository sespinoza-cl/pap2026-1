%% ============================================================================
%  quantify_v1_contamination.m
%  Cuantifica la contaminación EMG en el análisis tiempo-frecuencia de la v1
%  (preprocesamiento M6). Da DOS números para la carta al revisor:
%    (1) Correlación espacial entre la topografía del "efecto θ" (Chew, late)
%        y una topografía de ALTA FRECUENCIA (banda muscular). Alta corr => el
%        mapa θ está espacialmente impulsado por músculo (contaminación).
%    (2) Ratio de potencia θ en electrodos temporales/laterales vs frontocentrales.
%        >1 => el "θ" es máximo sobre los músculos, no en la línea media (neural).
%
%  Datos v1: D:\Exp2\Paper_plots\For plots\Casos_tbch1_tfraw.mat
%            tfraw = [actType(1=total,2=induced), chan, freq, time, subj]
%            frex = linspace(1,40,400); times = EEG.times
%  ============================================================================
clear; clc;
DIR = 'D:\Exp2\Paper_plots\For plots\';

% ── TF (Chew), componente inducida (dim1 = 2) ───────────────────────────────
T    = load([DIR 'Casos_tbch1_tfraw.mat']);
tfch = squeeze(T.tfraw(2,:,:,:,:));          % [chan x freq x time x subj]

% ── chanlocs + times (desde EEG_tbnc1.mat o chanlocs.mat) ───────────────────
EEG = [];
S = load([DIR 'EEG_tbnc1.mat']);
fn = fieldnames(S);
for k = 1:numel(fn)
    if isstruct(S.(fn{k})) && isfield(S.(fn{k}),'chanlocs'); EEG = S.(fn{k}); break; end
end
if isempty(EEG) || ~isfield(EEG,'times')
    error('No encontré EEG.times/chanlocs en EEG_tbnc1.mat — revisar nombre de variable.');
end
chanlocs = EEG.chanlocs;
times    = EEG.times;
frex     = linspace(1, 40, 400);             % eje de frecuencia v1 (ver Code_for_Figures.m)

% ── Índices de banda y ventana (late θ, igual que la figura v1) ─────────────
fTH  = dsearchn(frex', [4 7]');              % theta
fEMG = dsearchn(frex', [30 40]');            % banda muscular (alta frecuencia)
tLAT = dsearchn(times', [900 1300]');        % ventana late

% ── Topografías (Chew, late) promediando freq×time×sujetos ──────────────────
topo_theta = squeeze(mean(tfch(:, fTH(1):fTH(2),  tLAT(1):tLAT(2), :), [2 3 4]));
topo_emg   = squeeze(mean(tfch(:, fEMG(1):fEMG(2),tLAT(1):tLAT(2), :), [2 3 4]));

% ── (1) Correlación espacial θ-topo vs EMG-topo ─────────────────────────────
ok = ~isnan(topo_theta) & ~isnan(topo_emg);
[r_sp, p_sp] = corr(topo_theta(ok), topo_emg(ok), 'type', 'Spearman');

% ── (2) Ratio temporal/frontocentral en el mapa θ ───────────────────────────
labels = upper(string({chanlocs.labels}));
front  = ["AFZ","F1","FZ","F2","FC1","FCZ","FC2"];
temp   = ["T7","T8","FT7","FT8","TP7","TP8","F7","F8","FT9","FT10","TP9","TP10"];
iF = find(ismember(labels, front));
iT = find(ismember(labels, temp));
mF = mean(topo_theta(iF),'omitnan');
mT = mean(topo_theta(iT),'omitnan');
ratio_TF = mT / mF;

% ── Reporte ─────────────────────────────────────────────────────────────────
fprintf('\n=== CONTAMINACIÓN EMG — análisis TF v1 (Casos/Chew, θ 4-7 Hz, late 900-1300 ms) ===\n');
fprintf('  N electrodos = %d | N sujetos = %d\n', size(tfch,1), size(tfch,4));
fprintf('  (1) Corr espacial topo-θ vs topo-EMG(30-40 Hz): rho = %+.3f, p = %.3g\n', r_sp, p_sp);
fprintf('      -> rho alto/positivo = el mapa theta es espacialmente el mapa muscular.\n');
fprintf('  (2) Potencia θ:  frontocentral = %+.3f dB | temporal/lateral = %+.3f dB\n', mF, mT);
fprintf('      Ratio temporal/frontocentral = %.2f  (>1 = maximo sobre musculo)\n', ratio_TF);
fprintf('================================================================\n');
fprintf('  Pasame estos dos numeros y los agrego a la carta del revisor.\n\n');
