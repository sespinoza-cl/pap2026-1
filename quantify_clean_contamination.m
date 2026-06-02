%% ============================================================================
%  quantify_clean_contamination.m
%  Control: corre la MISMA métrica de contaminación de quantify_v1_contamination.m
%  pero sobre los datos LIMPIOS (pipeline nuevo). Si la limpieza funcionó, la
%  correlación espacial topo-θ vs topo-EMG y el ratio temporal/frontocentral
%  deben CAER respecto a v1 (ρ=+0.62, ratio=1.72).
%
%  Datos limpios: <DATA_ROOT>\EEG\TF\Group_Cases_Ch_tfraw.mat
%    tfraw_pre_g = [actType(1=total,2=NPL/induced), chan, freq, time, subj, metric(1=pwr,2=itpc)]
%  Correr desde la carpeta Analysis/ (usa S0_paths para las rutas).
%  ============================================================================
clear; clc;
P = S0_paths();

% ── TF limpio (Chew), potencia inducida: actType=2, metric=1 ────────────────
T = load(fullfile(P.dir_tf, 'Group_Cases_Ch_tfraw.mat'));
fn = fieldnames(T); G = [];
for k = 1:numel(fn)
    if ndims(T.(fn{k})) >= 6; G = T.(fn{k}); break; end
end
assert(~isempty(G), 'No encontré la matriz 6D tfraw_pre_g en el archivo.');
tfch = squeeze(G(2,:,:,:,:,1));              % [chan x freq x time x subj]

% ── chanlocs + times (desde un epoch .set) ──────────────────────────────────
sf = dir(fullfile(P.dir_epoch, '*_Ch.set'));
if isempty(sf); sf = dir(fullfile(P.dir_epoch, '*.set')); end
EEG = pop_loadset('filename', sf(1).name, 'filepath', P.dir_epoch, 'loadmode','info');
chanlocs = EEG.chanlocs(1:size(tfch,1));
times    = EEG.times;
frex     = linspace(1, 40, size(tfch,2));    % eje de frecuencia TF (1-40 Hz)

% ── Mismos índices de banda/ventana que en v1 ───────────────────────────────
fTH  = dsearchn(frex', [4 7]');
fEMG = dsearchn(frex', [30 40]');
tLAT = dsearchn(times', [900 1300]');

% ── Topografías (Chew, late) ────────────────────────────────────────────────
topo_theta = squeeze(mean(tfch(:, fTH(1):fTH(2),  tLAT(1):tLAT(2), :), [2 3 4]));
topo_emg   = squeeze(mean(tfch(:, fEMG(1):fEMG(2),tLAT(1):tLAT(2), :), [2 3 4]));

% ── (1) Correlación espacial θ-topo vs EMG-topo ─────────────────────────────
ok = ~isnan(topo_theta) & ~isnan(topo_emg);
[r_sp, p_sp] = corr(topo_theta(ok), topo_emg(ok), 'type', 'Spearman');

% ── (2) Ratio temporal/frontocentral en el mapa θ ───────────────────────────
labels = upper(string({chanlocs.labels}));
front  = ["AFZ","F1","FZ","F2","FC1","FCZ","FC2"];
temp   = ["T7","T8","FT7","FT8","TP7","TP8","F7","F8","FT9","FT10","TP9","TP10"];
iF = find(ismember(labels, front)); iT = find(ismember(labels, temp));
mF = mean(topo_theta(iF),'omitnan'); mT = mean(topo_theta(iT),'omitnan');

fprintf('\n=== CONTROL: misma métrica sobre datos LIMPIOS (pipeline nuevo) ===\n');
fprintf('  N electrodos = %d | N sujetos = %d\n', size(tfch,1), size(tfch,4));
fprintf('  (1) Corr espacial topo-θ vs topo-EMG(30-40 Hz): rho = %+.3f, p = %.3g\n', r_sp, p_sp);
fprintf('  (2) θ frontocentral = %+.3f dB | temporal/lateral = %+.3f dB | ratio T/F = %.2f\n', mF, mT, mT/mF);
fprintf('  --- v1 (M6) era: rho=+0.62 (p<1e-7), ratio T/F=1.72 ---\n');
fprintf('  Si rho y ratio cayeron => la limpieza removió la contaminación muscular.\n\n');
