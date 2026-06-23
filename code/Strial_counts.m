%% Strial_counts.m — Verificar conteo de trials para el paper
% Extrae el total de correct trials usados en el análisis TF y el promedio
% por sujeto por condición, para verificar el ">4,200 trials; 64.4±17.7
% per participant" reportado en Results §TF.
%
% Fuentes (en orden de prioridad):
%   1. Variables n_trials_* en v1_S2_TF_data.mat (si existen)
%   2. Conteo directo desde archivos de épocas en Data_PAC/Epochs/
%
% Output: outputs/reviewer/Strial_counts_result.txt
%
% Correr desde Analysis_V1_Final/:
%   cd('.../Analysis_V1_Final'); run S0_config; run code/Strial_counts

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end

fprintf('\n=== Strial_counts.m — Verificación conteo de trials ===\n');

%% ── Intento 1: variables n_trials en FILE_TF ─────────────────────────────
fprintf('\nBuscando n_trials_* en %s ...\n', FILE_TF);
vars = whos('-file', FILE_TF);
vnames = {vars.name};
fprintf('Variables en archivo: %s\n', strjoin(vnames, ', '));

n_trial_vars = vnames(~cellfun(@isempty, regexp(vnames, '^n_trial')));
if ~isempty(n_trial_vars)
    fprintf('Variables de trials encontradas: %s\n', strjoin(n_trial_vars, ', '));
    loaded = load(FILE_TF, n_trial_vars{:});
    all_counts = [];
    for v = n_trial_vars
        arr = loaded.(v{1})(:);
        fprintf('  %s: media=%.1f ± %.1f, total=%d, n=%d\n', ...
            v{1}, mean(arr,'omitnan'), std(arr,'omitnan'), sum(arr,'omitnan'), numel(arr));
        all_counts = [all_counts; arr]; %#ok<AGROW>
    end
    grand_total = sum(all_counts, 'omitnan');
    fprintf('\nTotal grand (todas condiciones): %d\n', grand_total);
    fprintf('>4200 check: %s\n', tern(grand_total>4200,'OK','FALLO'));
    per_cond = reshape(all_counts, [], numel(n_trial_vars));
    grand_mean = mean(per_cond(:), 'omitnan');
    grand_std  = std(per_cond(:), 'omitnan');
    fprintf('Media por condición por sujeto: %.1f ± %.1f (paper dice 64.4 ± 17.7)\n', ...
        grand_mean, grand_std);
else
    fprintf('No se encontraron n_trials_* en FILE_TF.\n');
    fprintf('Activando Intento 2: conteo desde archivos de épocas.\n');
end

%% ── Intento 2: contar EEG.trials directo desde epoch .set ───────────────
all_subj = [CASES, CONTROLS];
suffixes = {EP_SUFFIX, EP_SUFFIX_NC};  % _Ch_ep.set, _Nc_ep.set
cond_labels = {'Ch', 'Nc'};
n_subj = numel(all_subj);
n_cond = numel(suffixes);
trial_matrix = nan(n_subj, n_cond);

fprintf('\n--- Intento 2: conteo desde epoch files ---\n');
for s = 1:n_subj
    for c = 1:n_cond
        ep_file = fullfile(DIR_EPOCHS, [all_subj{s} suffixes{c}]);
        if exist(ep_file, 'file')
            info = pop_loadset('filename', [all_subj{s} suffixes{c}], ...
                               'filepath', DIR_EPOCHS, 'loadmode', 'info');
            trial_matrix(s, c) = info.trials;
        end
    end
    fprintf('  %d/%d %s: Ch=%d Nc=%d\n', s, n_subj, all_subj{s}, ...
        trial_matrix(s,1), trial_matrix(s,2));
end

fprintf('\n===== RESUMEN TRIAL COUNTS =====\n');
for c = 1:n_cond
    col = trial_matrix(:, c);
    valid_col = col(~isnan(col));
    fprintf('Condición %-4s: media=%.1f ± %.1f | mediana=%.0f | min=%d max=%d | n=%d\n', ...
        cond_labels{c}, mean(valid_col), std(valid_col), ...
        median(valid_col), min(valid_col), max(valid_col), numel(valid_col));
end
all_vals = trial_matrix(:);
grand_total2 = sum(all_vals, 'omitnan');
per_subj_mean = mean(all_vals,'omitnan');
per_subj_std  = std(all_vals,'omitnan');
fprintf('\nTotal ALL condiciones + sujetos: %d  → >4200? %s\n', ...
    grand_total2, tern(grand_total2>4200,'SÍ (OK)','NO (revisar texto)'));
fprintf('Media por condición-sujeto: %.1f ± %.1f  (paper dice 64.4 ± 17.7)\n', ...
    per_subj_mean, per_subj_std);

%% ── Single-trial CSV: verificar 2127 hits Cases-chewing ─────────────────
st_csv = fullfile(OUT_STATS, 'rev_singletrial.csv');
if exist(st_csv, 'file')
    T_st = readtable(st_csv);
    n_st = height(T_st);
    fprintf('\nSingle-trial CSV (Cases-chewing hits): %d filas  (paper dice 2,127)\n', n_st);
    if n_st ~= 2127
        fprintf('DISCREPANCIA: encontrado=%d vs paper=2127 → revisar\n', n_st);
    else
        fprintf('Coincide con el paper ✓\n');
    end
    % Trials por sujeto
    [grp, ids] = findgroups(T_st.subj);
    n_per_s = splitapply(@numel, T_st.RT, grp);
    fprintf('Single-trial por sujeto: media=%.1f ± %.1f (min=%d max=%d)\n', ...
        mean(n_per_s), std(n_per_s), min(n_per_s), max(n_per_s));
end

%% ── Guardar reporte ──────────────────────────────────────────────────────
out_file = fullfile(OUT_REVIEWER, 'Strial_counts_result.txt');
fid = fopen(out_file, 'w');
fprintf(fid, '=== Trial Counts — P2V1 ===\n');
fprintf(fid, 'Fecha: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf(fid, 'Fuente: epoch files (%s)\n\n', DIR_EPOCHS);
for c = 1:n_cond
    col = trial_matrix(:, c);
    valid_col = col(~isnan(col));
    fprintf(fid, 'Condición %-4s: media=%.1f ± %.1f | total=%d | n_subj=%d\n', ...
        cond_labels{c}, mean(valid_col), std(valid_col), sum(valid_col,'omitnan'), numel(valid_col));
end
fprintf(fid, '\nTotal ALL: %d  → >4200? %s\n', grand_total2, tern(grand_total2>4200,'SÍ','NO'));
fprintf(fid, 'Media por condición-sujeto: %.1f ± %.1f (paper dice 64.4 ± 17.7)\n', ...
    per_subj_mean, per_subj_std);
if exist(st_csv,'file')
    fprintf(fid, '\nSingle-trial CSV (Cases-Ch hits): %d (paper dice 2,127)\n', n_st);
end
fclose(fid);
fprintf('\nGuardado: %s\n', out_file);

function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end
