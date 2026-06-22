% figA_topo.m — Topoplots theta (3 contrastes) estilo figures_ok: JET ±1.5, contornos,
% electrodos SIGNIFICATIVOS POR CONTRASTE (FDR), etiqueta de banda, colorbar, CUADRADO.
% Mapas + significancia desde S2c_TF_GroupFigure.mat (topo_* + sig_*_fdr).
%   Cases Ch−Nc: sig_cas_fdr (Fpz,AFz) · Controls: sig_ctr_fdr (0) · Interacción: sig_int_fdr (18).
clear; clc;
HERE = fileparts(mfilename('fullpath')); ROOT = fileparts(HERE); cd(ROOT);
run(fullfile(ROOT,'S0_config.m'));
if exist('topoplot','file')~=2, addpath('D:\EEGLAB'); eeglab nogui; end

EEGr = pop_loadset('filename',[CASES{1} '_Nc_ep.set'],'filepath',DIR_EPOCHS);
chl  = EEGr.chanlocs(1:EEG_N);

S = load(fullfile(ROOT,'outputs','stats','S2c_TF_GroupFigure.mat'), ...
    'topo_cas','topo_ctr','topo_int','sig_cas_fdr','sig_ctr_fdr','sig_int_fdr');

fprintf('== figA_topo (jet, sig por contraste) ==\n');
save_topo(S.topo_cas(:), find(S.sig_cas_fdr(:)), '', 'topoS_theta_cases',    chl, EEG_N);
save_topo(S.topo_ctr(:), find(S.sig_ctr_fdr(:)), '', 'topoS_theta_controls', chl, EEG_N);
save_topo(S.topo_int(:), find(S.sig_int_fdr(:)), '', 'topo_theta_interaction',chl, EEG_N);
fprintf('Listo (cases=%d, ctrl=%d, int=%d electrodos sig).\n', ...
    nnz(S.sig_cas_fdr), nnz(S.sig_ctr_fdr), nnz(S.sig_int_fdr));

% ───────────────────────── local functions ─────────────────────────
function save_topo(data, sigidx, bandlbl, fname, chl, EEG_N)
    data = data(:); if numel(data) > EEG_N, data = data(1:EEG_N); end
    OUT = fullfile(fileparts(fileparts(mfilename('fullpath'))),'outputs','figures');
    f = figure('Visible','off','Color','w','Units','inches','Position',[1 1 3.5 3.5]);
    args = {'electrodes','off','whitebk','on','numcontour',6,'maplimits',[-1.5 1.5]};
    if ~isempty(sigidx)
        args = [args, {'emarker2',{sigidx(:)','.','k',14,1}}];
    end
    topoplot(data, chl, args{:});
    colormap(jet);
    cb = colorbar; cb.Label.String='Chew - No-chew (dB)'; cb.FontName='Arial'; cb.FontSize=9;
    cb.Ticks = -1.5:0.5:1.5;
    text(-0.55, -0.62, bandlbl, 'FontName','Arial','FontSize',12,'FontWeight','bold');
    set(gca,'FontName','Arial');
    set(f,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 3.5 3.5]);
    print(f, fullfile(OUT,[fname '.png']), '-dpng', '-r300');
    close(f); fprintf('  saved %s.png (%d sig)\n', fname, numel(sigidx));
end
