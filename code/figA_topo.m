% figA_topo.m — Topoplots theta (interacción + por grupo) Paper2 V1, estándar cuadrado/PNG.
% A7 topo_theta_interaction (marca Fpz/AFz FDR) + B8 theta cases/controls.
% Fuente: v1_S2_TF_data.mat (theta_topo_* por canal/sujeto). EEGLAB para chanlocs+topoplot.
clear; clc;
HERE = fileparts(mfilename('fullpath'));
ROOT = fileparts(HERE);
cd(ROOT);
run(fullfile(ROOT,'S0_config.m'));        % paths, CASES, EEG_N, ROI, etc.

if exist('topoplot','file')~=2
    addpath('D:\EEGLAB'); eeglab nogui;
end

EEGr = pop_loadset('filename',[CASES{1} '_Nc_ep.set'],'filepath',DIR_EPOCHS);
chl  = EEGr.chanlocs(1:EEG_N);

S = load(fullfile(ROOT,'data','computed','v1_S2_TF_data.mat'), ...
    'theta_topo_ch','theta_topo_nc','theta_topo_ch_ctr','theta_topo_nc_ctr');

tch = orient_chan(S.theta_topo_ch, EEG_N);  tnc = orient_chan(S.theta_topo_nc, EEG_N);
tchc= orient_chan(S.theta_topo_ch_ctr, EEG_N); tncc= orient_chan(S.theta_topo_nc_ctr, EEG_N);
fprintf('  dims: cases %dx%d  controls %dx%d (EEG_N=%d)\n', size(tch,1),size(tch,2),size(tchc,1),size(tchc,2),EEG_N);

dcas = mean(tch - tnc, 2);
dctr = mean(tchc - tncc, 2);
dint = dcas - dctr;

OUT  = fullfile(ROOT,'outputs','figures');
labels = {chl.labels};
sig  = find(ismember(upper(labels), {'FPZ','AFZ'}));

fprintf('== figA_topo ==\n');
save_topo(dint, '\theta interaction (Cases-Controls)', 'topo_theta_interaction', chl, sig, OUT);
save_topo(dcas, '\theta Cases: Chew-No-chew', 'topoS_theta_cases', chl, [], OUT);
save_topo(dctr, '\theta Controls: Chew-No-chew', 'topoS_theta_controls', chl, [], OUT);
fprintf('Listo topo theta.\n');

% ───────────────────────── local functions ─────────────────────────
function save_topo(data, ttl, fname, chl, sig, OUT)
    data = data(:);
    m = max(abs(data)); if m==0, m=1; end
    f = figure('Visible','off','Color','w','Units','inches','Position',[1 1 3.5 3.5]);
    if ~isempty(sig)
        topoplot(data, chl, 'electrodes','off','whitebk','on', ...
            'maplimits',[-m m], 'emarker2',{sig,'o','k',6,1});
    else
        topoplot(data, chl, 'electrodes','off','whitebk','on','maplimits',[-m m]);
    end
    colormap(rdbu_(64));
    cb = colorbar; cb.Label.String='Chew - No-chew (dB)'; cb.FontName='Arial'; cb.FontSize=9;
    title(ttl,'FontName','Arial','FontSize',12,'FontWeight','normal');
    set(gca,'FontName','Arial');
    set(f,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 3.5 3.5]);
    print(f, fullfile(OUT,[fname '.png']), '-dpng', '-r300');   % cuadrado 1050x1050
    close(f); fprintf('  saved %s.png\n', fname);
end

function Y = orient_chan(X, nch)
    if size(X,1)==nch, Y = X; elseif size(X,2)==nch, Y = X'; else, Y = X; end
end

function cmap = rdbu_(n)
    c = [103 0 31; 178 24 43; 214 96 77; 244 165 130; 253 219 199; ...
         247 247 247; 209 229 240; 146 197 222; 67 147 195; 33 102 172; 5 48 97]/255;
    c = flipud(c);  % azul(-) -> rojo(+)
    xi = linspace(1,size(c,1),n);
    cmap = [interp1(1:size(c,1),c(:,1),xi)' interp1(1:size(c,1),c(:,2),xi)' ...
            interp1(1:size(c,1),c(:,3),xi)'];
end
