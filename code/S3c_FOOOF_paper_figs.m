%% S3c_FOOOF_paper_figs.m  —  Figuras separadas FOOOF para el paper
%
%  S3c_FOOOF_broadband_PSD.png        — PSD + fit aperiódico punteado (4 grupos)
%  S3c_FOOOF_periodic_residual.png    — Componente periódico absoluto (4 grupos)
%  S3c_FOOOF_exponent_Cases.png       — Boxplot exponent Cases (NoChew vs Chew)
%  S3c_FOOOF_exponent_Controls.png    — Boxplot exponent Controls (NoChew vs Chew)
%
%  Reglas de estilo:
%   - Sin título en ninguna figura (el nombre de archivo describe el contenido)
%   - Sin N en etiquetas de grupo
%   - Bracket de significancia con margen suficiente (no se superpone)
%
%  Correr desde Analysis_V1_Final/ (requiere FOOOF_Workspace_V1.mat de S3b)

if ~exist('ROOT_FINAL','var')
    if exist('S0_config.m','file'), run('S0_config.m');
    else, error('Correr desde Analysis_V1_Final/'); end
end
fprintf('\n=== S3c_FOOOF_paper_figs.m ===\n');

%% ── Cargar resultados Python ─────────────────────────────────────────────
fw = fullfile(OUT_STATS,'FOOOF_Workspace_V1.mat');
assert(exist(fw,'file')==2,'Falta FOOOF_Workspace_V1.mat — correr S3b primero');
load(fw);

TMP = fullfile(OUT_STATS,'fooof_tmp');
R   = struct();
for gc = 1:2
    grp = {'Cases','Controls'}; g = grp{gc};
    for ic = 1:2
        conds = {'Ch','Nc'}; c = conds{ic};
        fn = fullfile(TMP,sprintf('%s_%s_psd_out.mat',g,c));
        assert(exist(fn,'file')==2,'Falta %s',fn);
        tmp = load(fn);
        R.(g).(c).exp       = tmp.exponents(:);
        R.(g).(c).ap_fits   = tmp.ap_fits;
        R.(g).(c).residuals = tmp.residuals;
        R.(g).(c).psd_log   = tmp.psd_log;
        R.(g).(c).freqs     = tmp.freqs_out(:);
    end
end

f      = R.Cases.Ch.freqs;
f_mask = f >= FIT_RANGE(1) & f <= FIT_RANGE(2);
fv     = f(f_mask);
scale  = 10;   % log10 → dB

mn  = @(M) nanmean(M,2);
sem = @(M) nanstd(M,0,2) ./ sqrt(sum(~isnan(M),2));

clr_cas = COL_CASE;
clr_ctr = COL_CTRL;
clr_nc  = COL_NC;
bds     = {BAND_THETA, BAND_ALPHA, BAND_BETA};
bdl     = {'\theta','\alpha','\beta'};

%% ════════════════════════════════════════════════════════════════════════
%%  FIG 1 — Broadband PSD + fit aperiódico punteado
%% ════════════════════════════════════════════════════════════════════════
psd_cas_ch = scale * R.Cases.Ch.psd_log(f_mask,:);
psd_cas_nc = scale * R.Cases.Nc.psd_log(f_mask,:);
psd_ctr_ch = scale * R.Controls.Ch.psd_log(f_mask,:);
psd_ctr_nc = scale * R.Controls.Nc.psd_log(f_mask,:);
ap_cas_ch  = scale * R.Cases.Ch.ap_fits(f_mask,:);
ap_cas_nc  = scale * R.Cases.Nc.ap_fits(f_mask,:);
ap_ctr_ch  = scale * R.Controls.Ch.ap_fits(f_mask,:);
ap_ctr_nc  = scale * R.Controls.Nc.ap_fits(f_mask,:);

all_mn = [mn(psd_cas_ch);mn(psd_cas_nc);mn(psd_ctr_ch);mn(psd_ctr_nc)];
ylo1 = floor(min(all_mn)-1.5);
yhi1 = ceil( max(all_mn)+1.0);

fig1 = figure('Units','inches','Position',[1 1 6.5 5],'Color','w');
ax1  = axes(fig1); hold(ax1,'on');

for bi=1:3
    xr=bds{bi};
    patch(ax1,[xr(1) xr(2) xr(2) xr(1)],[ylo1 ylo1 yhi1 yhi1],...
        [0.92 0.92 0.92],'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
    text(ax1,mean(xr),yhi1-0.7,bdl{bi},...
        'HorizontalAlignment','center','FontSize',FIG_FS,'Color',[0.5 0.5 0.5]);
end

fill_sem(ax1,fv,mn(psd_ctr_nc),sem(psd_ctr_nc),clr_ctr,0.13);
fill_sem(ax1,fv,mn(psd_ctr_ch),sem(psd_ctr_ch),clr_ctr,0.13);
fill_sem(ax1,fv,mn(psd_cas_nc),sem(psd_cas_nc),clr_nc, 0.18);
fill_sem(ax1,fv,mn(psd_cas_ch),sem(psd_cas_ch),clr_cas,0.18);

hL(1)=plot(ax1,fv,mn(psd_ctr_nc),'-', 'Color',clr_ctr,'LineWidth',FIG_LW);
hL(2)=plot(ax1,fv,mn(psd_ctr_ch),'--','Color',clr_ctr,'LineWidth',FIG_LW);
hL(3)=plot(ax1,fv,mn(psd_cas_nc),'-', 'Color',clr_nc, 'LineWidth',FIG_LW);
hL(4)=plot(ax1,fv,mn(psd_cas_ch),'--','Color',clr_cas,'LineWidth',FIG_LW+0.5);

% Fits aperiódicos — punteados, los 4 grupos
plot(ax1,fv,mn(ap_ctr_nc),':','Color',clr_ctr*0.6+0.4,'LineWidth',1.3,'HandleVisibility','off');
plot(ax1,fv,mn(ap_ctr_ch),':','Color',clr_ctr,         'LineWidth',1.3,'HandleVisibility','off');
plot(ax1,fv,mn(ap_cas_nc),':','Color',clr_nc*0.6+0.4,  'LineWidth',1.3,'HandleVisibility','off');
plot(ax1,fv,mn(ap_cas_ch),':','Color',clr_cas,          'LineWidth',1.3,'HandleVisibility','off');

hDot=plot(ax1,NaN,NaN,'k:','LineWidth',1.3);
legend(ax1,[hL hDot],{'Controls – No Chew','Controls – Chew',...
    'Cases – No Chew','Cases – Chew','Dotted: aperiodic fit'},...
    'Location','southwest','FontSize',FIG_FS-1,'Box','off');

xlim(ax1,FIT_RANGE); ylim(ax1,[ylo1 yhi1]);
xlabel(ax1,'Frequency (Hz)','FontSize',FIG_FS+1);
ylabel(ax1,'Power (dB)','FontSize',FIG_FS+1);
set(ax1,'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS);

exportgraphics(fig1,fullfile(OUT_FIGS,'S3c_FOOOF_broadband_PSD.png'),'Resolution',FIG_DPI);
close(fig1);
fprintf('Fig 1: S3c_FOOOF_broadband_PSD.png\n');

%% ════════════════════════════════════════════════════════════════════════
%%  FIG 2 — Componente periódico absoluto (sin título, eje Y con margen)
%% ════════════════════════════════════════════════════════════════════════
res_cas_ch = scale * R.Cases.Ch.residuals(f_mask,:);
res_cas_nc = scale * R.Cases.Nc.residuals(f_mask,:);
res_ctr_ch = scale * R.Controls.Ch.residuals(f_mask,:);
res_ctr_nc = scale * R.Controls.Nc.residuals(f_mask,:);

mn_r  = [mn(res_cas_ch);mn(res_cas_nc);mn(res_ctr_ch);mn(res_ctr_nc)];
sem_r = [sem(res_cas_ch);sem(res_cas_nc);sem(res_ctr_ch);sem(res_ctr_nc)];
ylo2  = floor((min(mn_r-sem_r)-0.4)*2)/2;
yhi2  = ceil( (max(mn_r+sem_r)+0.6)*2)/2;

fig2 = figure('Units','inches','Position',[1 1 6.5 5],'Color','w');
ax2  = axes(fig2); hold(ax2,'on');

for bi=1:3
    xr=bds{bi};
    patch(ax2,[xr(1) xr(2) xr(2) xr(1)],[ylo2 ylo2 yhi2 yhi2],...
        [0.92 0.92 0.92],'EdgeColor','none','FaceAlpha',0.6,'HandleVisibility','off');
    text(ax2,mean(xr),yhi2-0.15,bdl{bi},...
        'HorizontalAlignment','center','FontSize',FIG_FS,'Color',[0.5 0.5 0.5]);
end
yline(ax2,0,'k--','LineWidth',1.0,'HandleVisibility','off');

fill_sem(ax2,fv,mn(res_ctr_nc),sem(res_ctr_nc),clr_ctr,0.13);
fill_sem(ax2,fv,mn(res_ctr_ch),sem(res_ctr_ch),clr_ctr,0.13);
fill_sem(ax2,fv,mn(res_cas_nc),sem(res_cas_nc),clr_nc, 0.18);
fill_sem(ax2,fv,mn(res_cas_ch),sem(res_cas_ch),clr_cas,0.18);

plot(ax2,fv,mn(res_ctr_nc),'-', 'Color',clr_ctr,'LineWidth',FIG_LW,    'DisplayName','Controls – No Chew');
plot(ax2,fv,mn(res_ctr_ch),'--','Color',clr_ctr,'LineWidth',FIG_LW,    'DisplayName','Controls – Chew');
plot(ax2,fv,mn(res_cas_nc),'-', 'Color',clr_nc, 'LineWidth',FIG_LW,    'DisplayName','Cases – No Chew');
plot(ax2,fv,mn(res_cas_ch),'--','Color',clr_cas,'LineWidth',FIG_LW+0.5,'DisplayName','Cases – Chew');

legend(ax2,'Location','northeast','FontSize',FIG_FS-1,'Box','off');
xlim(ax2,FIT_RANGE); ylim(ax2,[ylo2 yhi2]);
xlabel(ax2,'Frequency (Hz)','FontSize',FIG_FS+1);
ylabel(ax2,'Periodic power above 1/f (dB)','FontSize',FIG_FS+1);
set(ax2,'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS);

exportgraphics(fig2,fullfile(OUT_FIGS,'S3c_FOOOF_periodic_residual.png'),'Resolution',FIG_DPI);
close(fig2);
fprintf('Fig 2: S3c_FOOOF_periodic_residual.png\n');

%% ════════════════════════════════════════════════════════════════════════
%%  FIG 3 & 4 — Exponente aperiódico: Cases y Controls por separado
%% ════════════════════════════════════════════════════════════════════════
exp_cas_ch = R.Cases.Ch.exp;   exp_cas_nc = R.Cases.Nc.exp;
exp_ctr_ch = R.Controls.Ch.exp;exp_ctr_nc = R.Controls.Nc.exp;

valid_cas = ~isnan(exp_cas_ch) & ~isnan(exp_cas_nc);
valid_ctr = ~isnan(exp_ctr_ch) & ~isnan(exp_ctr_nc);

[p_cas,~] = signrank(exp_cas_ch(valid_cas), exp_cas_nc(valid_cas));
[p_ctr,~] = signrank(exp_ctr_ch(valid_ctr), exp_ctr_nc(valid_ctr));

% Escala Y compartida (mismo rango en ambas figuras para comparabilidad)
all_exp = [exp_cas_ch; exp_cas_nc; exp_ctr_ch; exp_ctr_nc];
ylo_e   = floor(min(all_exp(~isnan(all_exp)))*10)/10 - 0.10;
% Dejar margen extra arriba para el bracket + texto de significancia
yhi_e   = ceil( max(all_exp(~isnan(all_exp)))*10)/10 + 0.45;

jit = 0.10;  bw = 0.40;

% Posición del bracket: fijo cerca del tope, dentro del ylim
yb_pct = 0.88;   % fracción del rango Y donde arranca el bracket

configs = {
    'Cases',    exp_cas_nc, exp_cas_ch, valid_cas, N_CASES,    clr_cas, p_cas, 'S3c_FOOOF_exponent_Cases.png';
    'Controls', exp_ctr_nc, exp_ctr_ch, valid_ctr, N_CONTROLS, clr_ctr, p_ctr, 'S3c_FOOOF_exponent_Controls.png';
};

for fi = 1:2
    grp_name  = configs{fi,1};
    exp_nc    = configs{fi,2};
    exp_ch    = configs{fi,3};
    valid     = configs{fi,4};
    nS        = configs{fi,5};
    clr_box   = configs{fi,6};
    p_val     = configs{fi,7};
    fname     = configs{fi,8};

    fig = figure('Units','inches','Position',[1 1 5 6.5],'Color','w');
    ax  = axes(fig); hold(ax,'on');

    draw_box(ax, exp_nc, 1, bw, clr_nc,  0.55);
    draw_box(ax, exp_ch, 2, bw, clr_box, 0.55);

    rng(RNG_SEED);
    for i = 1:nS
        if valid(i)
            line(ax,[1 2],[exp_nc(i) exp_ch(i)],...
                'Color',[0.5 0.5 0.5 0.28],'LineWidth',0.8);
        end
    end

    scatter(ax,1+randn(nS,1)*jit, exp_nc, 26, clr_nc, ...
        'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
    scatter(ax,2+randn(nS,1)*jit, exp_ch, 26, clr_box,...
        'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');

    errorbar(ax,1,nanmean(exp_nc),nanstd(exp_nc)/sqrt(sum(valid)),'k^',...
        'MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',5);
    errorbar(ax,2,nanmean(exp_ch),nanstd(exp_ch)/sqrt(sum(valid)),'k^',...
        'MarkerSize',8,'MarkerFaceColor','k','LineWidth',1.5,'CapSize',5);

    % Bracket de significancia — posicionado con fracción del rango Y
    yrange = yhi_e - ylo_e;
    yb = ylo_e + yb_pct * yrange;
    arm = yrange * 0.025;
    line(ax,[1 1 2 2],[yb-arm yb yb yb-arm],'Color','k','LineWidth',1.2);
    text(ax,1.5, yb + yrange*0.035, p2s(p_val),...
        'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');

    set(ax,'XTick',[1 2],'XTickLabel',{'No Chew','Chew'},...
        'XLim',[0.4 2.6],'YLim',[ylo_e yhi_e],...
        'YGrid','on','GridAlpha',0.22,'GridColor',[0.6 0.6 0.6],...
        'Box','off','TickDir','out','FontName',FIG_FONT,'FontSize',FIG_FS+2);

    ylabel(ax,'Aperiodic Exponent (\chi)','FontSize',FIG_FS+2);

    exportgraphics(fig,fullfile(OUT_FIGS,fname),'Resolution',FIG_DPI);
    close(fig);
    fprintf('Fig %d: %s\n', fi+2, fname);
end

fprintf('\n✓ S3c completado — 4 figuras en outputs/figures/\n\n');

%% ── FUNCIONES LOCALES ─────────────────────────────────────────────────────

function draw_box(ax, data, xc, bw, clr, alpha)
    data = data(~isnan(data));
    if isempty(data), return; end
    q  = quantile(data,[0.25 0.50 0.75]);
    iq = q(3)-q(1);
    wlo= min(data(data >= q(1)-1.5*iq));
    whi= max(data(data <= q(3)+1.5*iq));
    hw = bw/2;
    fill(ax,[xc-hw xc+hw xc+hw xc-hw xc-hw],[q(1) q(1) q(3) q(3) q(1)],...
        clr,'FaceAlpha',alpha,'EdgeColor',clr*0.65,'LineWidth',1.2);
    line(ax,[xc-hw xc+hw],[q(2) q(2)],'Color','w','LineWidth',2.5);
    line(ax,[xc xc],[q(1) wlo],'Color',clr*0.65,'LineWidth',1.2);
    line(ax,[xc xc],[q(3) whi],'Color',clr*0.65,'LineWidth',1.2);
    out = data(data<wlo | data>whi);
    if ~isempty(out)
        scatter(ax,repmat(xc,size(out)),out,18,clr*0.65,...
            'filled','MarkerEdgeColor','none','HandleVisibility','off');
    end
end

function fill_sem(ax, x, mn_v, se_v, clr, alp)
    xi = [x(:)', fliplr(x(:)')];
    yi = [mn_v(:)'-se_v(:)', fliplr(mn_v(:)'+se_v(:)')];
    fill(ax,xi,yi,clr,'FaceAlpha',alp,'EdgeColor','none','HandleVisibility','off');
end

function s = p2s(p)
    if isnan(p),    s='n/a';
    elseif p<0.001, s='***';
    elseif p<0.01,  s='**';
    elseif p<0.05,  s='*';
    else,           s='NS';
    end
end
