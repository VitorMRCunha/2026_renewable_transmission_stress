%% paper_figures_final_v4.m
% Reads CSV files only; no OPF is rerun.

clear; close all; clc;

rootdir = fileparts(mfilename('fullpath'));
if isempty(rootdir), rootdir = pwd; end

dir_intact = fullfile(rootdir,'results_paper','intact');
dir_n1det  = fullfile(rootdir,'results_paper','deterministic_n1');
dir_n1prob = fullfile(rootdir,'results_paper','probabilistic_n1');
dir_peak   = fullfile(rootdir,'results_paper','peak_sensitivity');
outdir     = fullfile(rootdir,'results_paper','final_figures');
if ~exist(outdir,'dir'), mkdir(outdir); end

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',9);
set(groot,'defaultAxesLineWidth',0.8);
set(groot,'defaultLineLineWidth',1.35);
set(groot,'defaultAxesBox','on');
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultLegendBox','off');

scenario_order = ["S1_Low_Conc","S2_High_Conc","S3_Low_Dist","S4_High_Dist"];
scenario_short = ["S1","S2","S3","S4"];

Trep  = readtable(required_file(dir_intact,'convergence_representative_branches.csv'));
Tint  = readtable(required_file(dir_intact,'convergence_factorial_interaction.csv'));
Tsys  = readtable(required_file(dir_intact,'convergence_system_metrics.csv'));
Trank = readtable(required_file(dir_intact,'convergence_rank_stability.csv'));

Tdet   = readtable(required_file(dir_n1det,'n1_contingency_summary_all.csv'));
Tshort = readtable(required_file(dir_n1det,'n1_critical_contingency_shortlist.csv'));

Tn1target = readtable(required_file(dir_n1prob,'n1_prob_target_branch_metrics.csv'));
Tn1pair   = readtable(required_file(dir_n1prob,'n1_prob_paired_scenario_differences.csv'));

Tpeak = readtable(required_file(dir_peak,'peak_sensitivity_target_branches.csv'));

%% Figure 1: Monte Carlo convergence
f = figure('Units','centimeters','Position',[2 2 17.5 8.2]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl,1); hold(ax,'on');
D = Trep(Trep.branch_idx==38 & abs(Trep.threshold-0.8)<1e-9,:);
markers = {'o','s','^','d'};
for i=1:numel(scenario_order)
    Q = D(string(D.scenario)==scenario_order(i),:);
    Q = sortrows(Q,'N_requested');
    errorbar(ax,Q.N_requested,Q.LCP_pct,Q.LCP_pct-Q.LCP_ci_low_pct, ...
        Q.LCP_ci_high_pct-Q.LCP_pct,'-','DisplayName',scenario_short(i), ...
        'Marker',markers{i},'MarkerSize',4);
end
xlabel(ax,'Number of trajectories, N');
ylabel(ax,'Branch 26--30 LCP(80%) [%]');
grid(ax,'on'); legend(ax,'Location','best');
xline(ax,1600,':','HandleVisibility','off');
title(ax,'(a) LCP convergence');

ax = nexttile(tl,2); hold(ax,'on');
Tint = sortrows(Tint,'N_requested');
errorbar(ax,Tint.N_requested,Tint.interaction_pp,Tint.interaction_pp-Tint.ci_low_pp, ...
    Tint.ci_high_pp-Tint.interaction_pp,'-o','MarkerSize',4);
yline(ax,0,'--');
xlabel(ax,'Number of trajectories, N');
ylabel(ax,'Interaction [percentage points]');
grid(ax,'on');
xline(ax,1600,':','HandleVisibility','off');
title(ax,'(b) Penetration–siting interaction');
export_figure(f,outdir,'Fig_MC_Convergence');

%% Figure 2: Structural versus renewable-sensitive branches
f = figure('Units','centimeters','Position',[2 2 17.5 8.2]);
tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for p=1:2
    if p==1
        br=155; brname='94--100'; ylims=[0 105];
    else
        br=38; brname='26--30'; ylims=[0 13];
    end
    ax=nexttile(tl,p); hold(ax,'on');
    Q=Trep(Trep.N_requested==2000 & Trep.branch_idx==br & abs(Trep.threshold-0.8)<1e-9,:);
    Q=order_scenarios(Q,scenario_order);
    B=bar(ax,categorical(scenario_short,scenario_short),[Q.base_LCP_pct Q.LCP_pct],'grouped');
    B(1).DisplayName='No-RES baseline'; B(2).DisplayName='Renewable case';
    for i=1:height(Q)
        errorbar(ax,i+0.15,Q.LCP_pct(i),Q.LCP_pct(i)-Q.LCP_ci_low_pct(i), ...
            Q.LCP_ci_high_pct(i)-Q.LCP_pct(i),'k','LineStyle','none', ...
            'CapSize',4,'HandleVisibility','off');
    end
    ylabel(ax,'LCP(80%) [%]'); ylim(ax,ylims); grid(ax,'on');
    title(ax,sprintf('(%c) Branch %s',char('a'+p-1),brname));
end
lgd=legend([B(1) B(2)],{'No-RES baseline','Renewable case'}, ...
    'Orientation','horizontal');
lgd.Layout.Tile='south';
export_figure(f,outdir,'Fig_Structural_vs_RES');

%% Figure 3: Losses and curtailment
f = figure('Units','centimeters','Position',[2 2 17.5 8.2]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl,1); hold(ax,'on');
Q=Tsys(Tsys.N_requested==2000 & string(Tsys.metric)=="daily_loss_MWh",:);
Q=order_scenarios(Q,scenario_order); x=1:4;
errorbar(ax,x,Q.estimate,Q.estimate-Q.ci_low,Q.ci_high-Q.estimate, ...
    'o','MarkerSize',6,'LineStyle','none','CapSize',6);
set(ax,'XTick',x,'XTickLabel',scenario_short);
ylabel(ax,'Daily losses [MWh/day]'); grid(ax,'on');
title(ax,'(a) Intact-network losses');
low_penalty = Q.estimate(1)-Q.estimate(3);
high_penalty = Q.estimate(2)-Q.estimate(4);
text(ax,0.05,0.94, ...
    {sprintf('S1-S3 = %.1f MWh/day',low_penalty), ...
     sprintf('S2-S4 = %.1f MWh/day',high_penalty)}, ...
    'Units','normalized','VerticalAlignment','top', ...
    'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5], ...
    'Margin',5,'FontSize',8);
ax=nexttile(tl,2); hold(ax,'on');
Qw=Tsys(Tsys.N_requested==2000 & string(Tsys.metric)=="wind_curt_MWh",:);
Qp=Tsys(Tsys.N_requested==2000 & string(Tsys.metric)=="trajectory_curt_prob_pct",:);
Qw=order_scenarios(Qw,scenario_order); Qp=order_scenarios(Qp,scenario_order);
yyaxis(ax,'left'); bw=bar(ax,x,Qw.estimate); ylabel(ax,'Wind curtailment [MWh/day]');
yyaxis(ax,'right'); hp=plot(ax,x,Qp.estimate,'-o','MarkerSize',5); ylabel(ax,'Curtailed trajectories [%]');
set(ax,'XTick',x,'XTickLabel',scenario_short); grid(ax,'on'); title(ax,'(b) Curtailment');
legend(ax,[bw hp],{'Mean wind curtailment','Trajectories curtailed'}, ...
    'Location','northoutside','Orientation','horizontal');
export_figure(f,outdir,'Fig_Losses_Curtailment');

%% Figure 4: Deterministic N-1 screening
f=figure('Units','centimeters','Position',[2 2 17.5 8.5]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl,1);
Q=Tshort(1:min(8,height(Tshort)),:); Q=sortrows(Q,'screening_severity_score','descend'); Q=flipud(Q);
barh(ax,categorical(string(Q.outage_branch),string(Q.outage_branch)),Q.screening_severity_score);
xlabel(ax,'Screening severity score'); grid(ax,'on'); title(ax,'(a) Selected connected outages');
xlim(ax,[0 1.15*max(Q.screening_severity_score)]);
for i=1:height(Q)
    text(ax,Q.screening_severity_score(i)+0.02*max(Q.screening_severity_score),i, ...
        sprintf('%.1f',Q.screening_severity_score(i)),'VerticalAlignment','middle','FontSize',8);
end
ax=nexttile(tl,2); hold(ax,'on');
Q=Tdet(Tdet.is_islanding==0 & Tdet.success_hours==24,:);
scatter(ax,Q.mean_loss_delta_pct,100*Q.max_abs_voltage_dev_from_intact_pu,20,Q.screening_severity_score,'filled');
xlabel(ax,'Mean loss change [%]'); ylabel(ax,'Maximum voltage deviation [10^{-2} p.u.]'); grid(ax,'on');
cb=colorbar(ax); cb.Label.String='Severity score';
labels=["26-30","8-5","26-25","23-25","94-100","94-96"];
offsets=[0.8 0.0; -4.0 0.2; 0.8 0.2; 0.8 -0.2; 0.8 0.2; 0.8 -0.2];
for i=1:numel(labels)
    k=find(string(Q.outage_branch)==labels(i),1);
    if ~isempty(k)
        text(ax,Q.mean_loss_delta_pct(k)+offsets(i,1),100*Q.max_abs_voltage_dev_from_intact_pu(k)+offsets(i,2), ...
            labels(i),'FontSize',8,'Clipping','off');
    end
end
xpad=0.08*range(Q.mean_loss_delta_pct); yv=100*Q.max_abs_voltage_dev_from_intact_pu;
xlim(ax,[min(Q.mean_loss_delta_pct)-xpad max(Q.mean_loss_delta_pct)+2*xpad]);
ylim(ax,[0 max(yv)+0.08*range(yv)]); title(ax,'(b) Physical screening dimensions');
export_figure(f,outdir,'Fig_N1_Deterministic');

%% Figure 5: Targeted probabilistic N-1
f=figure('Units','centimeters','Position',[2 2 17.5 8.5]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl,1); hold(ax,'on');
Q=Tn1target(string(Tn1target.scenario)=="S2_High_Conc" & string(Tn1target.target_branch)=="26-30" & Tn1target.target_is_outaged==0,:);
Q=sortrows(Q,'contingency_id'); x=1:height(Q);
B=bar(ax,x,[Q.base_LCP80_pct Q.LCP80_pct],'grouped');
for i=1:height(Q)
    errorbar(ax,i+0.15,Q.LCP80_pct(i),Q.LCP80_pct(i)-Q.ci_low_pct(i),Q.ci_high_pct(i)-Q.LCP80_pct(i), ...
        'k','LineStyle','none','CapSize',4,'HandleVisibility','off');
end
set(ax,'XTick',x,'XTickLabel',string(Q.outage_branch),'XTickLabelRotation',35);
ylabel(ax,'Branch 26--30 LCP(80%) [%]'); ylim(ax,[0 105]); grid(ax,'on');
legend(ax,[B(1) B(2)],{'Contingency baseline','S2 renewable'},'Location','northoutside','Orientation','horizontal');
title(ax,'(a) Absolute post-contingency LCP');
ax=nexttile(tl,2); hold(ax,'on');
delta=Q.LCP80_pct-Q.base_LCP80_pct;
bar(ax,x,delta); yline(ax,0,'--','HandleVisibility','off');
set(ax,'XTick',x,'XTickLabel',string(Q.outage_branch),'XTickLabelRotation',35);
ylabel(ax,'$\Delta \mathrm{LCP}(80\%)$ [percentage points]','Interpreter','latex'); grid(ax,'on'); title(ax,'(b) Renewable increment under S2');
for i=1:height(Q)
    if delta(i)>=0, va='bottom'; dy=0.4; else, va='top'; dy=-0.4; end
    text(ax,i,delta(i)+dy,sprintf('%.2f',delta(i)),'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',8);
end
export_figure(f,outdir,'Fig_N1_Probabilistic');

%% Figure 6: Peak-demand sensitivity
f=figure('Units','centimeters','Position',[2 2 17.5 8.2]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact'); markers={'o','s','^','d'};
for p=1:2
    if p==1, br=155; brname='94--100'; ylims=[0 105]; else, br=38; brname='26--30'; ylims=[0 35]; end
    ax=nexttile(tl,p); hold(ax,'on'); Q=Tpeak(Tpeak.branch_idx==br,:);
    baseQ=Q(string(Q.scenario)==scenario_order(1),:); baseQ=sortrows(baseQ,'peak_frac');
    plot(ax,baseQ.peak_frac,baseQ.base_LCP80_pct,'--','DisplayName','No-RES baseline','LineWidth',1.5);
    for i=1:numel(scenario_order)
        R=Q(string(Q.scenario)==scenario_order(i),:); R=sortrows(R,'peak_frac');
        plot(ax,R.peak_frac,R.res_LCP80_pct,'-','Marker',markers{i},'MarkerSize',4,'DisplayName',scenario_short(i));
    end
    xlabel(ax,'Peak-demand factor'); ylabel(ax,'LCP(80%) [%]'); ylim(ax,ylims); grid(ax,'on');
    title(ax,sprintf('(%c) Branch %s',char('a'+p-1),brname));
    if p==1, text(ax,0.03,0.08,'100% for all cases','Units','normalized','FontSize',8); else, legend(ax,'Location','northeast'); end
end
export_figure(f,outdir,'Fig_Peak_Sensitivity');

%% Supplementary figure: ranking stability
f=figure('Units','centimeters','Position',[2 2 12 8]);
ax=axes(f); hold(ax,'on'); markers={'o','s','^','d'}; offsets=[-18 -6 6 18];
for i=1:numel(scenario_order)
    Q=Trank(string(Trank.scenario)==scenario_order(i),:); Q=sortrows(Q,'N_requested');
    plot(ax,Q.N_requested+offsets(i),Q.top10_jaccard_pct,'-','Marker',markers{i},'MarkerSize',4,'DisplayName',scenario_short(i));
end
xlabel(ax,'Number of trajectories, N'); ylabel(ax,'Top-set overlap with N=2000 [%]');
ylim(ax,[0 105]); xlim(ax,[0 2050]); grid(ax,'on'); legend(ax,'Location','best'); title(ax,'Branch-ranking stability');
export_figure(f,outdir,'FigS_Rank_Stability');

fprintf('\nFinal figures written to:\n  %s\n',outdir);

function path = required_file(folder,name)
    path=fullfile(folder,name);
    if ~isfile(path), error('Required CSV not found: %s',path); end
end

function Q=order_scenarios(Q,scenario_order)
    [tf,ord]=ismember(string(Q.scenario),scenario_order);
    Q=Q(tf,:); ord=ord(tf);
    [~,idx]=sort(ord); Q=Q(idx,:);
end

function export_figure(f,outdir,stem)
    drawnow;
    exportgraphics(f,fullfile(outdir,stem+".pdf"),'ContentType','vector');
    exportgraphics(f,fullfile(outdir,stem+".png"),'Resolution',600);
    savefig(f,fullfile(outdir,stem+".fig"));
end
