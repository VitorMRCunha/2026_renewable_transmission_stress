%% =========================================================================
%% peak_demand_sensitivity_multisite_v1.m
%% Peak-demand scaling sensitivity — IEEE 118-bus, multi-site weather
%%
%% Purpose
%%   Test whether structural and renewable-sensitive branch classifications
%%   remain robust when the system peak-demand scaling factor is varied.
%%
%% Design
%%   f_peak = [0.60 0.65 0.70 0.75 0.80]
%%   4 factorial RES scenarios, Ns = 400 common trajectories, 24 hours
%%   trajectory-level bootstrap confidence intervals
%%
%% Required files
%%   case_pglib_opf_case118_ieee.m
%%   multisite_weather/multisite_weather_model.mat
%%   multisite_weather/solar_temporal_corr.csv
%%
%% Main outputs: results_paper/peak_sensitivity/
%%   peak_sensitivity_run_summary.csv
%%   peak_sensitivity_opf_summary.csv
%%   peak_sensitivity_baseline_all_branches.csv
%%   peak_sensitivity_target_branches.csv
%%   peak_sensitivity_system_metrics.csv
%%   peak_sensitivity_branch_robustness.csv
%%
%% Authors: Azevedo & Cunha, ISEP/IPP
%% =========================================================================

clear; clc; close all;
rng(42,'twister');

%% AUDIT PATCH: branch loading uses max apparent power at both ends
%% Sf = hypot(PF,QF); St = hypot(PT,QT); S = max(Sf,St).
%% Output folder is separate from the earlier from-end-only run.

%% -------------------------------------------------------------------------
%% 0. Configuration
%% -------------------------------------------------------------------------
T = 24;
Ns = 400;
B = 500;
ci_level = 0.95;
bootstrap_seed = 4201;

peak_grid = [0.60 0.65 0.70 0.75 0.80];
reference_peak = 0.70;

Tau = [0.70 0.80 0.90];
Vscreen_min = 0.95;
Vscreen_max = 1.05;
curt_tol_MW = 1e-4;

use_parallel = true;
requested_workers = 6;
save_compact_mat = true;

% Explicit screening flags used only for the robustness table. They are not
% claimed as universal definitions.
structural_LCP80_threshold_pct = 50;
res_sensitive_DeltaLCP80_threshold_pp = 2;

d_t = [0.67,0.63,0.60,0.59,0.59,0.60,0.64,0.71,0.78,0.84, ...
       0.88,0.90,0.89,0.88,0.87,0.87,0.88,0.91,0.94,0.97, ...
       1.00,0.96,0.87,0.75];

weather_dir = 'multisite_weather';
weather_mat = fullfile(weather_dir,'multisite_weather_model.mat');
solar_temporal_csv = fullfile(weather_dir,'solar_temporal_corr.csv');
case_file = './case_pglib_opf_case118_ieee.m';

if ~exist(weather_mat,'file'), error('Missing %s.',weather_mat); end
if ~exist(solar_temporal_csv,'file'), error('Missing %s.',solar_temporal_csv); end
if ~exist(case_file,'file'), error('Missing %s.',case_file); end

outdir = fullfile('results_paper','peak_sensitivity');
if ~exist(outdir,'dir'), mkdir(outdir); end

%% -------------------------------------------------------------------------
%% 1. Load weather calibration and generate one common master sample
%% -------------------------------------------------------------------------
fprintf('=== Peak-demand sensitivity, multi-site weather ===\n');
fprintf('f_peak values: %s\n',mat2str(peak_grid));
fprintf('Scenarios: 4 | trajectories: %d | hours: %d\n',Ns,T);
fprintf('Requested stochastic OPFs: %d\n',numel(peak_grid)*4*Ns*T);

W = load(weather_mat);
required_fields = {'solar_alpha_hourly','solar_beta_hourly','solar_envelope_norm', ...
    'solar_corr_gaussian','wind_k_weibull','wind_c_weibull_ms', ...
    'wind_rho_ar1','wind_corr_gaussian'};
for ii=1:numel(required_fields)
    if ~isfield(W,required_fields{ii})
        error('Field %s is absent from %s.',required_fields{ii},weather_mat);
    end
end

nzone = size(W.solar_corr_gaussian,1);
if nzone ~= 6, error('Expected six meteorological zones, found %d.',nzone); end

Tsol = readtable(solar_temporal_csv,'TextType','string');
solar_rho_ar1 = NaN(nzone,1);
for iz=1:nzone
    zid = sprintf('Z%d_',iz);
    mask = startsWith(Tsol.zone,zid) & Tsol.lag_hours==1;
    if nnz(mask)~=1
        error('Could not identify one lag-1 solar correlation for zone %d.',iz);
    end
    solar_rho_ar1(iz)=Tsol.correlation(mask);
end

fprintf('\n[1] Generating common spatially and temporally correlated trajectories...\n');
Xsolar = generate_vector_ar1(Ns,T,W.solar_corr_gaussian,solar_rho_ar1);
Xwind  = generate_vector_ar1(Ns,T,W.wind_corr_gaussian,W.wind_rho_ar1(:));
Usolar = min(max(normcdf(Xsolar),1e-8),1-1e-8);
Uwind  = min(max(normcdf(Xwind), 1e-8),1-1e-8);

v_ci=3.0; v_r=12.0; v_co=25.0;
G_solar_zone=zeros(Ns,T,nzone);
P_wind_zone=zeros(Ns,T,nzone);
for iz=1:nzone
    for t=1:T
        a=W.solar_alpha_hourly(iz,t);
        b=W.solar_beta_hourly(iz,t);
        env=W.solar_envelope_norm(iz,t);
        if env<=0 || ~isfinite(a) || ~isfinite(b) || a<=0 || b<=0
            G_solar_zone(:,t,iz)=0;
        else
            G_solar_zone(:,t,iz)=betainv(Usolar(:,t,iz),a,b)*env;
        end
    end
    v=wblinv(Uwind(:,:,iz),W.wind_c_weibull_ms(iz),W.wind_k_weibull(iz));
    P_wind_zone(:,:,iz)=wind_pc(v,v_ci,v_r,v_co,1.0);
end

%% -------------------------------------------------------------------------
%% 2. Network, scenarios and parallel pool
%% -------------------------------------------------------------------------
fprintf('[2] Loading network and constructing scenarios...\n');
mpc0=loadcase(case_file);
nb=size(mpc0.bus,1); nbr=size(mpc0.branch,1);
Ppk=sum(mpc0.bus(:,3));
Pd_frac=mpc0.bus(:,3)/Ppk;
Qd0=mpc0.bus(:,4);
Prat=mpc0.branch(:,6);
if any(~isfinite(Prat) | Prat<=0), error('All RATE_A values must be positive.'); end

P_solar_low=200; P_wind_low=300;
P_solar_high=400; P_wind_high=600;

CONC.solar_buses=45; CONC.wind_buses=75;
CONC.solar_frac=1.0; CONC.wind_frac=1.0;
CONC.solar_zones=4; CONC.wind_zones=3;

DIST.solar_buses=[9 45]; DIST.wind_buses=[50 75];
DIST.solar_frac=[0.4 0.6]; DIST.wind_frac=[0.35 0.65];
DIST.solar_zones=[2 5]; DIST.wind_zones=[1 6];

sc(1)=make_scenario('S1_Low_Conc','Low/Concentrated',P_solar_low,P_wind_low,CONC);
sc(2)=make_scenario('S2_High_Conc','High/Concentrated',P_solar_high,P_wind_high,CONC);
sc(3)=make_scenario('S3_Low_Dist','Low/Distributed',P_solar_low,P_wind_low,DIST);
sc(4)=make_scenario('S4_High_Dist','High/Distributed',P_solar_high,P_wind_high,DIST);

% Key branches requested by the reviewer narrative and prior results.
target_pairs=[94 100;26 30;23 25;94 96];
target_names=["94-100";"26-30";"23-25";"94-96"];
target_idx=zeros(size(target_pairs,1),1);
for k=1:size(target_pairs,1)
    target_idx(k)=find_branch(mpc0,target_pairs(k,1),target_pairs(k,2));
end

mpopt=mpoption('verbose',0,'out.all',0);

% Reproducibility metadata
try
    matlab_release = version('-release');
catch
    matlab_release = version;
end
try
    [mpv, ~] = mpver;
    matpower_version = string(mpv);
catch
    try
        matpower_version = string(strtrim(evalc('mpver')));
    catch
        matpower_version = "unknown";
    end
end
try
    resolved_ac_solver = string(mpopt.opf.ac.solver);
    if strlength(resolved_ac_solver)==0
        resolved_ac_solver = "MATPOWER default/automatic";
    end
catch
    resolved_ac_solver = "MATPOWER default/automatic";
end
try
    pool_now = gcp('nocreate');
    if isempty(pool_now), parallel_workers_at_setup = 0; else, parallel_workers_at_setup = pool_now.NumWorkers; end
catch
    parallel_workers_at_setup = NaN;
end
Tsoftware = table(string(matlab_release),string(matpower_version),string(resolved_ac_solver), ...
    parallel_workers_at_setup,string(computer('arch')), ...
    'VariableNames',{'MATLABRelease','MATPOWERVersion','ResolvedACOPFSolver', ...
    'ParallelWorkersAtSetup','ComputerArchitecture'});
writetable(Tsoftware,fullfile(outdir,'software_environment.csv'));

if use_parallel && license('test','Distrib_Computing_Toolbox')
    pool=gcp('nocreate');
    if isempty(pool)
        parpool('local',requested_workers);
    elseif pool.NumWorkers~=requested_workers
        fprintf('[POOL] Existing pool has %d workers; retaining it.\n',pool.NumWorkers);
    end
    mp_dir=fileparts(which('runopf'));
    if ~isempty(mp_dir)
        pctRunOnAll(sprintf('addpath(genpath(''%s''))',mp_dir));
    end
else
    use_parallel=false;
end

% Warm-up OPF outside all timing records.
mtw=mpc0;
mtw.bus(:,3)=Pd_frac*(Ppk*reference_peak*d_t(1));
mtw.bus(:,4)=Qd0*(reference_peak*d_t(1));
runopf(mtw,mpopt);

%% -------------------------------------------------------------------------
%% 3. Deterministic no-RES baseline at each peak-demand scaling
%% -------------------------------------------------------------------------
fprintf('\n[3] Running deterministic no-RES peak sensitivity...\n');
np=numel(peak_grid); nsc=numel(sc);
base=repmat(struct(),np,1);
base_rows={};

for ip=1:np
    fp=peak_grid(ip);
    bout=run_baseline_peak(mpc0,fp,Pd_frac,Qd0,Ppk,d_t,Prat,Tau, ...
        Vscreen_min,Vscreen_max,mpopt);
    if ip==1
        base=repmat(bout,np,1);
    else
        base(ip)=bout;
    end

    for br=1:nbr
        base_rows(end+1,:)={fp,br,mpc0.branch(br,1),mpc0.branch(br,2),Prat(br), ...
            bout.LCP(br,1),bout.LCP(br,2),bout.LCP(br,3), ...
            bout.max_loading_pct(br)}; %#ok<SAGROW>
    end
    fprintf('  f_peak=%.2f: success %.1f%% | mean cost %.2f | mean loss %.2f MW\n', ...
        fp,100*mean(bout.success),mean(bout.cost,'omitnan'),mean(bout.loss,'omitnan'));
end

%% -------------------------------------------------------------------------
%% 4. Stochastic multi-site factorial sensitivity
%% -------------------------------------------------------------------------
fprintf('\n[4] Running stochastic sensitivity experiments...\n');
res=repmat(struct(),np,nsc);

for ip=1:np
    fp=peak_grid(ip);
    fprintf('\n  Peak scaling %.2f (%d/%d)\n',fp,ip,np);

    for is=1:nsc
        cfg=sc(is);
        fprintf('    Scenario %d/4: %s\n',is,cfg.label);
        [mpc_sc,idx_s,idx_w]=add_res_generators(mpc0,cfg);
        [Psolar_units,Pwind_units]=scenario_available_power(cfg,G_solar_zone,P_wind_zone);

        exceed_count=zeros(Ns,nbr,numel(Tau),'uint8');
        daily_cost=NaN(Ns,1);
        daily_loss=NaN(Ns,1);
        daily_solar_curt=NaN(Ns,1);
        daily_wind_curt=NaN(Ns,1);
        daily_vmin=NaN(Ns,1);
        daily_vmax=NaN(Ns,1);
        daily_vviol_hours=NaN(Ns,1);
        successful_hours=zeros(Ns,1,'uint8');
        total_solve_time=NaN(Ns,1);

        tcase=tic;
        if use_parallel
            parfor s=1:Ns
                day=run_day_peak(mpc_sc,squeeze(Psolar_units(s,:,:)), ...
                    squeeze(Pwind_units(s,:,:)),idx_s,idx_w,fp, ...
                    Pd_frac,Qd0,Ppk,d_t,Prat,Tau,Vscreen_min,Vscreen_max, ...
                    curt_tol_MW,mpopt);
                exceed_count(s,:,:)=day.exceed_count;
                daily_cost(s)=day.daily_cost;
                daily_loss(s)=day.daily_loss_MWh;
                daily_solar_curt(s)=day.daily_solar_curt_MWh;
                daily_wind_curt(s)=day.daily_wind_curt_MWh;
                daily_vmin(s)=day.daily_vmin;
                daily_vmax(s)=day.daily_vmax;
                daily_vviol_hours(s)=day.daily_vviol_hours;
                successful_hours(s)=day.successful_hours;
                total_solve_time(s)=day.total_solve_time_s;
            end
        else
            for s=1:Ns
                day=run_day_peak(mpc_sc,squeeze(Psolar_units(s,:,:)), ...
                    squeeze(Pwind_units(s,:,:)),idx_s,idx_w,fp, ...
                    Pd_frac,Qd0,Ppk,d_t,Prat,Tau,Vscreen_min,Vscreen_max, ...
                    curt_tol_MW,mpopt);
                exceed_count(s,:,:)=day.exceed_count;
                daily_cost(s)=day.daily_cost;
                daily_loss(s)=day.daily_loss_MWh;
                daily_solar_curt(s)=day.daily_solar_curt_MWh;
                daily_wind_curt(s)=day.daily_wind_curt_MWh;
                daily_vmin(s)=day.daily_vmin;
                daily_vmax(s)=day.daily_vmax;
                daily_vviol_hours(s)=day.daily_vviol_hours;
                successful_hours(s)=day.successful_hours;
                total_solve_time(s)=day.total_solve_time_s;
            end
        end

        res(ip,is).peak_frac=fp;
        res(ip,is).name=cfg.name;
        res(ip,is).label=cfg.label;
        res(ip,is).exceed_count=exceed_count;
        res(ip,is).daily_cost=daily_cost;
        res(ip,is).daily_loss_MWh=daily_loss;
        res(ip,is).daily_solar_curt_MWh=daily_solar_curt;
        res(ip,is).daily_wind_curt_MWh=daily_wind_curt;
        res(ip,is).daily_vmin=daily_vmin;
        res(ip,is).daily_vmax=daily_vmax;
        res(ip,is).daily_vviol_hours=daily_vviol_hours;
        res(ip,is).successful_hours=successful_hours;
        res(ip,is).total_solve_time_s=total_solve_time;
        res(ip,is).wall_time_s=toc(tcase);

        fprintf('      completed in %.1f min | complete trajectories %.1f%%\n', ...
            res(ip,is).wall_time_s/60,100*mean(double(successful_hours)==T));
    end
end

%% -------------------------------------------------------------------------
%% 5. Statistical summaries and CSV rows
%% -------------------------------------------------------------------------
fprintf('\n[5] Calculating bootstrap intervals and robustness metrics...\n');
rng(bootstrap_seed,'twister');

opf_rows={}; target_rows={}; system_rows={};

% Baseline OPF summary rows
for ip=1:np
    fp=peak_grid(ip);
    opf_rows(end+1,:)={fp,"BASE_NoRES","No RES",T,nnz(base(ip).success), ...
        nnz(~base(ip).success),100*mean(base(ip).success), ...
        mean(base(ip).solve_time,'omitnan'),median(base(ip).solve_time,'omitnan'), ...
        max(base(ip).solve_time,[],'omitnan'),sum(base(ip).solve_time,'omitnan')}; %#ok<SAGROW>
end

for ip=1:np
    fp=peak_grid(ip);
    for is=1:nsc
        valid=double(res(ip,is).successful_hours)==T;
        nvalid=nnz(valid);
        requested=Ns*T;
        successful=sum(double(res(ip,is).successful_hours));
        failed=requested-successful;
        per_opf_time=res(ip,is).total_solve_time_s/T;
        opf_rows(end+1,:)={fp,string(res(ip,is).name),string(res(ip,is).label), ...
            requested,successful,failed,100*successful/requested, ...
            mean(per_opf_time,'omitnan'),median(per_opf_time,'omitnan'), ...
            max(per_opf_time,[],'omitnan'),sum(res(ip,is).total_solve_time_s,'omitnan')}; %#ok<SAGROW>

        metric_defs={ ...
            'daily_cost','currency/day'; ...
            'daily_loss_MWh','MWh/day'; ...
            'daily_solar_curt_MWh','MWh/day'; ...
            'daily_wind_curt_MWh','MWh/day'; ...
            'daily_vmin','p.u.'; ...
            'daily_vmax','p.u.'; ...
            'daily_vviol_hours','bus-hours/day'};
        for im=1:size(metric_defs,1)
            field=metric_defs{im,1}; unit=metric_defs{im,2};
            x=res(ip,is).(field)(valid);
            [est,se,lo,hi]=bootstrap_mean(x,B,ci_level);
            system_rows(end+1,:)={fp,string(res(ip,is).name),string(res(ip,is).label), ...
                nvalid,string(field),est,se,lo,hi,string(unit)}; %#ok<SAGROW>
        end

        % Curtailment occurrence probability by complete trajectory.
        x=(res(ip,is).daily_solar_curt_MWh(valid)+res(ip,is).daily_wind_curt_MWh(valid))>curt_tol_MW;
        [est,se,lo,hi]=bootstrap_mean(double(x),B,ci_level);
        system_rows(end+1,:)={fp,string(res(ip,is).name),string(res(ip,is).label), ...
            nvalid,"curtailment_trajectory_probability",100*est,100*se,100*lo,100*hi,"%"}; %#ok<SAGROW>

        for kt=1:numel(target_idx)
            br=target_idx(kt);
            counts=double(res(ip,is).exceed_count(valid,br,2));
            lcp=mean(counts)/T*100;
            [se,lo,hi]=bootstrap_lcp_vector(counts,T,B,ci_level);
            base_lcp=base(ip).LCP(br,2);
            target_rows(end+1,:)={fp,string(res(ip,is).name),string(res(ip,is).label), ...
                nvalid,br,target_pairs(kt,1),target_pairs(kt,2),target_names(kt),Prat(br), ...
                base_lcp,lcp,se,lo,hi,lcp-base_lcp,lo-base_lcp,hi-base_lcp}; %#ok<SAGROW>
        end
    end
end

Topf=cell2table(opf_rows,'VariableNames', ...
    {'peak_frac','scenario','scenario_label','requested_opfs','successful_opfs', ...
     'failed_opfs','success_rate_pct','mean_time_per_opf_s','median_time_per_opf_s', ...
     'max_mean_daily_opf_time_s','total_solver_time_s'});
writetable(Topf,fullfile(outdir,'peak_sensitivity_opf_summary.csv'));

Tbase=cell2table(base_rows,'VariableNames', ...
    {'peak_frac','branch_idx','from_bus','to_bus','rating_MVA', ...
     'base_LCP70_pct','base_LCP80_pct','base_LCP90_pct','base_max_loading_pct'});
writetable(Tbase,fullfile(outdir,'peak_sensitivity_baseline_all_branches.csv'));

Ttarget=cell2table(target_rows,'VariableNames', ...
    {'peak_frac','scenario','scenario_label','N_valid','branch_idx','from_bus','to_bus', ...
     'branch_name','rating_MVA','base_LCP80_pct','res_LCP80_pct','bootstrap_se_pp', ...
     'res_ci_low_pct','res_ci_high_pct','DeltaLCP80_pp','Delta_ci_low_pp','Delta_ci_high_pp'});
writetable(Ttarget,fullfile(outdir,'peak_sensitivity_target_branches.csv'));

Tsystem=cell2table(system_rows,'VariableNames', ...
    {'peak_frac','scenario','scenario_label','N_valid','metric','estimate','bootstrap_se', ...
     'ci_low','ci_high','unit'});
writetable(Tsystem,fullfile(outdir,'peak_sensitivity_system_metrics.csv'));

%% Branch robustness across f_peak values
robust_rows={};
for br=1:nbr
    base80=zeros(np,1);
    maxdelta=zeros(np,1);
    mindelta=zeros(np,1);
    for ip=1:np
        base80(ip)=base(ip).LCP(br,2);
        ds=zeros(nsc,1);
        for is=1:nsc
            valid=double(res(ip,is).successful_hours)==T;
            lcp=mean(double(res(ip,is).exceed_count(valid,br,2)))/T*100;
            ds(is)=lcp-base80(ip);
        end
        maxdelta(ip)=max(ds);
        mindelta(ip)=min(ds);
    end
    structural_all=all(base80>=structural_LCP80_threshold_pct);
    structural_any=any(base80>=structural_LCP80_threshold_pct);
    res_sensitive_any=any(maxdelta>=res_sensitive_DeltaLCP80_threshold_pp);
    robust_rows(end+1,:)={br,mpc0.branch(br,1),mpc0.branch(br,2),Prat(br), ...
        min(base80),max(base80),base80(peak_grid==reference_peak), ...
        min(mindelta),max(maxdelta),structural_all,structural_any,res_sensitive_any}; %#ok<SAGROW>
end

Trobust=cell2table(robust_rows,'VariableNames', ...
    {'branch_idx','from_bus','to_bus','rating_MVA','min_base_LCP80_pct', ...
     'max_base_LCP80_pct','base_LCP80_at_0p70_pct','min_DeltaLCP80_pp', ...
     'max_DeltaLCP80_pp','structural_flag_all_peak_levels', ...
     'structural_flag_any_peak_level','res_sensitive_flag_any_peak_level'});
writetable(Trobust,fullfile(outdir,'peak_sensitivity_branch_robustness.csv'));

%% Run-level summary
elapsed_solver=sum(Topf.total_solver_time_s,'omitnan');
run_summary=table(Ns,T,B,numel(peak_grid),nsc, ...
    numel(peak_grid)*nsc*Ns*T, sum(Topf.successful_opfs),sum(Topf.failed_opfs), ...
    elapsed_solver,structural_LCP80_threshold_pct,res_sensitive_DeltaLCP80_threshold_pp, ...
    'VariableNames',{'Ns','hours','bootstrap_replicates','n_peak_levels','n_scenarios', ...
    'requested_stochastic_opfs','successful_all_rows_including_baselines', ...
    'failed_all_rows_including_baselines','total_accumulated_solver_time_s', ...
    'structural_screen_threshold_LCP80_pct','res_sensitive_screen_threshold_DeltaLCP80_pp'});
writetable(run_summary,fullfile(outdir,'peak_sensitivity_run_summary.csv'));

if save_compact_mat
    save(fullfile(outdir,'peak_sensitivity_compact.mat'),'res','base','peak_grid', ...
        'reference_peak','sc','target_pairs','target_idx','Tau','Ns','T','B', ...
        'G_solar_zone','P_wind_zone','W','solar_rho_ar1','mpc0','Prat','d_t','-v7.3');
end

fprintf('\n=== Peak-demand sensitivity completed ===\n');
fprintf('Outputs written to: %s\n',outdir);
fprintf('Please share these files first:\n');
fprintf('  peak_sensitivity_run_summary.csv\n');
fprintf('  peak_sensitivity_opf_summary.csv\n');
fprintf('  peak_sensitivity_target_branches.csv\n');
fprintf('  peak_sensitivity_system_metrics.csv\n');
fprintf('  peak_sensitivity_branch_robustness.csv\n');

%% =========================================================================
%% Local functions
%% =========================================================================
function sc=make_scenario(name,label,Ps,Pw,inj)
sc.name=name; sc.label=label; sc.P_solar_max=Ps; sc.P_wind_max=Pw; sc.inj=inj;
end

function P=wind_pc(v,vci,vr,vco,Prated)
P=zeros(size(v)); r1=(v>=vci)&(v<=vr);
P(r1)=Prated*(v(r1).^3-vci^3)/(vr^3-vci^3);
P(v>vr & v<vco)=Prated;
end

function [mpc,idx_s,idx_w]=add_res_generators(mpc0,cfg)
mpc=mpc0; idx_s=zeros(1,numel(cfg.inj.solar_buses)); idx_w=zeros(1,numel(cfg.inj.wind_buses));
for i=1:numel(cfg.inj.solar_buses)
    Pi=cfg.P_solar_max*cfg.inj.solar_frac(i); row=zeros(1,size(mpc.gen,2));
    row(1:10)=[cfg.inj.solar_buses(i),0,0,0.5*Pi,-0.5*Pi,1,100,1,Pi,0];
    mpc.gen=[mpc.gen;row]; mpc.gencost=[mpc.gencost;zero_cost_row(mpc.gencost)];
    idx_s(i)=size(mpc.gen,1);
end
for i=1:numel(cfg.inj.wind_buses)
    Pi=cfg.P_wind_max*cfg.inj.wind_frac(i); row=zeros(1,size(mpc.gen,2));
    row(1:10)=[cfg.inj.wind_buses(i),0,0,0.5*Pi,-0.5*Pi,1,100,1,Pi,0];
    mpc.gen=[mpc.gen;row]; mpc.gencost=[mpc.gencost;zero_cost_row(mpc.gencost)];
    idx_w(i)=size(mpc.gen,1);
end
end

function gc=zero_cost_row(gencost)
ncol=size(gencost,2); gc=zeros(1,ncol); gc(1)=2; gc(4)=ncol-4;
end

function [Psolar_units,Pwind_units]=scenario_available_power(cfg,Gzone,Wzone)
Ns=size(Gzone,1); T=size(Gzone,2); ns=numel(cfg.inj.solar_buses); nw=numel(cfg.inj.wind_buses);
Psolar_units=zeros(Ns,T,ns); Pwind_units=zeros(Ns,T,nw);
for i=1:ns
    Psolar_units(:,:,i)=cfg.P_solar_max*cfg.inj.solar_frac(i)*Gzone(:,:,cfg.inj.solar_zones(i));
end
for i=1:nw
    Pwind_units(:,:,i)=cfg.P_wind_max*cfg.inj.wind_frac(i)*Wzone(:,:,cfg.inj.wind_zones(i));
end
end

function X=generate_vector_ar1(Ns,T,R,rho)
R=nearest_corr((R+R')/2,1e-10); rho=rho(:); A=diag(rho);
Q=project_psd(R-A*R*A',1e-10); L0=chol(R,'lower'); Lq=chol(Q,'lower');
X=zeros(Ns,T,numel(rho)); X(:,1,:)=randn(Ns,numel(rho))*L0';
for t=2:T
    prev=squeeze(X(:,t-1,:)); eps=randn(Ns,numel(rho))*Lq';
    X(:,t,:)=prev*A'+eps;
end
end

function M=project_psd(M,tol)
M=(M+M')/2; [V,D]=eig(M); d=diag(D); d(d<tol)=tol;
M=V*diag(d)*V'; M=(M+M')/2;
end

function R=nearest_corr(R,tol)
R=project_psd(R,tol); s=sqrt(diag(R)); R=R./(s*s');
R=project_psd((R+R')/2,tol); s=sqrt(diag(R)); R=R./(s*s'); R=(R+R')/2;
end

function bout=run_baseline_peak(mpc,fp,Pd_frac,Qd0,Ppk,d_t,Prat,Tau,Vmin,Vmax,mpopt)
T=numel(d_t); nbr=size(mpc.branch,1); nb=size(mpc.bus,1);
bout.success=false(T,1); bout.solve_time=NaN(T,1); bout.cost=NaN(T,1); bout.loss=NaN(T,1);
bout.S=NaN(T,nbr); bout.V=NaN(T,nb);
for t=1:T
    mt=mpc; scale=fp*d_t(t); mt.bus(:,3)=Pd_frac*(Ppk*scale); mt.bus(:,4)=Qd0*scale;
    ts=tic; r=runopf(mt,mpopt); bout.solve_time(t)=toc(ts); bout.success(t)=logical(r.success);
    if r.success
        bout.cost(t)=r.f; bout.loss(t)=sum(r.branch(:,14)+r.branch(:,16));
        bout.S(t,:)=max(hypot(r.branch(:,14),r.branch(:,15)), hypot(r.branch(:,16),r.branch(:,17)))'; bout.V(t,:)=r.bus(:,8)';
    end
end
bout.LCP=NaN(nbr,numel(Tau));
for j=1:numel(Tau), bout.LCP(:,j)=mean(bout.S>Tau(j)*Prat',1,'omitnan')'*100; end
bout.max_loading_pct=max(bout.S./Prat',[],1,'omitnan')'*100;
bout.min_voltage=min(bout.V,[],'all','omitnan'); bout.max_voltage=max(bout.V,[],'all','omitnan');
bout.vviol_bus_hours=nnz(bout.V<Vmin | bout.V>Vmax);
end

function day=run_day_peak(mpc,Ps_units,Pw_units,idx_s,idx_w,fp, ...
    Pd_frac,Qd0,Ppk,d_t,Prat,Tau,Vmin,Vmax,curt_tol,mpopt)
T=numel(d_t); if isvector(Ps_units),Ps_units=Ps_units(:);end
if isvector(Pw_units),Pw_units=Pw_units(:);end
if size(Ps_units,1)~=T,Ps_units=Ps_units';end
if size(Pw_units,1)~=T,Pw_units=Pw_units';end
nbr=size(mpc.branch,1); exceed=false(T,nbr,numel(Tau));
costs=NaN(T,1); losses=NaN(T,1); curtS=NaN(T,1); curtW=NaN(T,1);
vmins=NaN(T,1); vmaxs=NaN(T,1); vviol=NaN(T,1); ok=false(T,1); stime=NaN(T,1);
for t=1:T
    mt=mpc; scale=fp*d_t(t); mt.bus(:,3)=Pd_frac*(Ppk*scale); mt.bus(:,4)=Qd0*scale;
    for i=1:numel(idx_s),mt.gen(idx_s(i),9)=Ps_units(t,i);end
    for i=1:numel(idx_w),mt.gen(idx_w(i),9)=Pw_units(t,i);end
    ts=tic; r=runopf(mt,mpopt); stime(t)=toc(ts); ok(t)=logical(r.success);
    if r.success
        S=max(hypot(r.branch(:,14),r.branch(:,15)), hypot(r.branch(:,16),r.branch(:,17))); V=r.bus(:,8);
        for j=1:numel(Tau),exceed(t,:,j)=(S>Tau(j)*Prat)';end
        costs(t)=r.f; losses(t)=sum(r.branch(:,14)+r.branch(:,16));
        curtS(t)=sum(max(0,Ps_units(t,:)'-r.gen(idx_s,2)));
        curtW(t)=sum(max(0,Pw_units(t,:)'-r.gen(idx_w,2)));
        if curtS(t)<curt_tol,curtS(t)=0;end
        if curtW(t)<curt_tol,curtW(t)=0;end
        vmins(t)=min(V); vmaxs(t)=max(V); vviol(t)=nnz(V<Vmin | V>Vmax);
    end
end
if all(ok)
    day.exceed_count=uint8(squeeze(sum(exceed,1))); day.daily_cost=sum(costs);
    day.daily_loss_MWh=sum(losses); day.daily_solar_curt_MWh=sum(curtS);
    day.daily_wind_curt_MWh=sum(curtW); day.daily_vmin=min(vmins);
    day.daily_vmax=max(vmaxs); day.daily_vviol_hours=sum(vviol);
else
    day.exceed_count=zeros(nbr,numel(Tau),'uint8'); day.daily_cost=NaN;
    day.daily_loss_MWh=NaN; day.daily_solar_curt_MWh=NaN; day.daily_wind_curt_MWh=NaN;
    day.daily_vmin=NaN; day.daily_vmax=NaN; day.daily_vviol_hours=NaN;
end
day.successful_hours=uint8(nnz(ok)); day.total_solve_time_s=sum(stime,'omitnan');
end

function [estimate,se,lo,hi]=bootstrap_mean(x,B,ci_level)
x=x(:); x=x(isfinite(x)); n=numel(x); estimate=mean(x);
if n<2,se=NaN;lo=NaN;hi=NaN;return;end
boot=zeros(B,1); for b=1:B,boot(b)=mean(x(randi(n,n,1)));end
se=std(boot,0); alpha=(1-ci_level)/2; q=prctile(boot,[100*alpha 100*(1-alpha)]); lo=q(1);hi=q(2);
end

function [se,lo,hi]=bootstrap_lcp_vector(counts,T,B,ci_level)
counts=counts(:); n=numel(counts); boot=zeros(B,1);
for b=1:B,boot(b)=mean(counts(randi(n,n,1)))/T*100;end
se=std(boot,0); alpha=(1-ci_level)/2; q=prctile(boot,[100*alpha 100*(1-alpha)]);lo=q(1);hi=q(2);
end

function idx=find_branch(mpc,fbus,tbus)
idx=find((mpc.branch(:,1)==fbus & mpc.branch(:,2)==tbus) | ...
         (mpc.branch(:,1)==tbus & mpc.branch(:,2)==fbus),1);
if isempty(idx),error('Branch %d-%d not found.',fbus,tbus);end
end
