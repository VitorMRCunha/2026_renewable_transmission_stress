%% AUDIT PATCH: branch loading uses max apparent power at both ends
%% Sf = hypot(PF,QF); St = hypot(PT,QT); S = max(Sf,St).
%% Output folder is separate from the earlier from-end-only run.

%% =========================================================================
%% n1_probabilistic_multisite_pilot_v1.m
%% Probabilistic multi-site N-1 AC-OPF pilot — PGLib IEEE 118-bus
%%
%% Purpose
%%   1) Reuse the calibrated six-zone PVGIS/ERA5 weather model.
%%   2) Evaluate four RES factorial scenarios under six critical branch outages.
%%   3) Compare each stochastic contingency case against its own deterministic
%%      no-RES N-1 baseline (not against the intact baseline).
%%   4) Export solver diagnostics, LCP/Delta-LCP, voltage, losses, cost,
%%      curtailment, confidence intervals, and severity rankings.
%%
%% Scope
%%   This is a contingency-conditioned AC-OPF assessment, not a simultaneous
%%   security-constrained AC-OPF. One branch is removed before each hourly OPF.
%%
%% Required
%%   - MATPOWER 7.1+
%%   - Statistics and Machine Learning Toolbox
%%   - Parallel Computing Toolbox (optional)
%%   - case_pglib_opf_case118_ieee.m
%%   - multisite_weather/multisite_weather_model.mat
%%   - multisite_weather/solar_temporal_corr.csv
%%
%% Main outputs
%%   results_paper/probabilistic_n1/
%%       n1_prob_opf_summary.csv
%%       n1_prob_contingency_baselines.csv
%%       n1_prob_target_branch_metrics.csv
%%       n1_prob_system_metrics.csv
%%       n1_prob_scenario_contingency_ranking.csv
%%       n1_prob_paired_scenario_differences.csv
%%       n1_prob_daily_metrics_<scenario>_<outage>.csv
%%
%% Authors: Azevedo & Cunha, ISEP/IPP
%% =========================================================================

clear; clc; close all;
rng(42,'twister');

%% =========================================================================
%% 0. USER CONFIGURATION
%% =========================================================================
T = 24;
Ns = 400;                 % pilot; increase only after reviewing pilot results
B = 500;                  % trajectory bootstrap replicates
ci_level = 0.95;
bootstrap_seed = 7421;

use_parallel = true;
requested_workers = 6;
export_daily_metrics = true;
save_compact_mat = true;

Tau = [0.70 0.80 0.90];
Vscreen_min = 0.95;
Vscreen_max = 1.05;
curt_tol_MW = 1e-4;

peak_frac = 0.70;
d_t = [0.67,0.63,0.60,0.59,0.59,0.60,0.64,0.71,0.78,0.84, ...
       0.88,0.90,0.89,0.88,0.87,0.87,0.88,0.91,0.94,0.97, ...
       1.00,0.96,0.87,0.75];

% Six core connected contingencies selected by deterministic N-1 screening.
% Each row is [from_bus to_bus].
core_outages = [
    26 30
     8  5
    26 25
    23 25
    94 100
    94 96
];

% Principal branches reported even when they are not the most stressed.
target_pairs = [
    26 30
    94 100
    23 25
    25 26
    92 100
];

weather_dir = 'multisite_weather';
weather_mat = fullfile(weather_dir,'multisite_weather_model.mat');
solar_temporal_csv = fullfile(weather_dir,'solar_temporal_corr.csv');
case_file = './case_pglib_opf_case118_ieee.m';

outdir = fullfile('results_paper','probabilistic_n1');
if ~exist(outdir,'dir'), mkdir(outdir); end

mpopt = mpoption('verbose',0,'out.all',0);

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


fprintf('=== Probabilistic multi-site N-1 pilot ===\n');
fprintf('Scenarios: 4 | contingencies: %d | trajectories: %d | hours: %d\n', ...
    size(core_outages,1),Ns,T);
fprintf('Requested stochastic OPFs: %d\n',4*size(core_outages,1)*Ns*T);

%% =========================================================================
%% 1. LOAD AND VALIDATE MULTI-SITE WEATHER MODEL
%% =========================================================================
if ~exist(weather_mat,'file'), error('Missing %s.',weather_mat); end
if ~exist(solar_temporal_csv,'file'), error('Missing %s.',solar_temporal_csv); end

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
if nzone ~= 6, error('Expected six weather zones, found %d.',nzone); end

Tsol = readtable(solar_temporal_csv,'TextType','string');
solar_rho_ar1 = NaN(nzone,1);
for iz=1:nzone
    mask = startsWith(Tsol.zone,sprintf('Z%d_',iz)) & Tsol.lag_hours==1;
    if nnz(mask)~=1, error('Expected one lag-1 solar correlation for zone %d.',iz); end
    solar_rho_ar1(iz) = Tsol.correlation(mask);
end

fprintf('\n[1] Generating common multi-site renewable master sample...\n');
Xsolar = generate_vector_ar1(Ns,T,W.solar_corr_gaussian,solar_rho_ar1);
Xwind  = generate_vector_ar1(Ns,T,W.wind_corr_gaussian,W.wind_rho_ar1(:));
Usolar = min(max(normcdf(Xsolar),1e-8),1-1e-8);
Uwind  = min(max(normcdf(Xwind),1e-8),1-1e-8);

G_solar_zone = zeros(Ns,T,nzone);
P_wind_zone  = zeros(Ns,T,nzone);
v_ci=3.0; v_r=12.0; v_co=25.0;
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

%% =========================================================================
%% 2. NETWORK, SCENARIOS, OUTAGE AND TARGET INDICES
%% =========================================================================
fprintf('[2] Loading network and validating outages...\n');
mpc0 = loadcase(case_file);
nb=size(mpc0.bus,1); nbr=size(mpc0.branch,1);
Ppk=sum(mpc0.bus(:,3)); Pd_frac=mpc0.bus(:,3)/Ppk;
Qd0=mpc0.bus(:,4); Ppk_scaled=Ppk*peak_frac;
Prat=mpc0.branch(:,6);
if any(~isfinite(Prat)|Prat<=0), error('All RATE_A values must be positive and finite.'); end

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
sc(3)=make_scenario('S3_Low_Dist','Low/Distributed-diversified',P_solar_low,P_wind_low,DIST);
sc(4)=make_scenario('S4_High_Dist','High/Distributed-diversified',P_solar_high,P_wind_high,DIST);

ncont=size(core_outages,1);
outage_idx=zeros(ncont,1);
outage_label=strings(ncont,1);
for c=1:ncont
    outage_idx(c)=find_branch(mpc0,core_outages(c,1),core_outages(c,2));
    outage_label(c)=sprintf('%d-%d',core_outages(c,1),core_outages(c,2));
    if mpc0.branch(outage_idx(c),11)==0
        error('Requested outage %s is already out of service.',outage_label(c));
    end
    mt=mpc0; mt.branch(outage_idx(c),11)=0;
    if is_islanded(mt)
        error('Requested outage %s islands the network; remove it from core_outages.',outage_label(c));
    end
end

target_idx=zeros(size(target_pairs,1),1);
target_label=strings(size(target_pairs,1),1);
for k=1:size(target_pairs,1)
    target_idx(k)=find_branch(mpc0,target_pairs(k,1),target_pairs(k,2));
    target_label(k)=sprintf('%d-%d',target_pairs(k,1),target_pairs(k,2));
end

% Parallel pool and MATPOWER path propagation.
if use_parallel && license('test','Distrib_Computing_Toolbox')
    p=gcp('nocreate');
    if isempty(p), parpool(requested_workers); end
    mp_dir=fileparts(which('runopf'));
    if ~isempty(mp_dir), pctRunOnAll(sprintf('addpath(genpath(''%s''))',mp_dir)); end
else
    use_parallel=false;
end

%% =========================================================================
%% 3. DETERMINISTIC NO-RES BASELINE FOR EACH CONTINGENCY
%% =========================================================================
fprintf('\n[3] Running contingency-specific deterministic no-RES baselines...\n');
base=[];
base_rows={};
for c=1:ncont
    mt0=mpc0; mt0.branch(outage_idx(c),11)=0;
    bout=run_baseline_contingency(mt0,Pd_frac,Qd0,Ppk_scaled,Ppk,d_t,Prat, ...
        Tau,Vscreen_min,Vscreen_max,mpopt);
    % MATLAB cannot assign a populated struct into struct([]) because the
    % structures have dissimilar field sets. Preallocate from the first
    % returned baseline structure, then fill the remaining elements.
    if c==1
        base=repmat(bout,ncont,1);
    else
        base(c)=bout;
    end

    for j=1:numel(Tau)
        [maxlcp,ib]=max(bout.LCP(:,j));
        base_rows(end+1,:)={c,outage_idx(c),outage_label(c),Tau(j), ...
            bout.success_hours,maxlcp,ib,mpc0.branch(ib,1),mpc0.branch(ib,2), ...
            mean(bout.cost,'omitnan'),mean(bout.loss,'omitnan'), ...
            max(bout.loading_pct(:),[],'omitnan'), ...
            min(bout.V(:),[],'omitnan'),max(bout.V(:),[],'omitnan')}; %#ok<AGROW>
    end
    fprintf('  outage %s: %d/%d successful hours\n',outage_label(c),bout.success_hours,T);
end
Tbase=cell2table(base_rows,'VariableNames', ...
    {'contingency_id','outage_branch_idx','outage_branch','threshold', ...
     'successful_hours','max_base_LCP_pct','top_branch_idx','top_from_bus','top_to_bus', ...
     'mean_cost','mean_loss_MW','max_loading_pct','min_voltage_pu','max_voltage_pu'});
writetable(Tbase,fullfile(outdir,'n1_prob_contingency_baselines.csv'));

%% =========================================================================
%% 4. STOCHASTIC MULTI-SITE N-1 SIMULATIONS
%% =========================================================================
fprintf('\n[4] Running stochastic scenario-contingency combinations...\n');
R=cell(4,ncont);
run_clock=tic;

for isc=1:4
    cfg=sc(isc);
    [mpc_res,idx_s,idx_w]=add_res_generators(mpc0,cfg);
    [Psolar_units,Pwind_units]=scenario_available_power(cfg,G_solar_zone,P_wind_zone);

    for c=1:ncont
        fprintf('[4] %s | outage %s (%d/%d)\n',cfg.name,outage_label(c),c,ncont);
        mpc_cont=mpc_res;
        mpc_cont.branch(outage_idx(c),11)=0;

        exceed_count=zeros(Ns,nbr,numel(Tau),'uint8');
        vviol_count=zeros(Ns,nb,'uint8');
        daily_cost=NaN(Ns,1); daily_loss=NaN(Ns,1);
        daily_curt_s=NaN(Ns,1); daily_curt_w=NaN(Ns,1);
        max_curt_s=NaN(Ns,1); max_curt_w=NaN(Ns,1);
        curtailed_hours=NaN(Ns,1); successful_hours=zeros(Ns,1,'uint8');
        solve_time_s=NaN(Ns,1); max_loading_pct=NaN(Ns,1);
        min_voltage_pu=NaN(Ns,1); max_voltage_pu=NaN(Ns,1);
        target_exceed80=zeros(Ns,numel(target_idx),'uint8');
        target_max_loading=NaN(Ns,numel(target_idx));

        if use_parallel
            parfor s=1:Ns
                day=run_day_n1(mpc_cont,squeeze(Psolar_units(s,:,:)), ...
                    squeeze(Pwind_units(s,:,:)),idx_s,idx_w,Pd_frac,Qd0, ...
                    Ppk_scaled,Ppk,d_t,Prat,Tau,Vscreen_min,Vscreen_max, ...
                    curt_tol_MW,target_idx,outage_idx(c),mpopt);
                exceed_count(s,:,:)=day.exceed_count;
                vviol_count(s,:)=day.vviol_count;
                daily_cost(s)=day.daily_cost; daily_loss(s)=day.daily_loss_MWh;
                daily_curt_s(s)=day.daily_solar_curt_MWh;
                daily_curt_w(s)=day.daily_wind_curt_MWh;
                max_curt_s(s)=day.max_solar_curt_MW;
                max_curt_w(s)=day.max_wind_curt_MW;
                curtailed_hours(s)=day.curtailed_hours;
                successful_hours(s)=day.successful_hours;
                solve_time_s(s)=day.total_solve_time_s;
                max_loading_pct(s)=day.max_loading_pct;
                min_voltage_pu(s)=day.min_voltage_pu;
                max_voltage_pu(s)=day.max_voltage_pu;
                target_exceed80(s,:)=day.target_exceed80;
                target_max_loading(s,:)=day.target_max_loading;
            end
        else
            for s=1:Ns
                day=run_day_n1(mpc_cont,squeeze(Psolar_units(s,:,:)), ...
                    squeeze(Pwind_units(s,:,:)),idx_s,idx_w,Pd_frac,Qd0, ...
                    Ppk_scaled,Ppk,d_t,Prat,Tau,Vscreen_min,Vscreen_max, ...
                    curt_tol_MW,target_idx,outage_idx(c),mpopt);
                exceed_count(s,:,:)=day.exceed_count;
                vviol_count(s,:)=day.vviol_count;
                daily_cost(s)=day.daily_cost; daily_loss(s)=day.daily_loss_MWh;
                daily_curt_s(s)=day.daily_solar_curt_MWh;
                daily_curt_w(s)=day.daily_wind_curt_MWh;
                max_curt_s(s)=day.max_solar_curt_MW;
                max_curt_w(s)=day.max_wind_curt_MW;
                curtailed_hours(s)=day.curtailed_hours;
                successful_hours(s)=day.successful_hours;
                solve_time_s(s)=day.total_solve_time_s;
                max_loading_pct(s)=day.max_loading_pct;
                min_voltage_pu(s)=day.min_voltage_pu;
                max_voltage_pu(s)=day.max_voltage_pu;
                target_exceed80(s,:)=day.target_exceed80;
                target_max_loading(s,:)=day.target_max_loading;
                if mod(s,50)==0, fprintf('    s=%d/%d\n',s,Ns); end
            end
        end

        valid=double(successful_hours)==T;
        out=struct();
        out.exceed_count=exceed_count;
        out.vviol_count=vviol_count;
        out.daily_cost=daily_cost; out.daily_loss=daily_loss;
        out.daily_curt_s=daily_curt_s; out.daily_curt_w=daily_curt_w;
        out.max_curt_s=max_curt_s; out.max_curt_w=max_curt_w;
        out.curtailed_hours=curtailed_hours;
        out.successful_hours=successful_hours; out.solve_time_s=solve_time_s;
        out.max_loading_pct=max_loading_pct;
        out.min_voltage_pu=min_voltage_pu; out.max_voltage_pu=max_voltage_pu;
        out.target_exceed80=target_exceed80;
        out.target_max_loading=target_max_loading;
        out.valid=valid;
        R{isc,c}=out;

        if export_daily_metrics
            Tdaily=table((1:Ns)',double(successful_hours),solve_time_s,daily_cost,daily_loss, ...
                daily_curt_s,daily_curt_w,max_curt_s,max_curt_w,curtailed_hours, ...
                max_loading_pct,min_voltage_pu,max_voltage_pu, ...
                'VariableNames',{'trajectory','successful_hours','solve_time_s','daily_cost', ...
                'daily_loss_MWh','solar_curt_MWh','wind_curt_MWh','max_solar_curt_MW', ...
                'max_wind_curt_MW','curtailed_hours','max_loading_pct','min_voltage_pu','max_voltage_pu'});
            writetable(Tdaily,fullfile(outdir,sprintf('n1_prob_daily_metrics_%s_out_%s.csv', ...
                cfg.name,strrep(outage_label(c),'-','_'))));
        end
    end
end
fprintf('[4] All stochastic runs completed in %.2f h.\n',toc(run_clock)/3600);

%% =========================================================================
%% 5. BOOTSTRAP SUMMARIES AND TARGET-BRANCH METRICS
%% =========================================================================
fprintf('\n[5] Computing trajectory-bootstrap confidence intervals...\n');
rng(bootstrap_seed,'twister');

opf_rows={}; system_rows={}; target_rows={}; rank_rows={}; paired_rows={};

for isc=1:4
    for c=1:ncont
        out=R{isc,c}; valid=out.valid; nv=nnz(valid);
        requested=Ns*T; successful=sum(double(out.successful_hours));
        opf_rows(end+1,:)={sc(isc).name,sc(isc).label,c,outage_idx(c),outage_label(c), ...
            requested,successful,requested-successful,100*successful/requested, ...
            mean(out.solve_time_s(valid),'omitnan')/T,sum(out.solve_time_s,'omitnan')}; %#ok<AGROW>

        metric_defs={ ...
            'daily_cost',out.daily_cost,'currency/day'; ...
            'daily_loss_MWh',out.daily_loss,'MWh/day'; ...
            'solar_curt_MWh',out.daily_curt_s,'MWh/day'; ...
            'wind_curt_MWh',out.daily_curt_w,'MWh/day'; ...
            'curtailed_hours',out.curtailed_hours,'h/day'; ...
            'curtailment_trajectory_probability_pct',100*double(out.curtailed_hours>0),'%'; ...
            'max_loading_pct',out.max_loading_pct,'%'; ...
            'min_voltage_pu',out.min_voltage_pu,'p.u.'; ...
            'max_voltage_pu',out.max_voltage_pu,'p.u.'};
        for im=1:size(metric_defs,1)
            x=metric_defs{im,2}; x=x(valid);
            [est,se,lo,hi]=bootstrap_mean(x,B,ci_level);
            system_rows(end+1,:)={sc(isc).name,sc(isc).label,c,outage_idx(c),outage_label(c), ...
                nv,metric_defs{im,1},est,se,lo,hi,metric_defs{im,3}}; %#ok<AGROW>
        end

        % Network-wide LCP and Delta-LCP at all thresholds, retain top branch.
        for jt=1:numel(Tau)
            counts=double(squeeze(out.exceed_count(valid,:,jt)));
            LCP=mean(counts,1)/T*100;
            [sevec,lvec,hvec]=bootstrap_lcp_matrix(counts,T,B,ci_level);
            dLCP=LCP'-base(c).LCP(:,jt);
            active=true(nbr,1); active(outage_idx(c))=false;
            tmp=dLCP; tmp(~active)=-Inf;
            [maxd,ib]=max(tmp);
            rank_rows(end+1,:)={sc(isc).name,sc(isc).label,c,outage_idx(c),outage_label(c), ...
                Tau(jt),ib,mpc0.branch(ib,1),mpc0.branch(ib,2),Prat(ib), ...
                base(c).LCP(ib,jt),LCP(ib),maxd,sevec(ib),lvec(ib),hvec(ib)}; %#ok<AGROW>
        end

        % Prespecified target branches at tau=0.80.
        counts80=double(squeeze(out.exceed_count(valid,:,2)));
        LCP80=mean(counts80,1)/T*100;
        [se80,lo80,hi80]=bootstrap_lcp_matrix(counts80,T,B,ci_level);
        for k=1:numel(target_idx)
            ib=target_idx(k);
            is_out=(ib==outage_idx(c));
            if is_out
                lcp=NaN; se=NaN; lo=NaN; hi=NaN; base_lcp=NaN; delta=NaN;
                meanmax=NaN; maxmax=NaN;
            else
                lcp=LCP80(ib); se=se80(ib); lo=lo80(ib); hi=hi80(ib);
                base_lcp=base(c).LCP(ib,2); delta=lcp-base_lcp;
                x=out.target_max_loading(valid,k);
                meanmax=mean(x,'omitnan'); maxmax=max(x,[],'omitnan');
            end
            target_rows(end+1,:)={sc(isc).name,sc(isc).label,c,outage_idx(c),outage_label(c), ...
                target_label(k),ib,is_out,base_lcp,lcp,delta,se,lo,hi,meanmax,maxmax}; %#ok<AGROW>
        end
    end
end

%% =========================================================================
%% 6. PAIRED FACTORIAL CONTRASTS WITHIN EACH CONTINGENCY
%% =========================================================================
contrast_pairs=[1 3;2 4;2 1;4 3];
contrast_names=["S1-S3_low_conc_minus_dist";"S2-S4_high_conc_minus_dist"; ...
                "S2-S1_high_minus_low_conc";"S4-S3_high_minus_low_dist"];

for c=1:ncont
    for ip=1:size(contrast_pairs,1)
        ia=contrast_pairs(ip,1); ib=contrast_pairs(ip,2);
        A=R{ia,c}; D=R{ib,c};
        defs={ ...
            'daily_cost',A.daily_cost,D.daily_cost,'currency/day'; ...
            'daily_loss_MWh',A.daily_loss,D.daily_loss,'MWh/day'; ...
            'solar_curt_MWh',A.daily_curt_s,D.daily_curt_s,'MWh/day'; ...
            'wind_curt_MWh',A.daily_curt_w,D.daily_curt_w,'MWh/day'; ...
            'curtailed_hours',A.curtailed_hours,D.curtailed_hours,'h/day'};
        for im=1:size(defs,1)
            [est,se,lo,hi]=bootstrap_paired_difference(defs{im,2},defs{im,3},B,ci_level);
            paired_rows(end+1,:)={c,outage_idx(c),outage_label(c),contrast_names(ip), ...
                defs{im,1},est,se,lo,hi,defs{im,4}}; %#ok<AGROW>
        end

        % Paired LCP contrast at branch 26-30, unless that branch is outaged.
        k26=find(target_idx==find_branch(mpc0,26,30),1);
        if outage_idx(c)~=target_idx(k26)
            xa=double(A.target_exceed80(:,k26))/T*100;
            xb=double(D.target_exceed80(:,k26))/T*100;
            [est,se,lo,hi]=bootstrap_paired_difference(xa,xb,B,ci_level);
            paired_rows(end+1,:)={c,outage_idx(c),outage_label(c),contrast_names(ip), ...
                'LCP80_branch_26_30_pp',est,se,lo,hi,'percentage points'}; %#ok<AGROW>
        end
    end
end

%% =========================================================================
%% 7. EXPORT TABLES
%% =========================================================================
Topf=cell2table(opf_rows,'VariableNames', ...
    {'scenario','scenario_label','contingency_id','outage_branch_idx','outage_branch', ...
     'requested_OPFs','successful_OPFs','failed_OPFs','success_rate_pct', ...
     'mean_solve_time_s_per_OPF','total_solve_time_s'});
writetable(Topf,fullfile(outdir,'n1_prob_opf_summary.csv'));

Tsys=cell2table(system_rows,'VariableNames', ...
    {'scenario','scenario_label','contingency_id','outage_branch_idx','outage_branch', ...
     'N_valid','metric','estimate','bootstrap_se','ci_low','ci_high','unit'});
writetable(Tsys,fullfile(outdir,'n1_prob_system_metrics.csv'));

Ttarget=cell2table(target_rows,'VariableNames', ...
    {'scenario','scenario_label','contingency_id','outage_branch_idx','outage_branch', ...
     'target_branch','target_branch_idx','target_is_outaged','base_LCP80_pct','LCP80_pct', ...
     'delta_LCP80_pp','bootstrap_se_pp','ci_low_pct','ci_high_pct', ...
     'mean_daily_max_loading_pct','maximum_observed_loading_pct'});
writetable(Ttarget,fullfile(outdir,'n1_prob_target_branch_metrics.csv'));

Trank=cell2table(rank_rows,'VariableNames', ...
    {'scenario','scenario_label','contingency_id','outage_branch_idx','outage_branch', ...
     'threshold','top_delta_branch_idx','top_from_bus','top_to_bus','rating_MVA', ...
     'base_LCP_pct','RES_LCP_pct','delta_LCP_pp','bootstrap_se_pp', ...
     'LCP_ci_low_pct','LCP_ci_high_pct'});

% Add a transparent severity rank based on Delta-LCP80, mean losses and
% curtailment probability. Individual metrics remain available separately.
mask80=abs(Trank.threshold-0.80)<1e-12;
Trank80=Trank(mask80,:);
loss_lookup=zeros(height(Trank80),1); curt_lookup=zeros(height(Trank80),1);
for i=1:height(Trank80)
    m1=strcmp(Tsys.scenario,Trank80.scenario(i)) & ...
       Tsys.contingency_id==Trank80.contingency_id(i) & strcmp(Tsys.metric,'daily_loss_MWh');
    m2=strcmp(Tsys.scenario,Trank80.scenario(i)) & ...
       Tsys.contingency_id==Trank80.contingency_id(i) & strcmp(Tsys.metric,'curtailment_trajectory_probability_pct');
    loss_lookup(i)=Tsys.estimate(find(m1,1));
    curt_lookup(i)=Tsys.estimate(find(m2,1));
end
Trank80.mean_daily_loss_MWh=loss_lookup;
Trank80.curtailment_probability_pct=curt_lookup;
Trank80.severity_score=mean([minmax01(Trank80.delta_LCP_pp), ...
    minmax01(Trank80.mean_daily_loss_MWh),minmax01(Trank80.curtailment_probability_pct)],2,'omitnan');
[~,ord]=sort(Trank80.severity_score,'descend');
Trank80=Trank80(ord,:);
Trank80.severity_rank=(1:height(Trank80))';
Trank80=movevars(Trank80,'severity_rank','Before','scenario');
writetable(Trank80,fullfile(outdir,'n1_prob_scenario_contingency_ranking.csv'));

Tpaired=cell2table(paired_rows,'VariableNames', ...
    {'contingency_id','outage_branch_idx','outage_branch','contrast', ...
     'metric','estimate','bootstrap_se','ci_low','ci_high','unit'});
writetable(Tpaired,fullfile(outdir,'n1_prob_paired_scenario_differences.csv'));

if save_compact_mat
    save(fullfile(outdir,'n1_probabilistic_multisite_pilot_compact.mat'), ...
        'R','base','sc','core_outages','outage_idx','outage_label', ...
        'target_pairs','target_idx','target_label','G_solar_zone','P_wind_zone', ...
        'Ns','T','B','Tau','peak_frac','d_t','W','solar_rho_ar1','-v7.3');
end

fprintf('\n=== N-1 probabilistic multi-site pilot complete ===\n');
fprintf('Outputs: %s\n',outdir);
fprintf('Total stochastic OPFs requested: %d\n',4*ncont*Ns*T);
fprintf('Overall success rate: %.4f%%\n',100*sum(Topf.successful_OPFs)/sum(Topf.requested_OPFs));

%% =========================================================================
%% LOCAL FUNCTIONS
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

function [Ps,Pw]=scenario_available_power(cfg,Gzone,Wzone)
Ns=size(Gzone,1); T=size(Gzone,2);
Ps=zeros(Ns,T,numel(cfg.inj.solar_buses)); Pw=zeros(Ns,T,numel(cfg.inj.wind_buses));
for i=1:size(Ps,3)
    Ps(:,:,i)=cfg.P_solar_max*cfg.inj.solar_frac(i)*Gzone(:,:,cfg.inj.solar_zones(i));
end
for i=1:size(Pw,3)
    Pw(:,:,i)=cfg.P_wind_max*cfg.inj.wind_frac(i)*Wzone(:,:,cfg.inj.wind_zones(i));
end
end

function X=generate_vector_ar1(Ns,T,R,rho)
R=nearest_corr(R,1e-10); rho=rho(:); A=diag(rho);
Q=project_psd(R-A*R*A',1e-10); L0=chol(R,'lower'); Lq=chol(Q,'lower');
X=zeros(Ns,T,numel(rho)); X(:,1,:)=randn(Ns,numel(rho))*L0';
for t=2:T
    prev=squeeze(X(:,t-1,:)); eps=randn(Ns,numel(rho))*Lq'; X(:,t,:)=prev*A'+eps;
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

function out=run_baseline_contingency(mpc,Pdf,Qd0,Ppks,Ppk,dt,Prat,Tau,Vmin,Vmax,mpopt)
T=numel(dt); nbr=size(mpc.branch,1); nb=size(mpc.bus,1);
out.S=NaN(T,nbr); out.V=NaN(T,nb); out.cost=NaN(T,1); out.loss=NaN(T,1); out.ok=false(T,1);
for t=1:T
    mt=mpc; Plt=Ppks*dt(t); mt.bus(:,3)=Pdf*Plt; mt.bus(:,4)=Qd0*(Plt/Ppk);
    r=runopf(mt,mpopt); out.ok(t)=logical(r.success);
    if r.success
        out.S(t,:)=max(hypot(r.branch(:,14),r.branch(:,15)), hypot(r.branch(:,16),r.branch(:,17)))'; out.V(t,:)=r.bus(:,8)';
        out.cost(t)=r.f; out.loss(t)=sum(r.branch(:,14)+r.branch(:,16));
    end
end
out.success_hours=nnz(out.ok); out.loading_pct=100*out.S./Prat';
out.LCP=NaN(nbr,numel(Tau));
for j=1:numel(Tau), out.LCP(:,j)=mean(out.S>Tau(j)*Prat',1,'omitnan')'*100; end
out.VVP=mean(out.V<Vmin | out.V>Vmax,1,'omitnan')'*100;
end

function day=run_day_n1(mpc,Ps,Pw,idx_s,idx_w,Pdf,Qd0,Ppks,Ppk,dt,Prat,Tau,Vmin,Vmax,curt_tol,target_idx,outage_idx,mpopt)
T=numel(dt); nbr=size(mpc.branch,1); nb=size(mpc.bus,1);
if isvector(Ps), Ps=Ps(:); end
if isvector(Pw), Pw=Pw(:); end
if size(Ps,1)~=T, Ps=Ps'; end
if size(Pw,1)~=T, Pw=Pw'; end
ex=false(T,nbr,numel(Tau)); vv=false(T,nb); ok=false(T,1);
cost=NaN(T,1); loss=NaN(T,1); cs=NaN(T,1); cw=NaN(T,1); st=NaN(T,1);
maxload=NaN(T,1); minV=NaN(T,1); maxV=NaN(T,1);
targ_ex=false(T,numel(target_idx)); targ_load=NaN(T,numel(target_idx));
for t=1:T
    mt=mpc; Plt=Ppks*dt(t); mt.bus(:,3)=Pdf*Plt; mt.bus(:,4)=Qd0*(Plt/Ppk);
    for i=1:numel(idx_s), mt.gen(idx_s(i),9)=Ps(t,i); end
    for i=1:numel(idx_w), mt.gen(idx_w(i),9)=Pw(t,i); end
    ts=tic; r=runopf(mt,mpopt); st(t)=toc(ts); ok(t)=logical(r.success);
    if r.success
        S=max(hypot(r.branch(:,14),r.branch(:,15)), hypot(r.branch(:,16),r.branch(:,17))); V=r.bus(:,8); loading=100*S./Prat;
        loading(outage_idx)=NaN;
        for j=1:numel(Tau), ex(t,:,j)=(S>Tau(j)*Prat)'; end
        ex(t,outage_idx,:)=false; vv(t,:)=(V<Vmin|V>Vmax)';
        cost(t)=r.f; loss(t)=sum(r.branch(:,14)+r.branch(:,16));
        cs(t)=sum(max(0,Ps(t,:)'-r.gen(idx_s,2))); cw(t)=sum(max(0,Pw(t,:)'-r.gen(idx_w,2)));
        if cs(t)<curt_tol, cs(t)=0; end
        if cw(t)<curt_tol, cw(t)=0; end
        maxload(t)=max(loading,[],'omitnan'); minV(t)=min(V); maxV(t)=max(V);
        for k=1:numel(target_idx)
            if target_idx(k)~=outage_idx
                targ_ex(t,k)=S(target_idx(k))>0.80*Prat(target_idx(k));
                targ_load(t,k)=loading(target_idx(k));
            end
        end
    end
end
if all(ok)
    day.exceed_count=uint8(squeeze(sum(ex,1))); day.vviol_count=uint8(sum(vv,1));
    day.daily_cost=sum(cost); day.daily_loss_MWh=sum(loss);
    day.daily_solar_curt_MWh=sum(cs); day.daily_wind_curt_MWh=sum(cw);
    day.max_solar_curt_MW=max(cs); day.max_wind_curt_MW=max(cw);
    day.curtailed_hours=nnz((cs+cw)>curt_tol); day.max_loading_pct=max(maxload);
    day.min_voltage_pu=min(minV); day.max_voltage_pu=max(maxV);
    day.target_exceed80=uint8(sum(targ_ex,1)); day.target_max_loading=max(targ_load,[],1,'omitnan');
else
    day.exceed_count=zeros(nbr,numel(Tau),'uint8'); day.vviol_count=zeros(1,nb,'uint8');
    day.daily_cost=NaN; day.daily_loss_MWh=NaN; day.daily_solar_curt_MWh=NaN; day.daily_wind_curt_MWh=NaN;
    day.max_solar_curt_MW=NaN; day.max_wind_curt_MW=NaN; day.curtailed_hours=NaN;
    day.max_loading_pct=NaN; day.min_voltage_pu=NaN; day.max_voltage_pu=NaN;
    day.target_exceed80=zeros(1,numel(target_idx),'uint8'); day.target_max_loading=NaN(1,numel(target_idx));
end
day.successful_hours=uint8(nnz(ok)); day.total_solve_time_s=sum(st,'omitnan');
end

function [estimate,se,lo,hi]=bootstrap_mean(x,B,ci_level)
x=x(:); x=x(isfinite(x)); n=numel(x); estimate=mean(x);
if n<2, se=NaN; lo=NaN; hi=NaN; return; end
boot=zeros(B,1); for b=1:B, boot(b)=mean(x(randi(n,n,1))); end
se=std(boot,0); a=(1-ci_level)/2; q=prctile(boot,[100*a 100*(1-a)]); lo=q(1); hi=q(2);
end

function [se,lo,hi]=bootstrap_lcp_matrix(counts,T,B,ci_level)
[n,nbr]=size(counts); boot=zeros(B,nbr,'single');
for b=1:B, boot(b,:)=single(mean(counts(randi(n,n,1),:),1)/T*100); end
se=std(double(boot),0,1); a=(1-ci_level)/2; q=prctile(double(boot),[100*a 100*(1-a)],1); lo=q(1,:); hi=q(2,:);
end

function [estimate,se,lo,hi]=bootstrap_paired_difference(xa,xb,B,ci_level)
xa=xa(:); xb=xb(:); valid=isfinite(xa)&isfinite(xb); d=xa(valid)-xb(valid); n=numel(d);
if n==0, estimate=NaN; se=NaN; lo=NaN; hi=NaN; return; end
estimate=mean(d); boot=zeros(B,1); for b=1:B, boot(b)=mean(d(randi(n,n,1))); end
se=std(boot,0); a=(1-ci_level)/2; q=prctile(boot,[100*a 100*(1-a)]); lo=q(1); hi=q(2);
end

function idx=find_branch(mpc,fbus,tbus)
idx=find((mpc.branch(:,1)==fbus & mpc.branch(:,2)==tbus) | ...
         (mpc.branch(:,1)==tbus & mpc.branch(:,2)==fbus),1);
if isempty(idx), error('Branch %d-%d not found.',fbus,tbus); end
end

function tf=is_islanded(mpc)
active=mpc.branch(:,11)~=0; G=graph(mpc.branch(active,1),mpc.branch(active,2));
% IEEE bus numbers are consecutive 1..118, so graph() includes all buses up to max ID.
bins=conncomp(G); tf=numel(unique(bins))>1;
end

function y=minmax01(x)
x=double(x(:)); finite=isfinite(x); y=zeros(size(x));
if nnz(finite)==0, y(:)=NaN; return; end
mn=min(x(finite)); mx=max(x(finite));
if mx>mn, y(finite)=(x(finite)-mn)/(mx-mn); else, y(finite)=0; end
y(~finite)=NaN;
end
