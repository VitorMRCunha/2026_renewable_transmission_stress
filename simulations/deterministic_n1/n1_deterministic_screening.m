%% =========================================================================
%% n1_deterministic_screening_v1.m
%% Deterministic N-1 branch-outage screening — PGLib IEEE 118-bus
%%
%% Purpose
%%   1) Solve the intact no-RES 24-hour AC-OPF baseline.
%%   2) Screen every in-service single-branch outage.
%%   3) Detect islanding before attempting AC-OPF.
%%   4) Record solver reliability, loading, voltage, losses and cost.
%%   5) Export transparent severity rankings and a critical-contingency shortlist.
%%
%% This is a contingency-screening study, not a full SC-AC-OPF formulation.
%% Branch thermal limits and bus voltage limits remain enforced by MATPOWER.
%% Therefore, for converged OPFs, ranking focuses on proximity to limits,
%% threshold exposure, redispatch cost, losses, voltage stress and infeasibility.
%%
%% Requires
%%   - MATPOWER 7.1+
%%   - case_pglib_opf_case118_ieee.m in the MATLAB path/current folder
%%   - Parallel Computing Toolbox is optional
%%
%% Output folder
%%   results_paper/deterministic_n1/
%%
%% Authors: Azevedo & Cunha, ISEP/IPP
%% =========================================================================

clear; clc; close all;
rng(42, 'twister');

%% =========================================================================
%% 0. USER CONFIGURATION
%% =========================================================================
case_file = 'case_pglib_opf_case118_ieee.m';

T = 24;
peak_frac = 0.70;
d_t = [0.67,0.63,0.60,0.59,0.59,0.60,0.64,0.71,0.78,0.84,0.88,0.90, ...
       0.89,0.88,0.87,0.87,0.88,0.91,0.94,0.97,1.00,0.96,0.87,0.75];

% Post-processing screening band. The OPF itself uses case-specific limits.
Vscreen_min = 0.95;
Vscreen_max = 1.05;

% Loading thresholds used for contingency severity screening.
load_thresholds_pct = [80 90 95 99];

% Number of contingencies retained in each metric-specific ranking.
top_per_metric = 10;

% Maximum size of union shortlist for later probabilistic N-1 analysis.
shortlist_max = 12;

% Parallel execution. Uses numeric MATPOWER column indices only.
use_parallel = true;
requested_workers = 6;  % set [] to use MATLAB default

% Output paths.
outdir = fullfile('results_paper', 'deterministic_n1');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% =========================================================================
%% 1. MATPOWER SETUP AND CASE VALIDATION
%% =========================================================================
assert(exist('runopf','file') == 2, ...
    'MATPOWER runopf was not found. Add MATPOWER to the MATLAB path.');
assert(exist(case_file,'file') == 2, ...
    'Case file not found: %s', case_file);

mpc0 = loadcase(case_file);

% Numeric MATPOWER indices (parfor-safe; no define_constants inside loops).
BUS_I     = 1;
BUS_TYPE  = 2;
PD        = 3;
QD        = 4;
VM        = 8;
VMAX      = 12;
VMIN      = 13;

F_BUS     = 1;
T_BUS     = 2;
RATE_A    = 6;
BR_STATUS = 11;
PF        = 14;
QF        = 15;
PT        = 16;
QT        = 17;

nb  = size(mpc0.bus,1);
nbr = size(mpc0.branch,1);
ng  = size(mpc0.gen,1);

assert(numel(d_t) == T, 'd_t must contain exactly T=%d values.', T);

% PGLib RateA validation.
ratings = mpc0.branch(:, RATE_A);
if any(ratings <= 0 | ~isfinite(ratings))
    bad = find(ratings <= 0 | ~isfinite(ratings));
    error('Invalid RateA on %d branch(es). First bad branch index: %d.', ...
        numel(bad), bad(1));
end

active_branches = find(mpc0.branch(:, BR_STATUS) > 0);
ncont = numel(active_branches);

% Demand scaling used in the earlier paper experiments.
Ppk = sum(mpc0.bus(:, PD));
assert(Ppk > 0, 'Total case demand must be positive.');
Pd_frac = mpc0.bus(:, PD) / Ppk;
Qd_base = mpc0.bus(:, QD);
Ppk_scaled = Ppk * peak_frac;

mpopt = mpoption('verbose',0, 'out.all',0);

% Start/reuse parallel pool.
parallel_used = false;
if use_parallel && license('test','Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if isempty(pool)
            if isempty(requested_workers)
                pool = parpool;
            else
                pool = parpool(requested_workers);
            end
        end
        mp_dir = fileparts(which('runopf'));
        if ~isempty(mp_dir)
            pctRunOnAll(sprintf('addpath(genpath(''%s''))', mp_dir));
        end
        parallel_used = true;
        fprintf('[PAR] Using %d workers.\n', pool.NumWorkers);
    catch ME
        warning('Parallel pool unavailable; falling back to serial: %s', ME.message);
        parallel_used = false;
    end
end

fprintf('\n=== Deterministic N-1 Screening ===\n');
fprintf('Case: %s | buses=%d | branches=%d | active outages=%d | hours=%d\n', ...
    case_file, nb, nbr, ncont, T);
fprintf('Peak-demand scaling: %.2f\n', peak_frac);

%% =========================================================================
%% 2. IDENTIFY TARGET CORRIDORS
%% =========================================================================
target_94100 = find_branch_idx(mpc0, 94, 100);
target_2630  = find_branch_idx(mpc0, 26, 30);

fprintf('Target branch 94-100 index: %d\n', target_94100);
fprintf('Target branch 26-30  index: %d\n', target_2630);

%% =========================================================================
%% 3. INTACT 24-HOUR NO-RES BASELINE
%% =========================================================================
fprintf('\n[BASE] Running intact 24-hour no-RES AC-OPF...\n');

base_success = false(T,1);
base_solve_time = nan(T,1);
base_cost = nan(T,1);
base_loss = nan(T,1);
base_vm = nan(T,nb);
base_loading_pct = nan(T,nbr);

for t = 1:T
    mt = apply_hourly_load(mpc0, Pd_frac, Qd_base, Ppk_scaled, Ppk, d_t(t), PD, QD);
    tic_opf = tic;
    r = runopf(mt, mpopt);
    base_solve_time(t) = toc(tic_opf);

    if r.success
        base_success(t) = true;
        base_cost(t) = r.f;
        base_loss(t) = sum(r.branch(:,PF) + r.branch(:,PT));
        base_vm(t,:) = r.bus(:,VM)';
        base_loading_pct(t,:) = branch_loading_pct(r.branch, ratings, PF, QF, PT, QT)';
    end
end

if ~all(base_success)
    error('The intact baseline failed in %d of %d hours. Resolve before N-1 screening.', ...
        nnz(~base_success), T);
end

base_mean_cost = mean(base_cost);
base_mean_loss = mean(base_loss);
base_max_loading = max(base_loading_pct(:));
base_min_vm = min(base_vm(:));
base_max_vm = max(base_vm(:));

fprintf('[BASE] 24/24 converged | max loading=%.2f%% | mean loss=%.3f MW | mean cost=%.2f\n', ...
    base_max_loading, base_mean_loss, base_mean_cost);

% Export hourly intact baseline.
Tbase = table((1:T)', d_t(:), base_success, base_solve_time, base_cost, base_loss, ...
    max(base_loading_pct,[],2), min(base_vm,[],2), max(base_vm,[],2), ...
    base_loading_pct(:,target_94100), base_loading_pct(:,target_2630), ...
    'VariableNames', {'hour','demand_factor','success','solve_time_s','cost','loss_MW', ...
    'max_loading_pct','min_VM_pu','max_VM_pu','loading_94_100_pct','loading_26_30_pct'});
writetable(Tbase, fullfile(outdir,'n1_intact_baseline_hourly.csv'));

%% =========================================================================
%% 4. PRE-SCREEN ISLANDING FOR EVERY BRANCH OUTAGE
%% =========================================================================
fprintf('\n[TOPOLOGY] Checking all outages for islanding...\n');

is_islanding = false(ncont,1);
n_islands = ones(ncont,1);
smallest_island_buses = nan(ncont,1);
largest_island_buses = nan(ncont,1);

for ic = 1:ncont
    br = active_branches(ic);
    mt = mpc0;
    mt.branch(br, BR_STATUS) = 0;
    [is_islanding(ic), n_islands(ic), smallest_island_buses(ic), largest_island_buses(ic)] = ...
        topology_islanding(mt, BUS_I, F_BUS, T_BUS, BR_STATUS);
end

fprintf('[TOPOLOGY] Islanding outages: %d/%d.\n', nnz(is_islanding), ncont);

%% =========================================================================
%% 5. RUN DETERMINISTIC N-1 AC-OPF SCREENING
%% =========================================================================
fprintf('\n[N-1] Running %d branch outages x %d hours = %d possible OPFs...\n', ...
    ncont, T, ncont*T);

% Compact contingency-level metrics.
success_hours = zeros(ncont,1);
failed_hours = zeros(ncont,1);
total_solve_time = zeros(ncont,1);
mean_solve_time = nan(ncont,1);
max_solve_time = nan(ncont,1);
mean_cost = nan(ncont,1);
max_cost = nan(ncont,1);
mean_loss = nan(ncont,1);
max_loss = nan(ncont,1);
max_loading_pct = nan(ncont,1);
mean_system_max_loading_pct = nan(ncont,1);
max_loading_branch_idx = nan(ncont,1);
max_loading_hour = nan(ncont,1);

% Threshold counts/exposure.
hours_any_gt80 = zeros(ncont,1);
hours_any_gt90 = zeros(ncont,1);
hours_any_gt95 = zeros(ncont,1);
hours_any_gt99 = zeros(ncont,1);
branch_hours_gt80 = zeros(ncont,1);
branch_hours_gt90 = zeros(ncont,1);
branch_hours_gt95 = zeros(ncont,1);
branch_hours_gt99 = zeros(ncont,1);
stress_above80 = nan(ncont,1);
stress_above90 = nan(ncont,1);

% Voltage metrics.
min_vm = nan(ncont,1);
max_vm = nan(ncont,1);
max_case_voltage_violation_pu = nan(ncont,1);
case_voltage_violation_bus_hours = zeros(ncont,1);
screen_voltage_violation_bus_hours = zeros(ncont,1);
mean_abs_voltage_dev_from_intact = nan(ncont,1);
max_abs_voltage_dev_from_intact = nan(ncont,1);

% Target corridor metrics.
mean_loading_94100 = nan(ncont,1);
max_loading_94100 = nan(ncont,1);
mean_loading_2630 = nan(ncont,1);
max_loading_2630 = nan(ncont,1);

% Hourly detail arrays for reproducibility.
h_success = false(ncont,T);
h_solve_time = nan(ncont,T);
h_cost = nan(ncont,T);
h_loss = nan(ncont,T);
h_max_loading = nan(ncont,T);
h_min_vm = nan(ncont,T);
h_max_vm = nan(ncont,T);
h_loading_94100 = nan(ncont,T);
h_loading_2630 = nan(ncont,T);

% Run one contingency per parfor iteration to avoid nested parallel loops.
t_start = tic;
if parallel_used
    parfor ic = 1:ncont
        br = active_branches(ic);
        if is_islanding(ic)
            continue;
        end
        out = evaluate_contingency(mpc0, br, T, d_t, Pd_frac, Qd_base, ...
            Ppk_scaled, Ppk, ratings, base_vm, target_94100, target_2630, ...
            Vscreen_min, Vscreen_max, mpopt, ...
            PD, QD, VM, VMAX, VMIN, BR_STATUS, PF, QF, PT, QT);

        success_hours(ic) = out.success_hours;
        failed_hours(ic) = out.failed_hours;
        total_solve_time(ic) = out.total_solve_time;
        mean_solve_time(ic) = out.mean_solve_time;
        max_solve_time(ic) = out.max_solve_time;
        mean_cost(ic) = out.mean_cost;
        max_cost(ic) = out.max_cost;
        mean_loss(ic) = out.mean_loss;
        max_loss(ic) = out.max_loss;
        max_loading_pct(ic) = out.max_loading_pct;
        mean_system_max_loading_pct(ic) = out.mean_system_max_loading_pct;
        max_loading_branch_idx(ic) = out.max_loading_branch_idx;
        max_loading_hour(ic) = out.max_loading_hour;
        hours_any_gt80(ic) = out.hours_any_gt80;
        hours_any_gt90(ic) = out.hours_any_gt90;
        hours_any_gt95(ic) = out.hours_any_gt95;
        hours_any_gt99(ic) = out.hours_any_gt99;
        branch_hours_gt80(ic) = out.branch_hours_gt80;
        branch_hours_gt90(ic) = out.branch_hours_gt90;
        branch_hours_gt95(ic) = out.branch_hours_gt95;
        branch_hours_gt99(ic) = out.branch_hours_gt99;
        stress_above80(ic) = out.stress_above80;
        stress_above90(ic) = out.stress_above90;
        min_vm(ic) = out.min_vm;
        max_vm(ic) = out.max_vm;
        max_case_voltage_violation_pu(ic) = out.max_case_voltage_violation_pu;
        case_voltage_violation_bus_hours(ic) = out.case_voltage_violation_bus_hours;
        screen_voltage_violation_bus_hours(ic) = out.screen_voltage_violation_bus_hours;
        mean_abs_voltage_dev_from_intact(ic) = out.mean_abs_voltage_dev_from_intact;
        max_abs_voltage_dev_from_intact(ic) = out.max_abs_voltage_dev_from_intact;
        mean_loading_94100(ic) = out.mean_loading_94100;
        max_loading_94100(ic) = out.max_loading_94100;
        mean_loading_2630(ic) = out.mean_loading_2630;
        max_loading_2630(ic) = out.max_loading_2630;
        h_success(ic,:) = out.h_success;
        h_solve_time(ic,:) = out.h_solve_time;
        h_cost(ic,:) = out.h_cost;
        h_loss(ic,:) = out.h_loss;
        h_max_loading(ic,:) = out.h_max_loading;
        h_min_vm(ic,:) = out.h_min_vm;
        h_max_vm(ic,:) = out.h_max_vm;
        h_loading_94100(ic,:) = out.h_loading_94100;
        h_loading_2630(ic,:) = out.h_loading_2630;
    end
else
    for ic = 1:ncont
        br = active_branches(ic);
        if ~is_islanding(ic)
            out = evaluate_contingency(mpc0, br, T, d_t, Pd_frac, Qd_base, ...
                Ppk_scaled, Ppk, ratings, base_vm, target_94100, target_2630, ...
                Vscreen_min, Vscreen_max, mpopt, ...
                PD, QD, VM, VMAX, VMIN, BR_STATUS, PF, QF, PT, QT);

            success_hours(ic) = out.success_hours;
            failed_hours(ic) = out.failed_hours;
            total_solve_time(ic) = out.total_solve_time;
            mean_solve_time(ic) = out.mean_solve_time;
            max_solve_time(ic) = out.max_solve_time;
            mean_cost(ic) = out.mean_cost;
            max_cost(ic) = out.max_cost;
            mean_loss(ic) = out.mean_loss;
            max_loss(ic) = out.max_loss;
            max_loading_pct(ic) = out.max_loading_pct;
            mean_system_max_loading_pct(ic) = out.mean_system_max_loading_pct;
            max_loading_branch_idx(ic) = out.max_loading_branch_idx;
            max_loading_hour(ic) = out.max_loading_hour;
            hours_any_gt80(ic) = out.hours_any_gt80;
            hours_any_gt90(ic) = out.hours_any_gt90;
            hours_any_gt95(ic) = out.hours_any_gt95;
            hours_any_gt99(ic) = out.hours_any_gt99;
            branch_hours_gt80(ic) = out.branch_hours_gt80;
            branch_hours_gt90(ic) = out.branch_hours_gt90;
            branch_hours_gt95(ic) = out.branch_hours_gt95;
            branch_hours_gt99(ic) = out.branch_hours_gt99;
            stress_above80(ic) = out.stress_above80;
            stress_above90(ic) = out.stress_above90;
            min_vm(ic) = out.min_vm;
            max_vm(ic) = out.max_vm;
            max_case_voltage_violation_pu(ic) = out.max_case_voltage_violation_pu;
            case_voltage_violation_bus_hours(ic) = out.case_voltage_violation_bus_hours;
            screen_voltage_violation_bus_hours(ic) = out.screen_voltage_violation_bus_hours;
            mean_abs_voltage_dev_from_intact(ic) = out.mean_abs_voltage_dev_from_intact;
            max_abs_voltage_dev_from_intact(ic) = out.max_abs_voltage_dev_from_intact;
            mean_loading_94100(ic) = out.mean_loading_94100;
            max_loading_94100(ic) = out.max_loading_94100;
            mean_loading_2630(ic) = out.mean_loading_2630;
            max_loading_2630(ic) = out.max_loading_2630;
            h_success(ic,:) = out.h_success;
            h_solve_time(ic,:) = out.h_solve_time;
            h_cost(ic,:) = out.h_cost;
            h_loss(ic,:) = out.h_loss;
            h_max_loading(ic,:) = out.h_max_loading;
            h_min_vm(ic,:) = out.h_min_vm;
            h_max_vm(ic,:) = out.h_max_vm;
            h_loading_94100(ic,:) = out.h_loading_94100;
            h_loading_2630(ic,:) = out.h_loading_2630;
        end
        if mod(ic,20)==0 || ic==ncont
            fprintf('  contingency %d/%d | elapsed %.1f min\n', ic, ncont, toc(t_start)/60);
        end
    end
end
fprintf('[N-1] Screening completed in %.2f min.\n', toc(t_start)/60);

% Ensure failed-hour count is explicit for islanding cases.
failed_hours(is_islanding) = T;
success_rate_pct = 100 * success_hours / T;

%% =========================================================================
%% 6. DERIVED DELTAS, LABELS AND TRANSPARENT SEVERITY SCORE
%% =========================================================================
outage_branch_idx = active_branches(:);
outage_from_bus = mpc0.branch(outage_branch_idx,F_BUS);
outage_to_bus = mpc0.branch(outage_branch_idx,T_BUS);
outage_rating_MVA = ratings(outage_branch_idx);
outage_label = compose('%d-%d', outage_from_bus, outage_to_bus);

mean_cost_delta = mean_cost - base_mean_cost;
mean_cost_delta_pct = 100 * mean_cost_delta / base_mean_cost;
mean_loss_delta = mean_loss - base_mean_loss;
mean_loss_delta_pct = 100 * mean_loss_delta / base_mean_loss;

% Identify the branch responsible for maximum loading.
max_loading_from_bus = nan(ncont,1);
max_loading_to_bus = nan(ncont,1);
valid_idx = isfinite(max_loading_branch_idx) & max_loading_branch_idx >= 1;
max_loading_from_bus(valid_idx) = mpc0.branch(max_loading_branch_idx(valid_idx),F_BUS);
max_loading_to_bus(valid_idx) = mpc0.branch(max_loading_branch_idx(valid_idx),T_BUS);
max_loading_branch = strings(ncont,1);
max_loading_branch(valid_idx) = compose('%d-%d', ...
    max_loading_from_bus(valid_idx), max_loading_to_bus(valid_idx));

% Transparent normalized screening score. It is used only to organize the
% shortlist; individual physical metrics remain the primary reported results.
% Infeasible/non-convergent and islanding outages receive priority flags.
feasible_mask = ~is_islanding & success_hours > 0;
score_loading = minmax01(max_loading_pct, feasible_mask);
score_stress  = minmax01(stress_above90, feasible_mask);
score_voltage = minmax01(max_abs_voltage_dev_from_intact, feasible_mask);
score_cost    = minmax01(max(mean_cost_delta,0), feasible_mask);
score_loss    = minmax01(max(mean_loss_delta,0), feasible_mask);

severity_score = 100 * (0.35*score_loading + 0.25*score_stress + ...
    0.15*score_voltage + 0.15*score_cost + 0.10*score_loss);

% Explicit penalties only distinguish contingencies that could not be
% represented by a fully converged connected AC-OPF trajectory.
severity_score = severity_score + 200*(failed_hours > 0 & ~is_islanding) + 300*is_islanding;

%% =========================================================================
%% 7. MASTER CONTINGENCY SUMMARY CSV
%% =========================================================================
Tsummary = table( ...
    outage_branch_idx, outage_from_bus, outage_to_bus, outage_label, outage_rating_MVA, ...
    is_islanding, n_islands, smallest_island_buses, largest_island_buses, ...
    success_hours, failed_hours, success_rate_pct, total_solve_time, mean_solve_time, max_solve_time, ...
    mean_cost, mean_cost_delta, mean_cost_delta_pct, max_cost, ...
    mean_loss, mean_loss_delta, mean_loss_delta_pct, max_loss, ...
    max_loading_pct, mean_system_max_loading_pct, max_loading_branch_idx, max_loading_branch, max_loading_hour, ...
    hours_any_gt80, hours_any_gt90, hours_any_gt95, hours_any_gt99, ...
    branch_hours_gt80, branch_hours_gt90, branch_hours_gt95, branch_hours_gt99, ...
    stress_above80, stress_above90, ...
    min_vm, max_vm, max_case_voltage_violation_pu, case_voltage_violation_bus_hours, ...
    screen_voltage_violation_bus_hours, mean_abs_voltage_dev_from_intact, max_abs_voltage_dev_from_intact, ...
    mean_loading_94100, max_loading_94100, mean_loading_2630, max_loading_2630, ...
    severity_score, ...
    'VariableNames', { ...
    'outage_branch_idx','outage_from_bus','outage_to_bus','outage_branch','outage_rating_MVA', ...
    'is_islanding','n_islands','smallest_island_buses','largest_island_buses', ...
    'success_hours','failed_hours','success_rate_pct','total_solve_time_s','mean_solve_time_s','max_solve_time_s', ...
    'mean_cost','mean_cost_delta','mean_cost_delta_pct','max_cost', ...
    'mean_loss_MW','mean_loss_delta_MW','mean_loss_delta_pct','max_loss_MW', ...
    'max_loading_pct','mean_hourly_system_max_loading_pct','max_loading_branch_idx','max_loading_branch','max_loading_hour', ...
    'hours_any_gt80','hours_any_gt90','hours_any_gt95','hours_any_gt99', ...
    'branch_hours_gt80','branch_hours_gt90','branch_hours_gt95','branch_hours_gt99', ...
    'stress_above80_pp_sum','stress_above90_pp_sum', ...
    'min_VM_pu','max_VM_pu','max_case_voltage_violation_pu','case_voltage_violation_bus_hours', ...
    'screen_voltage_violation_bus_hours','mean_abs_voltage_dev_from_intact_pu','max_abs_voltage_dev_from_intact_pu', ...
    'mean_loading_94_100_pct','max_loading_94_100_pct','mean_loading_26_30_pct','max_loading_26_30_pct', ...
    'screening_severity_score'});

Tsummary = sortrows(Tsummary, {'screening_severity_score','max_loading_pct'}, {'descend','descend'});
Tsummary.screening_rank = (1:height(Tsummary))';
Tsummary = movevars(Tsummary, 'screening_rank', 'Before', 'outage_branch_idx');
writetable(Tsummary, fullfile(outdir,'n1_contingency_summary_all.csv'));

%% =========================================================================
%% 8. HOURLY DIAGNOSTIC CSV
%% =========================================================================
% One row per outage-hour. For islanding cases all operating metrics are NaN.
Nrows = ncont*T;
contingency_index_col = repelem((1:ncont)', T);
outage_branch_idx_col = repelem(outage_branch_idx, T);
outage_from_col = repelem(outage_from_bus, T);
outage_to_col = repelem(outage_to_bus, T);
hour_col = repmat((1:T)', ncont, 1);
islanding_col = repelem(is_islanding, T);

Thourly = table( ...
    contingency_index_col, outage_branch_idx_col, outage_from_col, outage_to_col, ...
    hour_col, islanding_col, reshape(h_success',Nrows,1), reshape(h_solve_time',Nrows,1), ...
    reshape(h_cost',Nrows,1), reshape(h_loss',Nrows,1), reshape(h_max_loading',Nrows,1), ...
    reshape(h_min_vm',Nrows,1), reshape(h_max_vm',Nrows,1), ...
    reshape(h_loading_94100',Nrows,1), reshape(h_loading_2630',Nrows,1), ...
    'VariableNames', {'contingency_index','outage_branch_idx','outage_from_bus','outage_to_bus', ...
    'hour','is_islanding','success','solve_time_s','cost','loss_MW','max_loading_pct', ...
    'min_VM_pu','max_VM_pu','loading_94_100_pct','loading_26_30_pct'});
writetable(Thourly, fullfile(outdir,'n1_contingency_hourly_diagnostics.csv'));

%% =========================================================================
%% 9. METRIC-SPECIFIC RANKINGS AND SHORTLIST
%% =========================================================================
connected_idx = find(~is_islanding);

rank_max_loading = top_indices(max_loading_pct, connected_idx, top_per_metric);
rank_stress90 = top_indices(stress_above90, connected_idx, top_per_metric);
rank_voltage = top_indices(max_abs_voltage_dev_from_intact, connected_idx, top_per_metric);
rank_cost = top_indices(mean_cost_delta, connected_idx, top_per_metric);
rank_loss = top_indices(mean_loss_delta, connected_idx, top_per_metric);
rank_target94100 = top_indices(max_loading_94100, connected_idx, top_per_metric);
rank_target2630 = top_indices(max_loading_2630, connected_idx, top_per_metric);
rank_failures = find(~is_islanding & failed_hours > 0);
rank_islanding = find(is_islanding);

% Export individual ranking tables.
write_ranking(Tsummary, outage_branch_idx, rank_max_loading, outdir, 'n1_rank_max_loading.csv');
write_ranking(Tsummary, outage_branch_idx, rank_stress90, outdir, 'n1_rank_stress_above90.csv');
write_ranking(Tsummary, outage_branch_idx, rank_voltage, outdir, 'n1_rank_voltage_deviation.csv');
write_ranking(Tsummary, outage_branch_idx, rank_cost, outdir, 'n1_rank_cost_increase.csv');
write_ranking(Tsummary, outage_branch_idx, rank_loss, outdir, 'n1_rank_loss_increase.csv');
write_ranking(Tsummary, outage_branch_idx, rank_target94100, outdir, 'n1_rank_target_94_100.csv');
write_ranking(Tsummary, outage_branch_idx, rank_target2630, outdir, 'n1_rank_target_26_30.csv');

% Candidate union. Put infeasible connected cases first, then physical metrics.
candidate_local_idx = unique([rank_failures(:); rank_max_loading(:); rank_stress90(:); ...
    rank_voltage(:); rank_cost(:); rank_loss(:); rank_target94100(:); rank_target2630(:)], 'stable');

% Remove islanding cases from the probabilistic shortlist. They are reported
% separately because standard AC-OPF cannot represent disconnected islands
% without additional island-balancing/load-shedding modelling.
candidate_local_idx = candidate_local_idx(~is_islanding(candidate_local_idx));

% Sort candidates using the transparent screening score.
[~,ord_c] = sort(severity_score(candidate_local_idx), 'descend', 'MissingPlacement','last');
candidate_local_idx = candidate_local_idx(ord_c);
candidate_local_idx = candidate_local_idx(1:min(shortlist_max,numel(candidate_local_idx)));

% Add reason flags.
reason_max_loading = ismember((1:ncont)', rank_max_loading);
reason_stress90 = ismember((1:ncont)', rank_stress90);
reason_voltage = ismember((1:ncont)', rank_voltage);
reason_cost = ismember((1:ncont)', rank_cost);
reason_loss = ismember((1:ncont)', rank_loss);
reason_target94100 = ismember((1:ncont)', rank_target94100);
reason_target2630 = ismember((1:ncont)', rank_target2630);
reason_opf_failure = failed_hours > 0 & ~is_islanding;

Tshort = rows_by_branch_idx(Tsummary, outage_branch_idx(candidate_local_idx));
Tshort.shortlist_rank = (1:height(Tshort))';
Tshort = movevars(Tshort, 'shortlist_rank', 'Before', 'screening_rank');

% Attach flags in the same order as Tshort.
[~,loc_short] = ismember(Tshort.outage_branch_idx, outage_branch_idx);
Tshort.reason_opf_failure = reason_opf_failure(loc_short);
Tshort.reason_max_loading = reason_max_loading(loc_short);
Tshort.reason_stress90 = reason_stress90(loc_short);
Tshort.reason_voltage = reason_voltage(loc_short);
Tshort.reason_cost = reason_cost(loc_short);
Tshort.reason_loss = reason_loss(loc_short);
Tshort.reason_target_94_100 = reason_target94100(loc_short);
Tshort.reason_target_26_30 = reason_target2630(loc_short);

writetable(Tshort, fullfile(outdir,'n1_critical_contingency_shortlist.csv'));

% Separate islanding table.
Tisland = rows_by_branch_idx(Tsummary, outage_branch_idx(rank_islanding));
writetable(Tisland, fullfile(outdir,'n1_islanding_outages.csv'));

%% =========================================================================
%% 10. RUN SUMMARY AND COMPACT MAT FILE
%% =========================================================================
run_summary = table( ...
    string(case_file), nb, nbr, ng, ncont, T, peak_frac, parallel_used, ...
    nnz(is_islanding), nnz(~is_islanding & failed_hours>0), ...
    sum(success_hours), ncont*T - nnz(is_islanding)*T, ...
    base_mean_cost, base_mean_loss, base_max_loading, base_min_vm, base_max_vm, ...
    toc(t_start), ...
    'VariableNames', {'case_file','n_buses','n_branches','n_generators','n_contingencies','n_hours', ...
    'peak_frac','parallel_used','n_islanding_outages','n_connected_outages_with_OPF_failure', ...
    'successful_contingency_OPFs','attempted_connected_contingency_OPFs', ...
    'baseline_mean_cost','baseline_mean_loss_MW','baseline_max_loading_pct', ...
    'baseline_min_VM_pu','baseline_max_VM_pu','screening_wall_time_s'});
writetable(run_summary, fullfile(outdir,'n1_run_summary.csv'));

save(fullfile(outdir,'n1_screening_compact.mat'), ...
    'Tsummary','Tshort','Tisland','Tbase','run_summary', ...
    'base_success','base_cost','base_loss','base_vm','base_loading_pct', ...
    'h_success','h_cost','h_loss','h_max_loading','h_min_vm','h_max_vm', ...
    'h_loading_94100','h_loading_2630', ...
    'active_branches','is_islanding','n_islands', ...
    'target_94100','target_2630','ratings','d_t','peak_frac', ...
    'Vscreen_min','Vscreen_max','load_thresholds_pct','-v7.3');

fprintf('\n=== N-1 SCREENING COMPLETE ===\n');
fprintf('Connected outages with at least one OPF failure: %d\n', nnz(~is_islanding & failed_hours>0));
fprintf('Islanding outages: %d\n', nnz(is_islanding));
fprintf('Probabilistic shortlist size: %d\n', height(Tshort));
fprintf('Outputs: %s\n', outdir);

if ~isempty(Tshort)
    fprintf('\nTop shortlisted contingencies:\n');
    disp(Tshort(:, {'shortlist_rank','outage_branch_idx','outage_branch', ...
        'success_rate_pct','max_loading_pct','stress_above90_pp_sum', ...
        'max_abs_voltage_dev_from_intact_pu','mean_cost_delta','mean_loss_delta_MW'}));
end

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================
function mt = apply_hourly_load(mpc, Pd_frac, Qd_base, Ppk_scaled, Ppk, dt, PD, QD)
    mt = mpc;
    Plt = Ppk_scaled * dt;
    mt.bus(:,PD) = Pd_frac * Plt;
    mt.bus(:,QD) = Qd_base * (Plt / Ppk);
end

function loading_pct = branch_loading_pct(branch, ratings, PF, QF, PT, QT)
    Sf = hypot(branch(:,PF), branch(:,QF));
    St = hypot(branch(:,PT), branch(:,QT));
    loading_pct = 100 * max(Sf, St) ./ ratings;
end

function idx = find_branch_idx(mpc, fbus, tbus)
    idx = find((mpc.branch(:,1)==fbus & mpc.branch(:,2)==tbus) | ...
               (mpc.branch(:,1)==tbus & mpc.branch(:,2)==fbus), 1);
    if isempty(idx)
        error('Target branch %d-%d was not found in the case.', fbus, tbus);
    end
end

function [is_island, ncomp, smallest_n, largest_n] = topology_islanding(mpc, BUS_I, F_BUS, T_BUS, BR_STATUS)
    bus_numbers = mpc.bus(:,BUS_I);
    active = mpc.branch(:,BR_STATUS) > 0;
    f = mpc.branch(active,F_BUS);
    t = mpc.branch(active,T_BUS);

    % Convert external bus numbers to row indices for graph construction.
    [tf_f, fi] = ismember(f, bus_numbers);
    [tf_t, ti] = ismember(t, bus_numbers);
    assert(all(tf_f & tf_t), 'Branch references unknown bus numbers.');

    G = graph(fi, ti, [], numel(bus_numbers));
    bins = conncomp(G);
    counts = accumarray(bins(:), 1);
    ncomp = numel(counts);
    is_island = ncomp > 1;
    smallest_n = min(counts);
    largest_n = max(counts);
end

function out = evaluate_contingency(mpc0, outage_br, T, d_t, Pd_frac, Qd_base, ...
        Ppk_scaled, Ppk, ratings, base_vm, target_94100, target_2630, ...
        Vscreen_min, Vscreen_max, mpopt, ...
        PD, QD, VM, VMAX, VMIN, BR_STATUS, PF, QF, PT, QT)

    nb = size(mpc0.bus,1);
    nbr = size(mpc0.branch,1);

    hs = false(1,T);
    htime = nan(1,T);
    hcost = nan(1,T);
    hloss = nan(1,T);
    hmaxload = nan(1,T);
    hminvm = nan(1,T);
    hmaxvm = nan(1,T);
    h94100 = nan(1,T);
    h2630 = nan(1,T);

    load_all = nan(T,nbr);
    vm_all = nan(T,nb);
    case_vviol = nan(T,nb);

    mpc_out = mpc0;
    mpc_out.branch(outage_br,BR_STATUS) = 0;

    for t = 1:T
        mt = apply_hourly_load(mpc_out, Pd_frac, Qd_base, Ppk_scaled, Ppk, d_t(t), PD, QD);
        tt = tic;
        r = runopf(mt,mpopt);
        htime(t) = toc(tt);

        if r.success
            hs(t) = true;
            hcost(t) = r.f;
            hloss(t) = sum(r.branch(:,PF) + r.branch(:,PT));
            vm_all(t,:) = r.bus(:,VM)';

            lpct = branch_loading_pct(r.branch, ratings, PF, QF, PT, QT);
            lpct(outage_br) = NaN;
            load_all(t,:) = lpct';
            hmaxload(t) = max(lpct,[],'omitnan');
            hminvm(t) = min(r.bus(:,VM));
            hmaxvm(t) = max(r.bus(:,VM));

            if outage_br ~= target_94100
                h94100(t) = lpct(target_94100);
            end
            if outage_br ~= target_2630
                h2630(t) = lpct(target_2630);
            end

            below = max(r.bus(:,VMIN) - r.bus(:,VM), 0);
            above = max(r.bus(:,VM) - r.bus(:,VMAX), 0);
            case_vviol(t,:) = max(below,above)';
        end
    end

    valid = hs(:);
    out.h_success = hs;
    out.h_solve_time = htime;
    out.h_cost = hcost;
    out.h_loss = hloss;
    out.h_max_loading = hmaxload;
    out.h_min_vm = hminvm;
    out.h_max_vm = hmaxvm;
    out.h_loading_94100 = h94100;
    out.h_loading_2630 = h2630;

    out.success_hours = nnz(valid);
    out.failed_hours = T - out.success_hours;
    out.total_solve_time = sum(htime,'omitnan');
    out.mean_solve_time = mean(htime,'omitnan');
    out.max_solve_time = max(htime,[],'omitnan');
    out.mean_cost = mean(hcost,'omitnan');
    out.max_cost = max(hcost,[],'omitnan');
    out.mean_loss = mean(hloss,'omitnan');
    out.max_loss = max(hloss,[],'omitnan');

    if any(valid)
        [out.max_loading_pct, linear_idx] = max(load_all(:),[],'omitnan');
        [out.max_loading_hour, out.max_loading_branch_idx] = ind2sub(size(load_all), linear_idx);
        out.mean_system_max_loading_pct = mean(hmaxload,'omitnan');

        out.hours_any_gt80 = nnz(any(load_all > 80,2));
        out.hours_any_gt90 = nnz(any(load_all > 90,2));
        out.hours_any_gt95 = nnz(any(load_all > 95,2));
        out.hours_any_gt99 = nnz(any(load_all > 99,2));
        out.branch_hours_gt80 = nnz(load_all > 80);
        out.branch_hours_gt90 = nnz(load_all > 90);
        out.branch_hours_gt95 = nnz(load_all > 95);
        out.branch_hours_gt99 = nnz(load_all > 99);
        out.stress_above80 = sum(max(load_all(:)-80,0),'omitnan');
        out.stress_above90 = sum(max(load_all(:)-90,0),'omitnan');

        out.min_vm = min(vm_all(:),[],'omitnan');
        out.max_vm = max(vm_all(:),[],'omitnan');
        out.max_case_voltage_violation_pu = max(case_vviol(:),[],'omitnan');
        out.case_voltage_violation_bus_hours = nnz(case_vviol > 1e-6);
        out.screen_voltage_violation_bus_hours = nnz(vm_all < Vscreen_min | vm_all > Vscreen_max);

        vdev = abs(vm_all - base_vm);
        out.mean_abs_voltage_dev_from_intact = mean(vdev(:),'omitnan');
        out.max_abs_voltage_dev_from_intact = max(vdev(:),[],'omitnan');

        out.mean_loading_94100 = mean(h94100,'omitnan');
        out.max_loading_94100 = max(h94100,[],'omitnan');
        out.mean_loading_2630 = mean(h2630,'omitnan');
        out.max_loading_2630 = max(h2630,[],'omitnan');
    else
        out.max_loading_pct = NaN;
        out.mean_system_max_loading_pct = NaN;
        out.max_loading_branch_idx = NaN;
        out.max_loading_hour = NaN;
        out.hours_any_gt80 = 0;
        out.hours_any_gt90 = 0;
        out.hours_any_gt95 = 0;
        out.hours_any_gt99 = 0;
        out.branch_hours_gt80 = 0;
        out.branch_hours_gt90 = 0;
        out.branch_hours_gt95 = 0;
        out.branch_hours_gt99 = 0;
        out.stress_above80 = NaN;
        out.stress_above90 = NaN;
        out.min_vm = NaN;
        out.max_vm = NaN;
        out.max_case_voltage_violation_pu = NaN;
        out.case_voltage_violation_bus_hours = 0;
        out.screen_voltage_violation_bus_hours = 0;
        out.mean_abs_voltage_dev_from_intact = NaN;
        out.max_abs_voltage_dev_from_intact = NaN;
        out.mean_loading_94100 = NaN;
        out.max_loading_94100 = NaN;
        out.mean_loading_2630 = NaN;
        out.max_loading_2630 = NaN;
    end
end

function y = minmax01(x, mask)
    y = zeros(size(x));
    valid = mask & isfinite(x);
    if ~any(valid), return; end
    xmin = min(x(valid));
    xmax = max(x(valid));
    if xmax > xmin
        y(valid) = (x(valid)-xmin)/(xmax-xmin);
    end
end

function idx = top_indices(metric, eligible_idx, k)
    vals = metric(eligible_idx);
    valid = isfinite(vals);
    eligible_valid = eligible_idx(valid);
    vals = vals(valid);
    if isempty(vals)
        idx = zeros(0,1);
        return;
    end
    [~,ord] = sort(vals,'descend');
    idx = eligible_valid(ord(1:min(k,numel(ord))));
end

function Tout = rows_by_branch_idx(Tsummary, branch_indices)
    if isempty(branch_indices)
        Tout = Tsummary([],:);
        return;
    end
    [tf,loc] = ismember(branch_indices, Tsummary.outage_branch_idx);
    loc = loc(tf);
    Tout = Tsummary(loc,:);
end

function write_ranking(Tsummary, active_branch_indices, local_idx, outdir, filename)
    if isempty(local_idx)
        Tout = Tsummary([],:);
    else
        Tout = rows_by_branch_idx(Tsummary, active_branch_indices(local_idx));
        Tout.metric_rank = (1:height(Tout))';
        Tout = movevars(Tout, 'metric_rank', 'Before', 'screening_rank');
    end
    writetable(Tout, fullfile(outdir,filename));
end
