%% threshold_classification_sensitivity.m
%
% Post-processes the definitive peak-demand robustness CSV.
% No AC-OPF simulations are performed.
%
% Inputs:
%   peak_sensitivity_branch_robustness.csv
%
% Outputs:
%   threshold_sensitivity_all_branches.csv
%   threshold_sensitivity_summary.csv
%   threshold_sensitivity_headline_branches.csv

clear;
clc;

%% Configuration

rootdir = fileparts(mfilename('fullpath'));
dir_peak=fullfile(rootdir,'results_paper','peak_sensitivity');
outdir = fullfile('results_paper','peak_sensitivity');
if ~exist(outdir,'dir'), mkdir(outdir); end

structural_thresholds_pct = [40 50 60];
renewable_thresholds_pp   = [1 2 3];

headline_pairs = [
    94 100;
    26  30
];

%% Read definitive robustness table

T = readtable(required_file(dir_peak,'peak_sensitivity_branch_robustness.csv'));

required_vars = {
    'branch_idx'
    'from_bus'
    'to_bus'
    'min_base_LCP80_pct'
    'max_base_LCP80_pct'
    'max_DeltaLCP80_pp'
};

missing_vars = setdiff(required_vars, T.Properties.VariableNames);

assert(isempty(missing_vars), ...
    'Missing required columns: %s', strjoin(missing_vars, ', '));

%% Preallocate output rows

n_branches = height(T);
n_cases = numel(structural_thresholds_pct) * ...
          numel(renewable_thresholds_pp);

branch_rows = cell(n_branches * n_cases, 13);
summary_rows = cell(n_cases, 9);

row_branch = 0;
row_summary = 0;

%% Threshold combinations

for structural_threshold = structural_thresholds_pct

    for renewable_threshold = renewable_thresholds_pp

        row_summary = row_summary + 1;

        persistent_structural = ...
            T.min_base_LCP80_pct >= structural_threshold;

        demand_dependent_structural = ...
            T.max_base_LCP80_pct >= structural_threshold & ...
            T.min_base_LCP80_pct < structural_threshold;

        renewable_sensitive = ...
            T.max_DeltaLCP80_pp >= renewable_threshold;

        for i = 1:n_branches

            row_branch = row_branch + 1;

            class_labels = strings(0,1);

            if persistent_structural(i)
                class_labels(end+1) = "persistent structural";
            elseif demand_dependent_structural(i)
                class_labels(end+1) = "demand-dependent structural";
            end

            if renewable_sensitive(i)
                class_labels(end+1) = "renewable-sensitive";
            end

            if isempty(class_labels)
                classification = "unflagged";
            else
                classification = strjoin(class_labels, " + ");
            end

            branch_rows(row_branch,:) = {
                structural_threshold, ...
                renewable_threshold, ...
                T.branch_idx(i), ...
                T.from_bus(i), ...
                T.to_bus(i), ...
                T.min_base_LCP80_pct(i), ...
                T.max_base_LCP80_pct(i), ...
                T.max_DeltaLCP80_pp(i), ...
                persistent_structural(i), ...
                demand_dependent_structural(i), ...
                renewable_sensitive(i), ...
                classification, ...
                sprintf('%d--%d', T.from_bus(i), T.to_bus(i))
            };
        end

        persistent_pairs = branch_pair_list( ...
            T, persistent_structural);

        demand_dependent_pairs = branch_pair_list( ...
            T, demand_dependent_structural);

        renewable_pairs = branch_pair_list( ...
            T, renewable_sensitive);

        summary_rows(row_summary,:) = {
            structural_threshold, ...
            renewable_threshold, ...
            sum(persistent_structural), ...
            persistent_pairs, ...
            sum(demand_dependent_structural), ...
            demand_dependent_pairs, ...
            sum(renewable_sensitive), ...
            renewable_pairs, ...
            n_branches
        };
    end
end

%% Write detailed branch-level table

Tall = cell2table(branch_rows, ...
    'VariableNames', {
        'structural_threshold_pct'
        'renewable_threshold_pp'
        'branch_idx'
        'from_bus'
        'to_bus'
        'min_base_LCP80_pct'
        'max_base_LCP80_pct'
        'max_DeltaLCP80_pp'
        'persistent_structural'
        'demand_dependent_structural'
        'renewable_sensitive'
        'classification'
        'branch'
    });

writetable(Tall, fullfile(outdir,'threshold_sensitivity_all_branches.csv'));

%% Write compact summary

Tsummary = cell2table(summary_rows, ...
    'VariableNames', {
        'structural_threshold_pct'
        'renewable_threshold_pp'
        'n_persistent_structural'
        'persistent_structural_branches'
        'n_demand_dependent_structural'
        'demand_dependent_structural_branches'
        'n_renewable_sensitive'
        'renewable_sensitive_branches'
        'n_total_branches'
    });

writetable(Tsummary, ...
    fullfile(outdir,'threshold_sensitivity_summary.csv'));

%% Extract headline branches

is_headline = false(height(Tall),1);

for p = 1:size(headline_pairs,1)

    a = headline_pairs(p,1);
    b = headline_pairs(p,2);

    is_headline = is_headline | ...
        (Tall.from_bus == a & Tall.to_bus == b) | ...
        (Tall.from_bus == b & Tall.to_bus == a);
end

Theadline = Tall(is_headline,:);

writetable(Theadline, ...
    fullfile(outdir,'threshold_sensitivity_headline_branches.csv'));

%% Display compact results

disp(Tsummary);

fprintf('\nHeadline branch classifications:\n');
disp(Theadline(:, {
    'structural_threshold_pct'
    'renewable_threshold_pp'
    'branch'
    'min_base_LCP80_pct'
    'max_base_LCP80_pct'
    'max_DeltaLCP80_pp'
    'classification'
}));

fprintf('\nFiles written:\n');
fprintf('  threshold_sensitivity_all_branches.csv\n');
fprintf('  threshold_sensitivity_summary.csv\n');
fprintf('  threshold_sensitivity_headline_branches.csv\n');

%% Local function

function txt = branch_pair_list(T, mask)

    idx = find(mask);

    if isempty(idx)
        txt = "";
        return;
    end

    labels = strings(numel(idx),1);

    for k = 1:numel(idx)
        labels(k) = sprintf('%d--%d', ...
            T.from_bus(idx(k)), ...
            T.to_bus(idx(k)));
    end

    txt = strjoin(labels, ', ');
end

function path = required_file(folder,name)
    path=fullfile(folder,name);
    if ~isfile(path), error('Required CSV not found: %s',path); end
end