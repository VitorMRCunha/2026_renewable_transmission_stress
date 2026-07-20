# Expected outputs and validation targets

## Definitive computational counts

| Study | Successful OPFs | Failed OPFs |
|---|---:|---:|
| Intact multi-site convergence | 192,024 | 0 |
| Deterministic connected N-1 screening | 4,248 | 0 |
| Targeted probabilistic N-1 | 230,400 | 0 |
| Peak-demand sensitivity | 192,120 | 0 |
| **Total** | **618,792** | **0** |

## Intact-network mandatory files

- `daily_metrics_S1_Low_Conc.csv`
- `daily_metrics_S2_High_Conc.csv`
- `daily_metrics_S3_Low_Dist.csv`
- `daily_metrics_S4_High_Dist.csv`
- `convergence_factorial_interaction.csv`
- `convergence_paired_scenario_differences.csv`
- `convergence_rank_stability.csv`
- `convergence_representative_branches.csv`
- `convergence_system_metrics.csv`
- `opf_convergence_summary.csv`

### Principal spot checks

Branch 26--30 LCP(80%) at \(N=2000\):

| Scenario | Expected value |
|---|---:|
| S1 | 5.517% |
| S2 | 10.610% |
| S3 | 2.529% |
| S4 | 6.298% |

Penetration-siting interaction:

```text
1.325 percentage points
95% CI: 0.519 to 2.122 percentage points
```

Branch 94--100:

```text
No-renewable LCP(80%) = 100%
S1--S4 LCP(80%) = 100%
Delta LCP(80%) = 0 percentage points
```

Mean daily losses:

| Scenario | MWh/day |
|---|---:|
| S1 | 1304.75 |
| S2 | 1364.54 |
| S3 | 1297.86 |
| S4 | 1296.90 |

Curtailment:

```text
Only S2 has non-zero mean curtailment.
Mean S2 wind curtailment ~= 0.181 MWh/day.
Trajectories with any curtailment ~= 0.95%.
```

## Deterministic N-1 mandatory files

- `n1_contingency_summary_all.csv`
- `n1_critical_contingency_shortlist.csv`
- `n1_run_summary.csv`

Checks:

```text
186 in-service branch outages assessed
9 islanding outages
177 connected outages
4,248 connected contingency OPFs
0 failed connected contingency OPFs
```

## Probabilistic N-1 mandatory files

- `n1_prob_contingency_baselines.csv`
- `n1_prob_opf_summary.csv`
- `n1_prob_paired_scenario_differences.csv`
- `n1_prob_scenario_contingency_ranking.csv`
- `n1_prob_system_metrics.csv`
- `n1_prob_target_branch_metrics.csv`

Principal S2 branch 26--30 checks:

```text
Outage 94--100:
  baseline LCP80 = 0%
  renewable LCP80 = 12.552%
  Delta LCP80 = 12.552 pp

Outage 94--96:
  baseline LCP80 = 8.333%
  renewable LCP80 = 19.958%
  Delta LCP80 = 11.625 pp

Outage 26--25:
  baseline LCP80 = 100%
  renewable LCP80 ~= 99.823%

Outage 23--25:
  baseline LCP80 = 100%
  renewable LCP80 ~= 99.750%
```

Total probabilistic N-1 OPFs:

```text
230,400 successful
0 failed
```

## Peak-demand sensitivity mandatory files

- `peak_sensitivity_branch_robustness.csv`
- `peak_sensitivity_opf_summary.csv`
- `peak_sensitivity_run_summary.csv`
- `peak_sensitivity_system_metrics.csv`
- `peak_sensitivity_target_branches.csv`

Checks:

```text
Demand factors = 0.60, 0.65, 0.70, 0.75, 0.80
Branch 94--100 baseline and renewable LCP80 = 100% at all levels
```

Branch 26--30 S2 baseline / renewable / increment:

| Peak factor | Baseline | S2 | Delta |
|---:|---:|---:|---:|
| 0.60 | 25.00 | 31.34 | 6.34 |
| 0.65 | 16.67 | 23.94 | 7.27 |
| 0.70 | 0.00 | 10.38 | 10.38 |
| 0.75 | 0.00 | 5.28 | 5.28 |
| 0.80 | 0.00 | 2.98 | 2.98 |

Total peak-sensitivity OPFs:

```text
192,120 successful
0 failed
```

## Figure outputs

- `Fig_MC_Convergence.pdf`
- `Fig_Structural_vs_RES.pdf`
- `Fig_Losses_Curtailment.pdf`
- `Fig_N1_Deterministic.pdf`
- `Fig_N1_Probabilistic.pdf`
- `Fig_Peak_Sensitivity.pdf`
- `FigS_Rank_Stability.pdf`

## Validation tolerance

Floating-point spot checks use small numerical tolerances because CSV
formatting, MATLAB versions, and bootstrap implementation details may alter
the last displayed decimal. Solver counts and branch identities must match
exactly.
