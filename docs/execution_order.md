# Execution order

## 0. Prerequisites

1. Install Python 3 and packages from `environment/requirements.txt`.
2. Install MATLAB R2024b.
3. Install MATPOWER 8.1 and confirm that `runopf` is on the MATLAB path.
4. Obtain the rated PGLib IEEE 118-bus case:
   `case_pglib_opf_case118_ieee`.
5. Confirm writable output directories and sufficient disk space.
6. For stochastic simulations, start a parallel pool with six workers if
   reproducing the reported configuration.

## 1. Weather calibration

Run from `weather_calibration/scripts/`:

```bash
python calibrate_multisite_weather.py
```

Expected calibration products include:

- `zone_metadata.csv`
- `solar_hourly_fitted_params.csv`
- `wind_site_fitted_params.csv`
- `solar_spatial_corr_gaussian.csv`
- `wind_spatial_corr_gaussian.csv`
- `solar_temporal_corr.csv`
- `wind_temporal_corr.csv`
- `solar_wind_cross_corr_gaussian.csv`
- `multisite_weather_model.mat`

The final MAT file must be placed where the MATLAB scripts expect it, or the
script path must be updated explicitly.

## 2. Intact-network convergence study

Run:

```matlab
run('simulations/intact/mcs_convergence_multisite.m')
```

This is the authoritative intact-network script. It uses:

- \(N_{\max}=2000\) complete 24-hour trajectories;
- common trajectories for S1--S4;
- bootstrap seed 4201;
- maximum apparent flow at either branch end.

Expected publication CSV files are listed in
`docs/expected_outputs.md`.

## 3. Deterministic N-1 screening

Run:

```matlab
run('simulations/deterministic_n1/n1_deterministic_screening.m')
```

The script evaluates all 186 single-branch outages, identifies islanding
outages, and solves 24 no-renewable AC-OPFs for each connected outage.

## 4. Targeted probabilistic N-1 study

Run:

```matlab
run('simulations/probabilistic_n1/n1_probabilistic_multisite.m')
```

It uses:

- six selected contingencies;
- 400 complete trajectories;
- common random numbers;
- bootstrap seed 7421;
- topology-matched no-renewable baselines.

## 5. Peak-demand sensitivity

Run:

```matlab
run('simulations/peak_sensitivity/peak_demand_sensitivity_multisite.m')
```

Demand factors:

```text
0.60, 0.65, 0.70, 0.75, 0.80
```

The same 400 renewable trajectories are reused at each demand level.

## 6. Validate frozen outputs

```bash
python tools/validate_frozen_results.py
python tools/verify_checksums.py
```

## 7. Regenerate publication figures

Open `figures/source/paper_figures_final.m` and set its input directories to:

```text
results/frozen/intact
results/frozen/deterministic_n1
results/frozen/probabilistic_n1
results/frozen/peak_sensitivity
```

Then run the script in MATLAB.
