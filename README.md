# Multi-site renewable transmission-stress analysis

Reproducibility package for:

**Structural and Renewable-Driven Transmission Stress Using Multi-Site Weather, Monte Carlo AC-OPF, and Targeted N-1 Screening**

## Purpose

This repository contains the definitive scripts, frozen result tables, figures,
configuration files, and documentation used to distinguish:

- persistent structural branch stress;
- renewable-sensitive branch stress;
- intact-network and post-contingency behaviour;
- concentrated and distributed, meteorologically diversified renewable siting.

The analysis uses a calibrated six-zone PVGIS/ERA5 weather model, common
24-hour renewable trajectories, nonlinear MATPOWER AC-OPF, deterministic
single-branch outage screening, targeted probabilistic N-1 analysis, and
peak-demand sensitivity.

## Publication status

This is release candidate **v1.0.0-rc1**. 

All publication results use the branch-loading definition based on
the maximum apparent-power magnitude at either branch end:

\[
\lambda_k = \frac{\max(|S_{f,k}|, |S_{t,k}|)}{S_k^{\max}}.
\]

## Repository contents

```text
environment/                 software and solver information
config/                      scenarios, seeds, checksums, metadata
weather_calibration/         Python calibration scripts and fitted parameters
simulations/                 authoritative MATLAB simulation scripts
results/frozen/              definitive publication CSV files
figures/source/              MATLAB post-processing script
figures/final/               final publication figures
docs/                        execution, validation, provenance, and limitations
tools/                       validation and checksum utilities
```

## Two reproduction modes

### 1. Verification-only mode — recommended for reviewers

This mode uses the frozen CSV files and does not rerun the 618,792 AC-OPFs.

1. Install MATLAB R2024b and MATPOWER 8.1.
2. Place the repository at a path with write premissions.
3. Open `figures/source/paper_figures_final.m`.
4. Point its four input directories to the corresponding folders under
   `results/frozen/`.
5. Run the script.
6. Compare the generated figures with `figures/final/`.
7. Run:

```bash
python tools/validate_frozen_results.py
python tools/verify_checksums.py
```

### 2. Full reproduction mode

Full reproduction requires the externally distributed PGLib IEEE 118-bus case,
the generated `weather_calibration/model/multisite_weather_model.mat`,
MATLAB/MATPOWER, and the Python packages listed in `environment/requirements.txt`.

Execute the stages in `docs/execution_order.md`. The reported computational
scope is:

| Study | Successful AC-OPFs | Failed |
|---|---:|---:|
| Intact multi-site convergence | 192,024 | 0 |
| Deterministic connected N-1 screening | 4,248 | 0 |
| Targeted probabilistic N-1 | 230,400 | 0 |
| Peak-demand sensitivity | 192,120 | 0 |
| **Total** | **618,792** | **0** |

Runtime depends strongly on CPU, solver build, filesystem, and parallel-pool
configuration. The reported stochastic runs used six MATLAB workers.

## Software environment

- MATLAB R2024b
- MATPOWER 8.1
- MIPS 1.5.2
- MP-Opt-Model 5.0
- MOST 1.3.1
- 64-bit Windows
- six parallel workers for stochastic runs

See `environment/` for the captured environment and solver settings.

## Data provenance

- Solar irradiance: PVGIS-SARAH2.
- Wind: ERA5 100-m zonal and meridional components.
- Calibration period: 2019-01-01 through 2020-12-31.
- Common hourly timestamps: 17,544.
- Weather representation: six generic mainland-Portuguese zones.
- Network: rated PGLib-OPF IEEE 118-bus implementation.

The weather zones are not a geographical mapping of the IEEE buses.

## Scenario definition

The four scenarios cross renewable level and siting pattern:

| Scenario | Renewable level | Siting |
|---|---|---|
| S1 | Low | Concentrated |
| S2 | High | Concentrated |
| S3 | Low | Distributed |
| S4 | High | Distributed |

Exact bus, zone, fraction, and installed-capacity definitions are stored in
`config/scenario_definitions.csv`.

## Randomness and pairing

- Master renewable trajectory seed: 42.
- Intact and peak-sensitivity bootstrap seed: 4201.
- Probabilistic N-1 bootstrap seed: 7421.
- Generator: MATLAB Mersenne Twister.
- Common random numbers are used for paired scenario comparisons.
- The same renewable trajectories are reused across peak-demand levels.
- Complete 24-hour trajectories are bootstrap units.

See `config/random_seeds.md`.

## Frozen results

The simulation output CSV data used in the paper are located in `results/frozen/`.
Their package checksums uniquely identify the included versions.

## Expected outputs and numerical validation

`docs/expected_outputs.md` lists solver counts, mandatory files, and numerical
spot checks. `tools/validate_frozen_results.py` checks the principal published
values against the frozen CSV files.

## External dependencies not redistributed here

The package does not currently contain:

- the PGLib `case_pglib_opf_case118_ieee` case file;
- raw PVGIS or ERA5 downloads;

See `docs/external_dependencies.md`.

## Integrity verification

SHA-256 hashes are stored in `config/file_checksums_sha256.txt`.

```bash
python tools/verify_checksums.py
```
Generate a new checksum manifest only after intentionally changing the release:

```bash
python tools/generate_checksums.py
```

## Citation

A preliminary `CITATION.cff` is included. 

## Licence

Original MATLAB and Python source code is licensed under the BSD 3-Clause
License.

Original documentation, figures, and generated result tables are licensed
under Creative Commons Attribution 4.0 International, unless otherwise stated.

Third-party software, benchmark cases, and source datasets remain subject to
their respective licences and terms. See `LICENSE` and the `LICENSES/`
directory.

## Contact

Corresponding author email is `vrc@isep.ipp.pt`.
