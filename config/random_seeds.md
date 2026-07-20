# Random seeds and common-random-number design

The MATLAB simulations use the Mersenne Twister generator.

| Study | Master trajectory seed | Bootstrap seed | Notes |
|---|---:|---:|---|
| Intact multi-site convergence | 42 | 4201 | One master pool reused across S1--S4 |
| Deterministic N-1 screening | 42 | Not applicable | Deterministic study; seed retained for script consistency |
| Targeted probabilistic N-1 | 42 | 7421 | Common trajectories reused across scenarios and contingencies |
| Peak-demand sensitivity | 42 | 4201 | Same renewable trajectories reused across all demand factors |
| Dual-definition N=50 audit | Script-defined | Not used for publication inference | Validation-only study |

## Independence unit

A complete 24-hour trajectory is the Monte Carlo and bootstrap sampling unit.
Hours within one trajectory are not resampled independently.

## Important implementation check

Each authoritative stochastic script explicitly calls `rng(..., 'twister')`.
Do not remove or relocate those calls when reproducing the published results.
