# Known limitations

- The IEEE 118-bus network is not geographically mapped to the used coordinates.
- Distributed scenarios combine electrical relocation, capacity splitting,
  and meteorological diversification.
- Most hourly solar Beta fits fail formal goodness-of-fit tests and are used
  as bounded empirical approximations.
- Solar-wind cross-correlation is not imposed in the principal generator.
- Demand uncertainty, forecast error, unit commitment, ramping, reserves,
  dynamic security, remedial actions, and load shedding are not modelled.
- The probabilistic N-1 stage evaluates six screened contingencies rather
  than all connected outages.
- Islanding outages are identified but not solved.
- Hourly AC-OPFs are independent and do not represent chronological market
  schedules with inter-temporal generator constraints.
