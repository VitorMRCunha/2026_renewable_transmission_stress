# External dependencies and attribution

## MATPOWER

The MATLAB simulations require MATPOWER 8.1 and use `runopf`. MATPOWER is
distributed separately under its own licence. Do not redistribute it as
first-party repository code.

## PGLib-OPF

The simulations use the rated IEEE 118-bus PGLib-OPF case
`case_pglib_opf_case118_ieee`. Obtain the case from the official PGLib-OPF
Github https://github.com/power-grid-lib/pglib-opf/blob/master/pglib_opf_case118_ieee.m

## PVGIS

Solar data originate from PVGIS-SARAH2. Users performing a fresh download
are responsible for complying with PVGIS attribution and service terms.

## ERA5

Wind data originate from ERA5. Users performing a fresh retrieval are
responsible for complying with Copernicus Climate Data Store terms and
attribution requirements.

## Python packages

See `environment/requirements.txt`. Package versions should be frozen in a
final release after testing the clean environment.
