# Data provenance

## Solar

PVGIS-SARAH2 solar irradiance data were obtained for six representative
mainland-Portuguese locations. PVGIS timestamps with a fixed minute offset
were rounded to the nearest UTC hour before alignment.

## Wind

ERA5 100-m zonal and meridional wind components were extracted at the
nearest grid point to each zone. Wind speed is calculated as
`sqrt(u100^2 + v100^2)`.

## Calibration period

- Start: 2019-01-01 00:00 UTC
- End: 2020-12-31 23:00 UTC
- Common timestamps: 17,544 hourly records

## Network

The electrical network is the rated PGLib-OPF implementation of the IEEE
118-bus system, `case_pglib_opf_case118_ieee`.

## Interpretation

The meteorological zones are generic data-generating locations and are not
a geographic embedding of the IEEE 118-bus network.
