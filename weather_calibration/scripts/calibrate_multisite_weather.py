"""
calibrate_multisite_weather.py
================================
Builds a six-zone, spatially dependent renewable-weather model for the
IEEE 118-bus study. The zones are generic transmission regions; Portuguese
coordinates are used only as representative meteorological sampling points.

The script:
  1. downloads aligned PVGIS solar data for six representative locations;
  2. downloads ERA5 100-m wind components over a Portugal bounding box;
  3. aligns all sites on a common hourly UTC index;
  4. derives solar clear-sky indices and wind speeds;
  5. estimates site-specific marginals and temporal persistence;
  6. estimates spatial Gaussian-copula correlation matrices;
  7. exports CSV and MAT files for MATLAB trajectory generation.

Recommended common calibration period: 2019-2020, because it is covered by
both the selected PVGIS SARAH-2 dataset and ERA5.

Prerequisites
-------------
    pip install requests cdsapi xarray netCDF4 scipy numpy pandas matplotlib

CDS setup
---------
Create ~/.cdsapirc after accepting the ERA5 licence at the Copernicus CDS.

Outputs
-------
multisite_weather/
    zone_metadata.csv
    aligned_hourly_weather.csv
    solar_hourly_fitted_params.csv
    wind_site_fitted_params.csv
    solar_spatial_corr_gaussian.csv
    wind_spatial_corr_gaussian.csv
    solar_wind_cross_corr_gaussian.csv
    solar_temporal_corr.csv
    wind_temporal_corr.csv
    multisite_weather_model.mat
    multisite_weather_diagnostics.pdf

Notes
-----
* Solar dependence is estimated from daytime clear-sky indices.
* Wind dependence is estimated from hourly 100-m wind speeds.
* Rank-to-Gaussian transforms are used before computing spatial correlation,
  making the matrices suitable for Gaussian-copula simulation in MATLAB.
* Nearest-grid ERA5 extraction is used for transparency and reproducibility.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import requests
from scipy import stats
from scipy.io import savemat

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# CONFIGURATION
# =============================================================================

STARTYEAR = 2019
ENDYEAR = 2020
PVGIS_DATABASE = "PVGIS-SARAH2"
OUTDIR = Path("multisite_weather")
ERA5_FILE = OUTDIR / f"era5_multisite_{STARTYEAR}_{ENDYEAR}.nc"

MIN_GHI_THRESHOLD = 20.0       # W/m^2
CLEARSKY_PERCENTILE = 95.0
MIN_VALID_SOLAR_SAMPLES = 30
SOLAR_CLIP = (1e-4, 1.0 - 1e-4)
WIND_MIN_FIT_SPEED = 0.5       # m/s
MAX_LAG_HOURS = 24

# Generic network zones represented by geographically separated Portuguese sites.
# The labels deliberately avoid implying that IEEE-118 buses have Portuguese geography.
ZONES = [
    {"zone": "Z1_NW", "description": "Generic north-west zone",     "lat": 41.15, "lon": -8.60},
    {"zone": "Z2_NE", "description": "Generic north-east zone",     "lat": 41.50, "lon": -6.75},
    {"zone": "Z3_CW", "description": "Generic central-west zone",   "lat": 40.20, "lon": -8.40},
    {"zone": "Z4_CI", "description": "Generic central-interior zone", "lat": 39.75, "lon": -7.50},
    {"zone": "Z5_SW", "description": "Generic south-west zone",    "lat": 37.95, "lon": -8.65},
    {"zone": "Z6_SE", "description": "Generic south-east zone",    "lat": 38.05, "lon": -7.40},
]


# =============================================================================
# UTILITIES
# =============================================================================

def ensure_output_dir() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)


def zone_names() -> List[str]:
    return [z["zone"] for z in ZONES]


def gaussianize(series: pd.Series) -> pd.Series:
    """Rank-transform a series to approximately standard normal scores."""
    x = series.astype(float)
    valid = x.notna()
    out = pd.Series(np.nan, index=x.index, dtype=float)
    n = int(valid.sum())
    if n < 3:
        return out
    ranks = stats.rankdata(x.loc[valid].to_numpy(), method="average")
    u = (ranks - 0.5) / n
    u = np.clip(u, 1e-6, 1 - 1e-6)
    out.loc[valid] = stats.norm.ppf(u)
    return out


def nearest_psd_correlation(a: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Project a symmetric matrix to a numerically positive-semidefinite correlation matrix."""
    a = np.asarray(a, dtype=float)
    a = 0.5 * (a + a.T)
    vals, vecs = np.linalg.eigh(a)
    vals = np.maximum(vals, eps)
    b = (vecs * vals) @ vecs.T
    d = np.sqrt(np.diag(b))
    b = b / np.outer(d, d)
    b = np.clip(b, -1.0, 1.0)
    np.fill_diagonal(b, 1.0)
    return b


def lag_corr(x: pd.Series, lag: int) -> float:
    a = x.astype(float)
    return float(a.corr(a.shift(lag)))


def safe_ks_beta(x: np.ndarray, alpha: float, beta: float) -> Tuple[float, float]:
    if len(x) < 3:
        return np.nan, np.nan
    return stats.kstest(x, lambda q: stats.beta.cdf(q, alpha, beta))


def safe_ks_weibull(x: np.ndarray, shape: float, scale: float) -> Tuple[float, float]:
    if len(x) < 3:
        return np.nan, np.nan
    return stats.kstest(x, lambda q: stats.weibull_min.cdf(q, shape, 0, scale))


# =============================================================================
# PVGIS SOLAR
# =============================================================================

def download_pvgis_zone(zone: Dict[str, object]) -> pd.DataFrame:
    zone_id = str(zone["zone"])
    raw_path = OUTDIR / f"pvgis_{zone_id}_{STARTYEAR}_{ENDYEAR}.csv"

    if raw_path.exists():
        print(f"[PVGIS] Using cached {raw_path}")
        df = pd.read_csv(raw_path, parse_dates=["time"])
        return df

    url = "https://re.jrc.ec.europa.eu/api/v5_2/seriescalc"
    params = {
        "lat": float(zone["lat"]),
        "lon": float(zone["lon"]),
        "startyear": STARTYEAR,
        "endyear": ENDYEAR,
        "raddatabase": PVGIS_DATABASE,
        "pvcalculation": 0,
        "components": 1,
        "outputformat": "json",
        "browser": 0,
    }

    print(f"[PVGIS] Downloading {zone_id} ({zone['lat']}, {zone['lon']})")
    response = requests.get(url, params=params, timeout=180)
    if response.status_code != 200:
        raise RuntimeError(
            f"PVGIS failed for {zone_id}: HTTP {response.status_code}\n{response.text[:500]}"
        )

    payload = response.json()
    try:
        hourly = payload["outputs"]["hourly"]
    except KeyError as exc:
        raise RuntimeError(f"Unexpected PVGIS response for {zone_id}: {json.dumps(payload)[:500]}") from exc

    df = pd.DataFrame(hourly)
    df.columns = [c.strip() for c in df.columns]
    df["time"] = pd.to_datetime(df["time"], format="%Y%m%d:%H%M", utc=True)

    # PVGIS has used more than one naming convention across API versions.
    # Typical variants are Gb(i)/Gd(i)/Gr(i), Gb_i/Gd_i/Gr_i, or a direct
    # global-irradiance field G(i). Normalise these names before processing.
    alias_map = {
        "Gb(i)": "Gb_i", "Gd(i)": "Gd_i", "Gr(i)": "Gr_i", "G(i)": "G_i",
        "Gb_i": "Gb_i", "Gd_i": "Gd_i", "Gr_i": "Gr_i", "G_i": "G_i",
    }
    rename_map = {c: alias_map[c] for c in df.columns if c in alias_map}
    df = df.rename(columns=rename_map)

    # Prefer a direct global-irradiance column when supplied by PVGIS.
    if "G_i" in df.columns:
        df["G_i"] = pd.to_numeric(df["G_i"], errors="coerce")
        df["GHI_Wm2"] = df["G_i"]
    else:
        components = ["Gb_i", "Gd_i", "Gr_i"]
        missing = [c for c in components if c not in df.columns]
        if missing:
            raise RuntimeError(
                f"PVGIS response for {zone_id} lacks usable irradiance columns. "
                f"Missing {missing}; available columns are: {list(df.columns)}"
            )
        for col in components:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df["GHI_Wm2"] = df["Gb_i"] + df["Gd_i"] + df["Gr_i"]
    df[["time", "GHI_Wm2"]].to_csv(raw_path, index=False)
    print(f"[PVGIS] Saved {len(df):,} rows -> {raw_path}")
    return df[["time", "GHI_Wm2"]]


def process_solar(all_solar: Dict[str, pd.DataFrame]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return aligned solar data and hour/site Beta-parameter table."""
    frames = []
    params_rows = []

    for zone in ZONES:
        zid = str(zone["zone"])
        df = all_solar[zid].copy()
        df["time"] = pd.to_datetime(df["time"], utc=True)

        # PVGIS SARAH timestamps can be reported at a fixed minute offset
        # (commonly HH:10), whereas ERA5 uses exact top-of-hour timestamps.
        # Both records represent the same hourly interval, so normalise PVGIS
        # timestamps to the nearest UTC hour before cross-dataset alignment.
        original_minutes = sorted(df["time"].dt.minute.dropna().unique().tolist())
        if original_minutes != [0]:
            print(f"[TIME] {zid}: PVGIS minute offsets {original_minutes}; rounding to nearest hour")
        df["time"] = df["time"].dt.round("h")

        df = df.drop_duplicates("time").set_index("time").sort_index()
        df["hour"] = df.index.hour

        envelope = np.zeros(24)
        for hour in range(24):
            vals = df.loc[df["hour"] == hour, "GHI_Wm2"].dropna().to_numpy()
            vals = vals[vals > MIN_GHI_THRESHOLD]
            if len(vals) >= 10:
                envelope[hour] = np.percentile(vals, CLEARSKY_PERCENTILE)

        env_lookup = pd.Series(envelope, index=np.arange(24))
        df["solar_envelope_Wm2"] = df["hour"].map(env_lookup)
        df["solar_kt"] = df["GHI_Wm2"] / df["solar_envelope_Wm2"]
        df.loc[df["solar_envelope_Wm2"] < MIN_GHI_THRESHOLD, "solar_kt"] = np.nan
        df.loc[df["GHI_Wm2"] <= MIN_GHI_THRESHOLD, "solar_kt"] = np.nan
        df["solar_kt"] = df["solar_kt"].clip(*SOLAR_CLIP)

        for hour in range(24):
            x = df.loc[df["hour"] == hour, "solar_kt"].dropna().to_numpy()
            row = {
                "zone": zid,
                "hour_utc": hour,
                "envelope_Wm2": envelope[hour],
                "n_samples": len(x),
                "alpha": np.nan,
                "beta": np.nan,
                "mean_kt": np.nan,
                "std_kt": np.nan,
                "KS_stat": np.nan,
                "KS_pvalue": np.nan,
            }
            if len(x) >= MIN_VALID_SOLAR_SAMPLES:
                alpha, beta, _, _ = stats.beta.fit(x, floc=0, fscale=1)
                ks_stat, ks_p = safe_ks_beta(x, alpha, beta)
                row.update({
                    "alpha": alpha,
                    "beta": beta,
                    "mean_kt": float(np.mean(x)),
                    "std_kt": float(np.std(x, ddof=1)),
                    "KS_stat": ks_stat,
                    "KS_pvalue": ks_p,
                })
            params_rows.append(row)

        frames.append(
            df[["GHI_Wm2", "solar_envelope_Wm2", "solar_kt"]].rename(
                columns={
                    "GHI_Wm2": f"solar_ghi_{zid}",
                    "solar_envelope_Wm2": f"solar_env_{zid}",
                    "solar_kt": f"solar_kt_{zid}",
                }
            )
        )

    aligned = pd.concat(frames, axis=1, join="inner").sort_index()
    params = pd.DataFrame(params_rows)
    return aligned, params


# =============================================================================
# ERA5 WIND
# =============================================================================

def download_era5() -> None:
    """Download ERA5 in monthly chunks and merge them into one NetCDF file.

    The CDS cost-limit is based on the number of requested fields/grid points.
    Requesting two complete years over the full bounding box in one operation can
    exceed that limit. Monthly requests are small, resumable, and cached.
    """
    if ERA5_FILE.exists():
        print(f"[ERA5] Using cached {ERA5_FILE}")
        return

    try:
        import cdsapi
        import xarray as xr
    except ImportError as exc:
        raise RuntimeError("Install cdsapi and xarray before downloading ERA5") from exc

    north = max(float(z["lat"]) for z in ZONES) + 0.35
    south = min(float(z["lat"]) for z in ZONES) - 0.35
    west = min(float(z["lon"]) for z in ZONES) - 0.35
    east = max(float(z["lon"]) for z in ZONES) + 0.35

    chunk_dir = OUTDIR / "era5_monthly"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    client = cdsapi.Client()
    monthly_files = []

    print(
        f"[ERA5] Monthly download, bounding box "
        f"N={north}, W={west}, S={south}, E={east}"
    )

    for year in range(STARTYEAR, ENDYEAR + 1):
        for month in range(1, 13):
            month_file = chunk_dir / f"era5_wind_{year}_{month:02d}.nc"
            monthly_files.append(month_file)

            if month_file.exists() and month_file.stat().st_size > 10_000:
                print(f"[ERA5] Using cached {month_file.name}")
                continue

            request = {
                "product_type": ["reanalysis"],
                "variable": [
                    "100m_u_component_of_wind",
                    "100m_v_component_of_wind",
                ],
                "year": [str(year)],
                "month": [f"{month:02d}"],
                "day": [f"{d:02d}" for d in range(1, 32)],
                "time": [f"{h:02d}:00" for h in range(24)],
                "area": [north, west, south, east],
                "data_format": "netcdf",
                "download_format": "unarchived",
            }

            print(f"[ERA5] Downloading {year}-{month:02d} ...")
            try:
                client.retrieve(
                    "reanalysis-era5-single-levels",
                    request,
                    str(month_file),
                )
            except Exception as exc:
                # Remove a partial file so the next run retries cleanly.
                if month_file.exists():
                    month_file.unlink()
                raise RuntimeError(
                    f"ERA5 download failed for {year}-{month:02d}. "
                    "Previously completed months remain cached."
                ) from exc

            print(f"[ERA5] Saved {month_file.name}")

    print(f"[ERA5] Merging {len(monthly_files)} monthly files ...")
    datasets = []
    try:
        for file in monthly_files:
            ds = xr.open_dataset(file)
            # Load each small monthly file before closing it, avoiding a dask
            # dependency and too many simultaneously open Windows file handles.
            datasets.append(ds.load())
            ds.close()

        # Different CDS versions use either 'valid_time' or 'time'. Both are
        # handled later; combine_by_coords preserves whichever coordinate exists.
        merged = xr.combine_by_coords(datasets, combine_attrs="override")
        merged.to_netcdf(ERA5_FILE)
        merged.close()
    finally:
        for ds in datasets:
            try:
                ds.close()
            except Exception:
                pass

    print(f"[ERA5] Merged dataset saved -> {ERA5_FILE}")


def find_wind_variables(ds) -> Tuple[str, str]:
    candidates_u = ["u100", "100m_u_component_of_wind"]
    candidates_v = ["v100", "100m_v_component_of_wind"]
    u_var = next((name for name in candidates_u if name in ds.data_vars), None)
    v_var = next((name for name in candidates_v if name in ds.data_vars), None)

    if u_var is None:
        u_var = next((n for n in ds.data_vars if "u" in n.lower() and "100" in n.lower()), None)
    if v_var is None:
        v_var = next((n for n in ds.data_vars if "v" in n.lower() and "100" in n.lower()), None)
    if u_var is None or v_var is None:
        raise RuntimeError(f"Could not identify ERA5 u100/v100 variables: {list(ds.data_vars)}")
    return u_var, v_var


def find_coord_name(ds, options: Iterable[str]) -> str:
    for name in options:
        if name in ds.coords:
            return name
    raise RuntimeError(f"Could not identify coordinate among {list(options)}; coords={list(ds.coords)}")


def process_wind() -> Tuple[pd.DataFrame, pd.DataFrame]:
    import xarray as xr

    ds = xr.open_dataset(ERA5_FILE)
    u_var, v_var = find_wind_variables(ds)
    lat_name = find_coord_name(ds, ["latitude", "lat"])
    lon_name = find_coord_name(ds, ["longitude", "lon"])
    time_name = find_coord_name(ds, ["valid_time", "time"])

    frames = []
    rows = []

    for zone in ZONES:
        zid = str(zone["zone"])
        point = ds.sel(
            {lat_name: float(zone["lat"]), lon_name: float(zone["lon"])},
            method="nearest",
        )
        u = point[u_var].squeeze(drop=True)
        v = point[v_var].squeeze(drop=True)
        speed = np.sqrt(u ** 2 + v ** 2)

        idx = pd.to_datetime(speed[time_name].values, utc=True)
        values = np.asarray(speed.values, dtype=float).reshape(-1)
        series = pd.Series(values, index=idx, name=f"wind_speed_{zid}").sort_index()
        series = series[~series.index.duplicated(keep="first")]
        frames.append(series.to_frame())

        x = series.dropna().to_numpy()
        xfit = x[x > WIND_MIN_FIT_SPEED]
        shape, _, scale = stats.weibull_min.fit(xfit, floc=0)
        ks_stat, ks_p = safe_ks_weibull(xfit, shape, scale)
        rows.append({
            "zone": zid,
            "latitude_requested": float(zone["lat"]),
            "longitude_requested": float(zone["lon"]),
            "latitude_era5": float(point[lat_name].values),
            "longitude_era5": float(point[lon_name].values),
            "n_samples": len(x),
            "k_weibull": shape,
            "c_weibull_ms": scale,
            "mean_speed_ms": float(np.mean(x)),
            "std_speed_ms": float(np.std(x, ddof=1)),
            "rho_ar1": lag_corr(series, 1),
            "KS_stat": ks_stat,
            "KS_pvalue": ks_p,
        })

    ds.close()
    aligned = pd.concat(frames, axis=1, join="inner").sort_index()
    return aligned, pd.DataFrame(rows)


# =============================================================================
# DEPENDENCE ESTIMATION AND EXPORT
# =============================================================================

def dependence_outputs(aligned: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    zids = zone_names()
    solar_cols = [f"solar_kt_{z}" for z in zids]
    wind_cols = [f"wind_speed_{z}" for z in zids]

    solar_g = pd.DataFrame({z: gaussianize(aligned[f"solar_kt_{z}"]) for z in zids})
    wind_g = pd.DataFrame({z: gaussianize(aligned[f"wind_speed_{z}"]) for z in zids})

    solar_corr = pd.DataFrame(
        nearest_psd_correlation(solar_g.corr(min_periods=100).to_numpy()),
        index=zids, columns=zids,
    )
    wind_corr = pd.DataFrame(
        nearest_psd_correlation(wind_g.corr(min_periods=100).to_numpy()),
        index=zids, columns=zids,
    )
    cross_corr = solar_g.corrwith(wind_g, axis=0).to_frame("same_zone_solar_wind_corr")

    solar_temporal_rows = []
    wind_temporal_rows = []
    for zid in zids:
        for lag in range(1, MAX_LAG_HOURS + 1):
            solar_temporal_rows.append({
                "zone": zid,
                "lag_hours": lag,
                "correlation": lag_corr(solar_g[zid], lag),
            })
            wind_temporal_rows.append({
                "zone": zid,
                "lag_hours": lag,
                "correlation": lag_corr(wind_g[zid], lag),
            })

    return {
        "solar_corr": solar_corr,
        "wind_corr": wind_corr,
        "cross_corr": cross_corr,
        "solar_temporal": pd.DataFrame(solar_temporal_rows),
        "wind_temporal": pd.DataFrame(wind_temporal_rows),
    }


def make_diagnostics(aligned: pd.DataFrame, deps: Dict[str, pd.DataFrame]) -> None:
    zids = zone_names()
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    im0 = axes[0, 0].imshow(deps["solar_corr"].to_numpy(), vmin=-1, vmax=1, cmap="coolwarm")
    axes[0, 0].set_title("Solar Gaussian-copula spatial correlation")
    axes[0, 0].set_xticks(range(len(zids)), zids, rotation=45, ha="right")
    axes[0, 0].set_yticks(range(len(zids)), zids)
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046)

    im1 = axes[0, 1].imshow(deps["wind_corr"].to_numpy(), vmin=-1, vmax=1, cmap="coolwarm")
    axes[0, 1].set_title("Wind Gaussian-copula spatial correlation")
    axes[0, 1].set_xticks(range(len(zids)), zids, rotation=45, ha="right")
    axes[0, 1].set_yticks(range(len(zids)), zids)
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046)

    for zid in zids:
        s = deps["solar_temporal"]
        w = deps["wind_temporal"]
        axes[1, 0].plot(s.loc[s["zone"] == zid, "lag_hours"], s.loc[s["zone"] == zid, "correlation"], label=zid)
        axes[1, 1].plot(w.loc[w["zone"] == zid, "lag_hours"], w.loc[w["zone"] == zid, "correlation"], label=zid)

    axes[1, 0].set_title("Solar clear-sky-index temporal correlation")
    axes[1, 1].set_title("Wind-speed temporal correlation")
    for ax in axes[1, :]:
        ax.set_xlabel("Lag (hours)")
        ax.set_ylabel("Correlation")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)

    fig.tight_layout()
    fig.savefig(OUTDIR / "multisite_weather_diagnostics.pdf", bbox_inches="tight")
    plt.close(fig)


def export_matlab_model(
    solar_params: pd.DataFrame,
    wind_params: pd.DataFrame,
    deps: Dict[str, pd.DataFrame],
) -> None:
    zids = zone_names()
    solar_alpha = np.full((len(zids), 24), np.nan)
    solar_beta = np.full((len(zids), 24), np.nan)
    solar_env = np.zeros((len(zids), 24))

    for i, zid in enumerate(zids):
        p = solar_params[solar_params["zone"] == zid].set_index("hour_utc")
        solar_alpha[i, :] = p.reindex(range(24))["alpha"].to_numpy()
        solar_beta[i, :] = p.reindex(range(24))["beta"].to_numpy()
        solar_env[i, :] = p.reindex(range(24))["envelope_Wm2"].fillna(0).to_numpy()
        mx = solar_env[i, :].max()
        if mx > 0:
            solar_env[i, :] /= mx

    wp = wind_params.set_index("zone").reindex(zids)
    mdict = {
        "zone_names": np.array(zids, dtype=object),
        "zone_lat": np.array([z["lat"] for z in ZONES], dtype=float),
        "zone_lon": np.array([z["lon"] for z in ZONES], dtype=float),
        "solar_alpha_hourly": solar_alpha,
        "solar_beta_hourly": solar_beta,
        "solar_envelope_norm": solar_env,
        "solar_corr_gaussian": deps["solar_corr"].to_numpy(),
        "wind_k_weibull": wp["k_weibull"].to_numpy(dtype=float),
        "wind_c_weibull_ms": wp["c_weibull_ms"].to_numpy(dtype=float),
        "wind_rho_ar1": wp["rho_ar1"].to_numpy(dtype=float),
        "wind_corr_gaussian": deps["wind_corr"].to_numpy(),
        "same_zone_solar_wind_corr": deps["cross_corr"]["same_zone_solar_wind_corr"].reindex(zids).to_numpy(),
        "start_year": STARTYEAR,
        "end_year": ENDYEAR,
    }
    savemat(OUTDIR / "multisite_weather_model.mat", mdict, do_compression=True)


def main() -> None:
    ensure_output_dir()

    pd.DataFrame(ZONES).to_csv(OUTDIR / "zone_metadata.csv", index=False)

    all_solar = {str(z["zone"]): download_pvgis_zone(z) for z in ZONES}
    solar_aligned, solar_params = process_solar(all_solar)

    download_era5()
    wind_aligned, wind_params = process_wind()

    aligned = solar_aligned.join(wind_aligned, how="inner").sort_index()
    if aligned.empty:
        solar_min = solar_aligned.index.min()
        solar_max = solar_aligned.index.max()
        wind_min = wind_aligned.index.min()
        wind_max = wind_aligned.index.max()
        solar_minutes = sorted(pd.Index(solar_aligned.index.minute).unique().tolist())
        wind_minutes = sorted(pd.Index(wind_aligned.index.minute).unique().tolist())
        raise RuntimeError(
            "No common timestamps remain after PVGIS/ERA5 alignment. "
            f"PVGIS range={solar_min}..{solar_max}, minutes={solar_minutes}; "
            f"ERA5 range={wind_min}..{wind_max}, minutes={wind_minutes}."
        )

    expected_start = pd.Timestamp(f"{STARTYEAR}-01-01 00:00:00", tz="UTC")
    expected_end = pd.Timestamp(f"{ENDYEAR}-12-31 23:00:00", tz="UTC")
    aligned = aligned.loc[(aligned.index >= expected_start) & (aligned.index <= expected_end)]
    aligned.index.name = "time_utc"

    deps = dependence_outputs(aligned)

    aligned.reset_index().to_csv(OUTDIR / "aligned_hourly_weather.csv", index=False)
    solar_params.to_csv(OUTDIR / "solar_hourly_fitted_params.csv", index=False)
    wind_params.to_csv(OUTDIR / "wind_site_fitted_params.csv", index=False)
    deps["solar_corr"].to_csv(OUTDIR / "solar_spatial_corr_gaussian.csv", index_label="zone")
    deps["wind_corr"].to_csv(OUTDIR / "wind_spatial_corr_gaussian.csv", index_label="zone")
    deps["cross_corr"].to_csv(OUTDIR / "solar_wind_cross_corr_gaussian.csv", index_label="zone")
    deps["solar_temporal"].to_csv(OUTDIR / "solar_temporal_corr.csv", index=False)
    deps["wind_temporal"].to_csv(OUTDIR / "wind_temporal_corr.csv", index=False)

    export_matlab_model(solar_params, wind_params, deps)
    make_diagnostics(aligned, deps)

    print("\n=== MULTI-SITE WEATHER CALIBRATION COMPLETE ===")
    print(f"Common hourly records: {len(aligned):,}")
    print(f"Period: {aligned.index.min()} to {aligned.index.max()}")
    print(f"Outputs: {OUTDIR.resolve()}")
    print("\nSolar spatial Gaussian correlation:")
    print(deps["solar_corr"].round(3).to_string())
    print("\nWind spatial Gaussian correlation:")
    print(deps["wind_corr"].round(3).to_string())


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted by user")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
