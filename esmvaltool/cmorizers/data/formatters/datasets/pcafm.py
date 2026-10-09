"""ESMValTool CMORizer for PC-AFM (and all CorrDiff-based) model output.

Converts raw NetCDF output from the PC-AFM generation script into
CMOR-compliant OBS6-style files for use with ESMValTool. Latitude,
longitude, and time coordinates are read from a reference high-resolution
training input file, since the generation script does not write these
coordinates correctly.

This CMORizer is invoked via the standard ESMValTool interface::

    esmvaltool data format --datasets PCAFM

Configuration is read from
``esmvaltool/cmorizers/data/cmor_config/PCAFM.yml``, which must include
the ``reference_file`` key pointing to a nextGEMS high-resolution target
file covering the same domain and time period as the model output.

The raw output file is expected under ``in_dir`` and is identified by the
``filename`` glob pattern in the YAML config.

Output directory structure
--------------------------
::

    <rootpath>/
        Tier1/
            PCAFM/
                OBS6_PCAFM_reanaly_<ens>-<region>_3hr_<var>.nc
            HIGHRES-REF-<region>/
                OBS6_HIGHRES-REF-<region>_reanaly_1_3hr_<var>.nc

Notes
-----
- The generation script stores spatial axes as (y, x) transposed relative
  to (lat, lon); a -90 degree rotation is applied to correct alignment.
- Precipitation stored as square root (``2rootpr``) or log precipitation 
  ("log(x+10e-6)) is back-transformed to  physical units (kg m-2 s-1) 
  before writing.
- Height scalar coordinates (2 m or 10 m) are added via ESMValTool
  utilities to satisfy CF conventions for near-surface variables.
"""

import glob
import logging
import os
from datetime import datetime

import iris
import netCDF4 as nc
import numpy as np
import xarray as xr
from esmvaltool.cmorizers.data import utilities
from ncdata.iris_xarray import cubes_from_xarray

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Variable metadata
# ---------------------------------------------------------------------------

# Maps CMOR variable name -> CF metadata written to output files.
VAR_METADATA = {
    "uas": {
        "standard_name": "eastward_wind",
        "long_name": "Eastward Near-Surface Wind",
        "units": "m s-1",
    },
    "vas": {
        "standard_name": "northward_wind",
        "long_name": "Northward Near-Surface Wind",
        "units": "m s-1",
    },
    "tas": {
        "standard_name": "air_temperature",
        "long_name": "Near-Surface Air Temperature",
        "units": "K",
    },
    "pr": {
        "standard_name": "precipitation_flux",
        "long_name": "Precipitation",
        "units": "kg m-2 s-1",
    },
    "huss": {
        "standard_name": "specific_humidity",
        "long_name": "Near-Surface Specific Humidity",
        "units": "1",
    },
    "ps": {
        "standard_name": "surface_air_pressure",
        "long_name": "Surface Air Pressure",
        "units": "Pa",
    },
}

# Maps raw variable names in the generation output to (cmor_name, transform).
# Variables not listed here are passed through with their original name.
TRANSFORM_MAP = {
    "2rootpr": ("pr", lambda x: np.power(np.clip(x, 0.0, None), 2)),
    "logpr":   ("pr", lambda x: np.clip(np.exp(x) - 1e-5, 0.0, None)),
    "hus_sfc": ("huss", None),
    "pres_sfc": ("ps", None),
}

# Variables requiring height scalar coordinates (CF convention).
# Note: in nextGEMS cycle 3, huss is taken at the lowest model level, not at
# 2 m; the 2 m coordinate is attached for CMOR compliance only.
HEIGHT_2M  = {"tas", "huss"}
HEIGHT_10M = {"uas", "vas"}


# ---------------------------------------------------------------------------
# Coordinate loading
# ---------------------------------------------------------------------------

def _load_reference_coords(ref_file, n_y, n_x):
    """Read lat, lon, and time coordinates from a reference training file.

    Parameters
    ----------
    ref_file : str
        Path to a nextGEMS high-resolution target NetCDF file covering the
        same spatial domain and time period as the model output.
    n_y, n_x : int
        Expected spatial dimensions of the model output. Used to validate
        that the reference grid matches.

    Returns
    -------
    lat : np.ndarray, shape (n_y,)
    lon : np.ndarray, shape (n_x,)
    time_vals : np.ndarray
    time_units : str
    time_calendar : str
    """
    with nc.Dataset(ref_file, "r") as ds:
        lat = np.array(ds.variables["lat"][:])
        lon = np.array(ds.variables["lon"][:])
        tvar = ds.variables["time"]
        time_vals = np.array(tvar[:])
        time_units = tvar.units
        time_calendar = getattr(tvar, "calendar", "standard")

    if len(lat) != n_y:
        raise ValueError(
            f"Reference file has {len(lat)} latitude points but model "
            f"output has {n_y}."
        )
    if len(lon) != n_x:
        raise ValueError(
            f"Reference file has {len(lon)} longitude points but model "
            f"output has {n_x}."
        )

    logger.info(
        "Coordinates loaded from reference file: "
        "lat [%.2f, %.2f], lon [%.2f, %.2f], %d time steps.",
        lat.min(), lat.max(), lon.min(), lon.max(), len(time_vals),
    )
    return lat, lon, time_vals, time_units, time_calendar


def _fix_time(time_vals, freq_hours):
    """Reconstruct a regularly spaced time axis from the first valid entry.

    The PC-AFM generation script does not always write time values
    correctly beyond the first step. This function rebuilds the full
    sequence at ``freq_hours`` intervals from the first value.

    Parameters
    ----------
    time_vals : np.ndarray
        Potentially corrupted time coordinate array from the reference file.
    freq_hours : int
        Output time step in hours.

    Returns
    -------
    np.ndarray
    """
    return time_vals[0] + np.arange(len(time_vals)) * freq_hours


# ---------------------------------------------------------------------------
# Variable extraction
# ---------------------------------------------------------------------------

def _extract_variable(raw_name, group_data):
    """Return (cmor_name, var_data, var_info) for a raw variable.

    Applies back-transformations for precipitation stored in transformed
    space (square root or log).
    """
    if raw_name in TRANSFORM_MAP:
        cmor_name, transform = TRANSFORM_MAP[raw_name]
        raw_data = np.array(group_data.variables[raw_name][:])
        var_data = transform(raw_data) if transform is not None else raw_data
    elif raw_name in VAR_METADATA:
        cmor_name = raw_name
        var_data = np.array(group_data.variables[raw_name][:])
    else:
        raise KeyError(f"Variable '{raw_name}' has no metadata entry.")

    return cmor_name, var_data, VAR_METADATA[cmor_name]


# ---------------------------------------------------------------------------
# File writer
# ---------------------------------------------------------------------------

def _write_cmor_file(
    filepath, var_name, var_data, lat, lon,
    time_vals, time_units, time_calendar,
    var_info, group_name, ensemble_member,
    cfg,
):
    """Write a single CMOR-compliant NetCDF file.

    Parameters
    ----------
    filepath : str
    var_name : str
        CMOR variable name.
    var_data : np.ndarray
        Shape (time, y, x). Spatial rotation is applied internally.
    lat, lon : np.ndarray
        1-D coordinate arrays from the reference file.
    time_vals : np.ndarray
    time_units, time_calendar : str
    var_info : dict
        CF metadata: standard_name, long_name, units.
    group_name : str
        Source group label ('prediction' or 'truth').
    ensemble_member : int or None
        1-based ensemble index for prediction files; None for truth.
    cfg : dict
        CMORizer configuration dict (provides institution, source, etc.).
    """
    # The generation script stores (y, x) transposed relative to (lat, lon).
    # A -90 degree rotation corrects the alignment.
    var_data = np.rot90(var_data, k=-1, axes=(1, 2))

    ds = xr.Dataset()
    ds[var_name] = xr.DataArray(
        var_data,
        dims=["time", "lon", "lat"],
        attrs={
            "standard_name": var_info["standard_name"],
            "long_name":     var_info["long_name"],
            "units":         var_info["units"],
            "coordinates":   "lon lat",
        },
    )
    ds["time"] = xr.DataArray(
        time_vals, dims=["time"],
        attrs={
            "standard_name": "time",
            "long_name":     "time",
            "axis":          "T",
            "calendar":      time_calendar,
            "units":         time_units,
        },
    )
    ds["lon"] = xr.DataArray(
        lon, dims=["lon"],
        attrs={
            "standard_name": "longitude",
            "long_name":     "longitude",
            "units":         "degrees_east",
            "axis":          "X",
        },
    )
    ds["lat"] = xr.DataArray(
        lat, dims=["lat"],
        attrs={
            "standard_name": "latitude",
            "long_name":     "latitude",
            "units":         "degrees_north",
            "axis":          "Y",
        },
    )

    ds.attrs.update({
        "Conventions":       "CF-1.8",
        "title":             f'{var_info["long_name"]} — {group_name}',
        "institution":       cfg.get("attributes", {}).get("institution", "DLR"),
        "source":            cfg.get("attributes", {}).get(
                                 "source",
                                 "PC-AFM Physics-Constrained Adaptive Flow Matching"
                             ),
        "history":           f"CMORized on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "frequency":         "3hrPt",
        "realm":             "atmos",
        "product":           "model-output",
        "variable_id":       var_name,
        "grid":              "native grid",
        "grid_label":        "gn",
        "nominal_resolution": "~6 km",
        "comment":           f"Data from {group_name} group",
    })

    if ensemble_member is not None:
        ds.attrs["realization_index"] = ensemble_member
        ds.attrs["variant_label"]     = f"r{ensemble_member}i1p1f1"
        ds.attrs["ensemble_member"]   = ensemble_member

    cube = cubes_from_xarray(ds)[0]
    if var_name in HEIGHT_2M:
        utilities.add_height2m(cube)
    elif var_name in HEIGHT_10M:
        utilities.add_height10m(cube)
    cube = utilities.fix_coords(cube)

    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    iris.save(cube, filepath, fill_value=1e20, unlimited_dimensions=["time"])
    logger.info("Written: %s", filepath)


# ---------------------------------------------------------------------------
# ESMValTool CMORizer entry point
# ---------------------------------------------------------------------------

def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """CMORize PC-AFM generation output.

    This function is called automatically by ESMValTool when running::

        esmvaltool data format --datasets PCAFM

    Parameters
    ----------
    in_dir : str
        Directory containing the raw PC-AFM output NetCDF file(s).
        The file is located using the ``filename`` glob pattern from
        ``PCAFM.yml``.
    out_dir : str
        Root output directory for CMOR-compliant files (typically the
        OBS6 rootpath configured in the ESMValTool user config).
    cfg : dict
        Contents of ``esmvaltool/cmorizers/data/cmor_config/PCAFM.yml``.
        Required keys:
          - ``filename``       : glob pattern for the raw output file
          - ``reference_file`` : path to the nextGEMS reference NetCDF
                                 from which lat/lon/time are read
          - ``method_name``    : dataset label (e.g. 'PCAFM')
          - ``region``         : region tag appended to filenames
          - ``groups``         : list of NetCDF groups to process
          - ``time_freq_hours``: output time step in hours (default 3)
          - ``fix_time``       : whether to reconstruct the time axis
                                 (default true)
    cfg_user : dict
        ESMValTool user configuration (not used directly here).
    start_date, end_date : datetime or None
        Optional date range filtering (not currently applied; all time
        steps in the file are written).
    """
    # ------------------------------------------------------------------
    # Locate the raw input file
    # ------------------------------------------------------------------
    pattern = os.path.join(in_dir, cfg["filename"])
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(
            f"No file matching '{pattern}' found in {in_dir}. "
            "Check the 'filename' key in PCAFM.yml."
        )
    if len(matches) > 1:
        logger.warning(
            "Multiple files match '%s'; using the first: %s",
            pattern, matches[0],
        )
    input_file = matches[0]
    logger.info("CMORizing: %s", input_file)

    # ------------------------------------------------------------------
    # Read configuration
    # ------------------------------------------------------------------
    ref_file       = cfg["reference_file"]
    method_name    = cfg.get("method_name", "PCAFM")
    region         = cfg.get("region", None)
    groups         = cfg.get("groups", ["prediction", "truth"])
    freq_hours     = int(cfg.get("time_freq_hours", 3))
    fix_time_flag  = bool(cfg.get("fix_time", True))

    if not os.path.isfile(ref_file):
        raise FileNotFoundError(
            f"Reference file not found: {ref_file}. "
            "Set 'reference_file' in PCAFM.yml to a valid path."
        )

    # ------------------------------------------------------------------
    # Open raw output and infer spatial dimensions
    # ------------------------------------------------------------------
    ds_raw = nc.Dataset(input_file, "r")

    first_group = ds_raw.groups[groups[0]]
    first_var   = next(iter(first_group.variables))
    raw_shape   = first_group.variables[first_var].shape
    # Prediction: (ensemble, time, y, x); truth: (time, y, x)
    n_y, n_x = raw_shape[-2], raw_shape[-1]

    # ------------------------------------------------------------------
    # Load and optionally reconstruct coordinates from the reference file
    # ------------------------------------------------------------------
    lat, lon, time_vals, time_units, time_calendar = _load_reference_coords(
        ref_file, n_y, n_x
    )
    if fix_time_flag:
        time_vals = _fix_time(time_vals, freq_hours)
        logger.info("Time axis reconstructed at %d-hour intervals.", freq_hours)

    # ------------------------------------------------------------------
    # Process each group and variable
    # ------------------------------------------------------------------
    for group_name in groups:
        if group_name not in ds_raw.groups:
            logger.warning(
                "Group '%s' not found in %s; skipping.",
                group_name, input_file,
            )
            continue
        logger.info("Processing group: %s", group_name)
        group = ds_raw.groups[group_name]

        processable = [
            v for v in group.variables
            if v in VAR_METADATA or v in TRANSFORM_MAP
        ]
        if not processable:
            logger.warning(
                "No recognised variables in group '%s'.", group_name
            )
            continue

        for raw_name in processable:
            logger.info("  Variable: %s", raw_name)
            try:
                cmor_name, var_data, var_info = _extract_variable(
                    raw_name, group
                )
            except Exception as exc:
                logger.error("  Skipping %s: %s", raw_name, exc)
                continue

            region_tag    = f"-{region}" if region else ""
            is_prediction = (group_name == "prediction" and var_data.ndim == 4)

            if is_prediction:
                n_ens = var_data.shape[0]
                dataset_label = method_name
                for ens_idx in range(n_ens):
                    filename = (
                        f"OBS6_{dataset_label}_reanaly_"
                        f"{ens_idx + 1}{region_tag}_3hr_{cmor_name}.nc"
                    )
                    filepath = os.path.join(
                        out_dir, "Tier1", dataset_label, filename
                    )
                    _write_cmor_file(
                        filepath     = filepath,
                        var_name     = cmor_name,
                        var_data     = var_data[ens_idx],
                        lat          = lat,
                        lon          = lon,
                        time_vals    = time_vals,
                        time_units   = time_units,
                        time_calendar= time_calendar,
                        var_info     = var_info,
                        group_name   = group_name,
                        ensemble_member = ens_idx + 1,
                        cfg          = cfg,
                    )
            else:
                dataset_label = f"HIGHRES-REF{region_tag}"
                filename = (
                    f"OBS6_{dataset_label}_reanaly_1_3hr_{cmor_name}.nc"
                )
                filepath = os.path.join(
                    out_dir, "Tier1", dataset_label, filename
                )
                _write_cmor_file(
                    filepath     = filepath,
                    var_name     = cmor_name,
                    var_data     = var_data,
                    lat          = lat,
                    lon          = lon,
                    time_vals    = time_vals,
                    time_units   = time_units,
                    time_calendar= time_calendar,
                    var_info     = var_info,
                    group_name   = group_name,
                    ensemble_member = None,
                    cfg          = cfg,
                )

    ds_raw.close()
    logger.info("CMORization complete. Output written to: %s", out_dir)