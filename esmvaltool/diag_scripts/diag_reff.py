"""Python diagnostic for plotting geographical maps."""

import logging
from copy import deepcopy
from pathlib import Path

import xarray as xr
from esmvalcore.preprocessor import extract_region, mask_landsea, regrid

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = f"{attributes['short_name']}."

    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_types": ["map"],
        "authors": [
            "bock_lisa",
        ],
        "references": [
            "bock24acp",
        ],
        "ancestors": ancestor_files,
    }
    return record


def compute_weighted_mean(xr, xr_weights, vert_dim="plev"):
    """Compute weighted mean of xr with weights xr_weights along vertical dimension."""

    # Ensure weights are non-negative and have the same shape as xr
    weights = xr_weights.where(xr_weights >= 0, 0)

    # Compute numerator and denominator for weighted mean
    numerator = (xr * weights).sum(dim=vert_dim)
    denominator = weights.sum(dim=vert_dim)

    # Safe division: result = numerator / denominator, mask where denominator is zero
    weighted_mean = numerator / denominator
    weighted_mean = weighted_mean.where(denominator != 0)

    # copy attributes from original xr
    weighted_mean.attrs = xr.attrs

    return weighted_mean


def extract_highest_nonzero_level(xr, vert_dim="plev"):
    """Extract the highest non-zero level from a 3D cube."""

    # define non-zero mask (treat NaN as zero)
    nonzero = (~xr.isnull()) & (xr != 0)

    # reverse vertical axis so argmax finds the highest (top-most) True
    rev = xr.isel({vert_dim: slice(None, None, -1)})
    mask_rev = nonzero.isel({vert_dim: slice(None, None, -1)})

    # index of first True in reversed axis; where all False argmax returns 0 so mask later
    idx_rev = mask_rev.argmax(dim=vert_dim)
    has_any = mask_rev.any(dim=vert_dim)

    # select values at that index and mask locations with no non-zero level
    xr_highest = rev.isel({vert_dim: idx_rev})
    xr_highest = xr_highest.where(has_any)

    return xr_highest


def compute_cloud_optical_depth(cwp, re):
    """
    Compute Cloud Optical Depth (COD) from CWP and Effective Radius.

    Formula: tau = (3 * CWP) / (2 * rho_water * r_e)
    Assumes spherical liquid water droplets.

    Parameters
    ----------
    cwp : xr.DataArray
        Cloud Water Path with dimensions (time, lat, lon).
    re : xr.DataArray
        Cloud Effective Radius with matching dimensions.

    Returns
    -------
    xr.DataArray
        Cloud Optical Depth (dimensionless).
    """
    rho_water = 1000.0  # kg/m^3

    # Compute optical depth (aligns automatically via xarray)
    cod = (3.0 * cwp) / (2.0 * rho_water * re)

    return cod


def preprocess_cube(cube_1x1):
    """Preprocess cube for causal inference."""
    # cube_reg = extract_region(cube_1x1, start_latitude=-30., end_latitude=-10., start_longitude=265., end_longitude=285.)
    # print(cube_reg)
    cube_mask = mask_landsea(cube_1x1, mask_out="land")
    cube_regrid = regrid(
        cube_mask, target_grid="5x5", scheme="linear"
    )  # area_weighted')
    cube = extract_region(
        cube_regrid,
        start_latitude=-30.0,
        end_latitude=-10.0,
        start_longitude=265.0,
        end_longitude=285.0,
    )
    return cube


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)

    input_data = list(cfg["input_data"].values())

    groups = group_metadata(input_data, "variable_group", sort="dataset")

    xr_clouds = {}

    for var in groups:
        logger.info(f"Processing variable group: {var}")
        group = groups[var]

        selection = select_metadata(group, short_name=var, dataset="ICON")
        logger.info(f"Processing variable {var} in dataset ICON")
        attributes = selection[0]
        infile = attributes["filename"]
        ds = xr.open_dataset(infile)
        xr_clouds[var] = ds[var]

    # 1. Extract reff at highest non-zero level

    logger.info("Extracting reff at highest non-zero level.")

    xr_clouds["reff_highest"] = extract_highest_nonzero_level(
        xr_clouds["reffclwc"], vert_dim="plev"
    )

    # 2. Save cloud top pressure

    vert_dim = "plev"
    xr_clouds["ctp"] = xr_clouds["reff_highest"][vert_dim]

    # 3. Compute weighted mean of reff with weights clw

    logger.info("Computing weighted mean of reff with weights clw.")

    xr_clouds["reff_mean"] = compute_weighted_mean(
        xr_clouds["reffclwc"], xr_clouds["clw"], vert_dim="plev"
    )

    # 4. Compute optical depth with reff and lwp (tau = 3 * lwp / (2 * reff * rho_water))

    xr_clouds["tau"] = compute_cloud_optical_depth(
        xr_clouds["clwvi"], xr_clouds["reff_mean"]
    )

    # Preprocess and save the variables

    for var in ["clt", "clwvi", "reff_highest", "reff_mean", "tau", "ctp"]:
        logger.info(f"Saving variable {var}.")

        cube_1x1 = xr_clouds[var].to_iris()
        cube = preprocess_cube(cube_1x1)

        # Write output
        basename = attributes["dataset"] + "_" + var
        caption = attributes  #'cloud_effective_radius'
        provenance_record = get_provenance_record(
            caption, ancestor_files=[infile]
        )
        save_data(basename, provenance_record, cfg, cube)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
