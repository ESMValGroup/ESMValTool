"""Diagnostic to caulculate climate indices using xclim.

Description
-----------
in xclim indicators are low level functions that can be used to calculate
a variety of climate indices. This diagnostic is a wrapper around xclim
by passing the input data and parameters to the xclim indicator function.
The output is saved to a netcdf file along with a metadata.yml to use it as
ancestor for further calculations.

Configuration options in recipe
-------------------------------
realm: str, required
    Name of the realm used as submodule in ``xclim.indicators``. This includes:
    'atmos', 'land', 'seaIce', 'ocean', 'cf', 'icclim', 'anuclim'. For a
    complete list and more details:
    https://xclim.readthedocs.io/en/stable/api_indicators.html
indicator: str, required
    The name of the xclim indicator to calculate. Must be a valid indicator
    from the xclim package in the provided realm. For a complete list and more
    details: https://xclim.readthedocs.io/en/stable/api_indicators.html
input: list, optional
    List of strings with short_names to be passed to the xclim function.
    By default all variables are passed.
kwargs: dict, optional
    Additional keyword arguments that are passed to the xclim indicator
    function. Even though this is an optional parameter for the diagnostic
    there might be additional required arguments for the indicator from xclim.
    See https://xclim.readthedocs.io/en/stable/api_indicators.html for details.
    By default: None
post_proc: dict, optional
    Dictionary with esmvalcore preprocessor configuration. Preprocessor
    functions are applied after calculating the index.
    NOTE: this is an experimental feature and is not tested for most
    preprocessor functions.
    By default: None
"""

import logging
import os

import esmvalcore.preprocessor as pp
import xarray
import xclim
import yaml

from esmvaltool.diag_scripts.shared import (  # ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

log = logging.getLogger()


def prepare_input_data(cfg, metas):
    """Load xarray datasets from netcdf files."""
    var_data = {}
    all_variables = group_metadata(metas, "short_name").keys()
    input_variables = cfg.get("input", all_variables)
    for var in input_variables:
        var_metas = select_metadata(metas, short_name=var)
        if len(var_metas) != 1:
            raise ValueError(f"Exactly one file per variable is required \
                found {len(var_metas)} for {var}")
        var_data[var] = xarray.open_dataset(
            var_metas[0]["filename"],
            engine="netcdf4",
            chunks={"time": 12},
        )[var]
    return var_data


def calculate_indicator(cfg, metas):
    """Calculate the index using xclim."""
    ind_module = getattr(xclim.indicators, cfg.get("realm", "atmos"))
    ind_func = getattr(ind_module, cfg["indicator"])
    input_kwargs = prepare_input_data(cfg, metas)
    input_kwargs.update(cfg.get("kwargs", {}))
    indicator = ind_func(**input_kwargs)
    print(indicator)
    return indicator


def save_results(cfg, meta, indicator, dataset, output):
    """Write output to netcdf file."""
    basename = f"{dataset}_{indicator.name}"
    fname = get_diagnostic_filename(basename, cfg)
    log.info("Saving results to %s", fname)
    indicator = indicator.rename({"lat": "latitude", "lon": "longitude"})
    indicator.coords["latitude"].attrs["units"] = "degrees"
    indicator.coords["longitude"].attrs["units"] = "degrees"
    if "standard_name" in indicator.attrs:
        del indicator.attrs["standard_name"]
    # indicator.to_netcdf(fname)
    cube = indicator.to_iris()
    cube = post_proc(cube, cfg.get("post_proc", {}))
    pp.save([cube], fname)
    meta = meta.copy()
    meta["filename"] = fname
    meta["short_name"] = indicator.name
    meta["long_name"] = indicator.attrs["long_name"]
    meta["standard_name"] = indicator.attrs.get("standard_name")
    output[fname] = meta


def post_proc(cube, procs):
    """Post processing of the indicator."""
    for coord in cube.coords():
        if not coord.has_bounds():
            coord.guess_bounds()
    for key, kwargs in procs.items():
        meth = getattr(pp, key)
        cube = meth(cube, **kwargs)
    return cube


def main(cfg):
    """Calculate xclim indicator per dataset."""
    grouped_meta = group_metadata(cfg["input_data"].values(), "dataset")
    output = {}
    for dataset, metas in grouped_meta.items():
        log.info("Calculating xclim indicator for dataset: %s", dataset)
        indicator = calculate_indicator(cfg, metas)
        save_results(cfg, metas[0], indicator, dataset, output)
    # write metadata.yml for output
    metadata_file = os.path.join(cfg["work_dir"], "metadata.yml")
    with open(metadata_file, "w", encoding="utf-8") as f:
        yaml.dump(output, f)
    log.info("Done!")


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
