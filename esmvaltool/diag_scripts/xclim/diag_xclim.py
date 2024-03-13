# Diagnostic to caulculate climate indices using xclim
# in xclim indicators are high level functions that can be used to calculate
# a variety of climate indices. This diagnostic is a wrapper around xclim
# 
# indicators: 
#   name: potential_evapotranspiration
#   realm: atmos
#   input: ['tasmax', 'tasmin', 'pr']
#   kwargs: {'method': 'hargreaves'}


import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared import (
    # ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    select_metadata,
)
import logging

import xarray
import xclim
log = logging.getLogger()


def prepare_input_data(cfg, metas):
    # create a dataset from multiple netcdf files
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
            # parallel=True,
            engine="netcdf4",
            chunks={"time": 12},
        )[var]
    return var_data


def calculate_indicator(cfg, metas, **kwargs):
    """Calculate the index using xclim."""
    # Get the function
    ind_module = getattr(xclim.indicators, cfg.get("realm", "atmos"))
    ind_func = getattr(ind_module, cfg["indicator"])
    input_kwargs = prepare_input_data(cfg, metas)
    input_kwargs.update(cfg.get("kwargs", {}))
    indicator = ind_func(**input_kwargs)
    print(indicator)
    # Calculate the index


def main(cfg):
    grouped_meta = group_metadata(cfg["input_data"].values(), "dataset")
    for dataset, metas in grouped_meta.items():
        log.info("Calulating xclim indicator for dataset: %s", dataset)
        calculate_indicator(cfg, metas)
    log.info("Done!")


if __name__ == "__main__":
    with e.run_diagnostic() as config:
        main(config)