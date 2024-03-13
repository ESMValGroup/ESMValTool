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
)

import xarray
import xclim


def prepare_input_data(cfg, metas):
    # create a dataset from multiple netcdf files
    datasets = {}
    input_variables = cfg("input", group_metadata(metas, "variable").keys())
    for var in input_variables:
        datasets[var] = xarray.open_dataset(
            cfg["input_files"],
            parallel=True,
            engine="h5netcdf",
            chunks={"time": 12},
        )
    return datasets


def calculate_indicator(cfg, realm, indicator, **kwargs):
    """Calculate the index using xclim."""
    # Get the function
    ind_module = getattr(xclim.indicators, realm)
    ind_func = getattr(ind_module, indicator)
    input_kwargs = prepare_input_data(cfg, metas)
    input_kwargs.update(cfg.get("kwargs", {}))
    indicator = ind_func(**input_kwargs)    
    print(indicator)
    # Calculate the index


def main(config):
    indicator = "potential_evapotranspiration"
    realm = "atmos"
    calculate_indicator(config, realm, indicator, method="")


if __name__ == "__main__":
    with e.run_diagnostic() as config:
        main(config)