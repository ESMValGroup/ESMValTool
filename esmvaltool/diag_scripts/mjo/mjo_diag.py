"""Entry point for the MJO Wheeler-Kiladis diagnostics."""

import logging
from pathlib import Path
from pprint import pformat

import iris
import spectra_compute

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


_VARIABLE_GROUPS = {
    "ua": (
        (85000.0, "x_wind_850hPa"),
        (20000.0, "x_wind_200hPa"),
    ),
    "pr": ((None, "Precipitation"),),
    "rlut": ((None, "toa_outgoing_longwave_flux"),),
}


def _load_cube(input_file, variable_group, pressure_level=None):
    """Load the input cube and derive the variable name for the spectra run."""
    if variable_group == "ua":
        cube = iris.load_cube(
            input_file, iris.Constraint(air_pressure=pressure_level)
        )
    else:
        cube = iris.load_cube(input_file)

    if variable_group == "ua" and pressure_level == 85000.0:
        varname = "x_wind_850hPa"
    elif variable_group == "ua" and pressure_level == 20000.0:
        varname = "x_wind_200hPa"
    elif variable_group == "pr":
        varname = "Precipitation"
    elif variable_group == "rlut":
        varname = "toa_outgoing_longwave_flux"
    else:
        msg = f"Unsupported variable group: {variable_group}"
        raise ValueError(msg)

    return cube, varname


def _run_spectra(cfg, attributes):
    """Run the space-time and seasonal spectra diagnostics."""
    spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
    spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()


def main(cfg):
    """Run the MJO spectra diagnostics for all configured input datasets."""
    input_data = cfg["input_data"].values()
    logger.info("INPUT DATA:\n%s", input_data)

    # Loop over variables/datasets in alphabetical order
    groups = group_metadata(input_data, "variable_group", sort="dataset")
    logger.info("GROUPS:\n%s", pformat(groups))
    for group_name, group_attributes in groups.items():
        logger.info("Processing variable %s", group_name)
        for attributes in group_attributes:
            logger.info("Processing dataset %s", attributes["dataset"])

            input_file = attributes["filename"]
            logger.info("Input data file: %s", input_file)

            output_basename = Path(input_file).stem
            attributes["plot_dir"] = cfg["plot_dir"]
            attributes["work_dir"] = cfg["work_dir"]

            for pressure_level, _ in _VARIABLE_GROUPS[
                attributes["variable_group"]
            ]:
                attributes["cube"], attributes["varname"] = _load_cube(
                    input_file,
                    attributes["variable_group"],
                    pressure_level=pressure_level,
                )
                logger.info("Attributes:\n%s", attributes)
                _run_spectra(cfg, attributes)

            if group_name != attributes["short_name"]:
                output_basename = group_name + "_" + output_basename
            if "caption" not in attributes:
                attributes["caption"] = input_file

            logger.info("output_basename: %s", output_basename)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
