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


def main(cfg):
    """Compute the time average for each input dataset."""
    # Describe the preprocessed data that we will use as input
    input_data = cfg["input_data"].values()
    logger.info("INPUT DATA:\n%s", input_data)

    # Loop over variables/datasets in alphabetical order
    groups = group_metadata(input_data, "variable_group", sort="dataset")
    logger.info("GROUPS:\n%s", pformat(groups))
    for group_name in groups:
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes["dataset"])

            input_file = attributes["filename"]
            logger.info("Input data file: ", input_file)

            output_basename = Path(input_file).stem
            attributes["plot_dir"] = cfg["plot_dir"]
            attributes["work_dir"] = cfg["work_dir"]

            if attributes["variable_group"] == "ua":
                for pressure_level in [85000.0, 20000.0]:
                    attributes["cube"] = iris.load_cube(
                        input_file,
                        iris.Constraint(air_pressure=pressure_level),
                    )

                    if pressure_level == 85000.0:
                        attributes["varname"] = "x_wind_850hPa"
                    elif pressure_level == 20000.0:
                        attributes["varname"] = "x_wind_200hPa"

                    logger.info("Attributes:\n%s", attributes)
                    # Call Spectra calculations
                    spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
                    spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()

            elif attributes["variable_group"] == "pr":
                attributes["cube"] = iris.load_cube(input_file)
                attributes["varname"] = "Precipitation"

                logger.info("Attributes:\n%s", attributes)
                # Call Spectra calculations
                spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
                spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()

            elif attributes["variable_group"] == "rlut":  # Check rlut TODO
                attributes["cube"] = iris.load_cube(input_file)
                attributes["varname"] = "toa_outgoing_longwave_flux"

                logger.info("Attributes:\n%s", attributes)
                # Call Spectra calculations
                spectra_compute.WKSpectra(cfg, attributes).wkSpaceTime()
                spectra_compute.WKSpectra(cfg, attributes).SpectraSeason()

            if group_name != attributes["short_name"]:
                output_basename = group_name + "_" + output_basename
            if "caption" not in attributes:
                attributes["caption"] = input_file

            logger.info("output_basename: %s", output_basename)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
