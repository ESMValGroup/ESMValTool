"""Write the global climatological seasonal radiation budget to a text file."""

import csv
import logging
import os

import iris

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

SEASONS = {0: "djf", 1: "mam", 2: "jja", 3: "son"}


def organise_seasonal_data(model_data):
    """Return the seasonal data from the cubes.

    Parameters
    ----------
    model_data : :class:`iris.cube.CubeList`
        The cubes containing seasonal data.

    Returns
    -------
    list of lists
        The seasonal data in the form ``[[<long_name + season>, value], ...]``.
    """
    seasonal_data = []
    for cube in model_data:
        long_name = cube.long_name
        for season in cube.slices_over("season_number"):
            season_name = SEASONS[season.coord("season_number").points[0]]
            value = season.data
            seasonal_data.append([f"{long_name} {season_name}", str(value)])
        average_value = cube.data.mean()
        seasonal_data.append([f'{long_name} {"ann"}', str(average_value)])
    return seasonal_data


def write_seasonal_data_output(output_dir, model_dataset, seasonal_data):
    """Write seasonal data to CSV file.

    The CSV file will have the name ``<model_dataset>_metrics.csv`` and can be
    used for the normalised metric assessment plot.

    Parameters
    ----------
    output_dir : string
        The full path to the directory in which the CSV file will be written.
    model_dataset : string
        The model name.
    seasonal_data : list of lists
        The seasonal data to write.
    """
    file_name = f"{model_dataset}_metrics.csv"
    file_path = os.path.join(output_dir, file_name)

    with open(file_path, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        for line in seasonal_data:
            csv_writer.writerow(line)


def main(config):
    """Seasonal radiation budget comparison for models defined in the
    radiation_budget recipe file.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    logger = logging.getLogger(__name__)

    input_data = config["input_data"]
    datasets = group_metadata(input_data.values(), "dataset")

    for model_dataset, group in datasets.items():
        # 'model_dataset' is the name of the model dataset.
        # 'group' is a list of dictionaries containing metadata.
        logger.info("Processing data for %s", model_dataset)
        filenames = [item["filename"] for item in group]
        model_data = iris.load(filenames)
        seasonal_data = organise_seasonal_data(model_data)

        write_seasonal_data_output(config["work_dir"], model_dataset,
                                   seasonal_data)


if __name__ == "__main__":
    with run_diagnostic() as CONFIG:
        main(CONFIG)
