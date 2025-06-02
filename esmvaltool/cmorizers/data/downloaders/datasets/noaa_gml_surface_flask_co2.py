"""Script to download NOAA GML surface CO2 flask data from its website."""

import logging

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite
):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    dataset_info : dict
         Dataset information from the datasets.yml file
    start_date : datetime
        Start of the interval to download
    end_date : datetime
        End of the interval to download
    overwrite : bool
        Overwrite already downloaded files
    """
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    path = "https://gml.noaa.gov/aftp/data/trace_gases/co2/flask/surface/"
    file = "co2_surface-flask_ccgg_text.tar.gz"
    downloader.download_file(
        path + file,
        wget_options=[],
    )
