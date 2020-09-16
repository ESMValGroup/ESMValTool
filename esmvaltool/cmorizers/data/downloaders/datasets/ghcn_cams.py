"""Script to download Duveiller2018 from its webpage."""
import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, _, __, overwrite):
    """
    Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
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
        overwrite=overwrite,
    )
    downloader.tier = 2
    os.makedirs(downloader.local_folder, exist_ok=True)
    downloader.download_file(
        "ftp://ftp.cdc.noaa.gov/Datasets/ghcncams/air.mon.mean.nc",
        wget_options=[]

    )
