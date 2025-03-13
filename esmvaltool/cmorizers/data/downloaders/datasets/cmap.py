"""Script to download CMAP (CPC Merged Analysis of Precipitation)."""

import logging

from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader

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
    downloader = FTPDownloader(
        config=config,
        server="ftp2.psl.noaa.gov",
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    downloader.download_file("/Datasets/cmap/enh/precip.mon.mean.nc")
