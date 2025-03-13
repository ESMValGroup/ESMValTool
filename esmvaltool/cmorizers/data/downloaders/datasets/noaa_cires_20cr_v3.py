"""Script to download NOAA-CIRES-20CR-V3."""

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
        server="ftp.cdc.noaa.gov",
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    downloader.set_cwd("Datasets/20thC_ReanV3/Monthlies/")
    downloader.download_file(
        "miscSI-MO/cldwtr.eatm.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "miscSI-MO/pr_wtr.eatm.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "prsSI-MO/shum.mon.mean.nc", sub_folder="pressure"
    )
    downloader.download_file(
        "miscMO/tcdc.eatm.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "ntatFlxSI-MO/ulwrf.ntat.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "ntatFlxSI-MO/uswrf.ntat.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "ntatFlxSI-MO/csulf.ntat.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "ntatFlxSI-MO/csusf.ntat.mon.mean.nc", sub_folder="surface"
    )
