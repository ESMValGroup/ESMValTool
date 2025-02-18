"""Script to download NOAA-CIRES-20CR-V2."""

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

    downloader.set_cwd("/Projects/20thC_ReanV2/Monthlies/")
    downloader.download_file(
        "monolevel/cldwtr.eatm.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "monolevel/pr_wtr.eatm.mon.mean.nc", sub_folder="surface"
    )
    downloader.download_file(
        "pressure/shum.mon.mean.nc", sub_folder="pressure"
    )
    downloader.download_file(
        "gaussian/monolevel/tcdc.eatm.mon.mean.nc", sub_folder="surface_gauss"
    )
    downloader.download_file(
        "gaussian/monolevel/ulwrf.ntat.mon.mean.nc", sub_folder="surface_gauss"
    )
    downloader.download_file(
        "gaussian/monolevel/uswrf.ntat.mon.mean.nc", sub_folder="surface_gauss"
    )
    downloader.download_file(
        "gaussian/monolevel/prate.mon.mean.nc", sub_folder="surface_gauss"
    )
    downloader.download_file(
        "gaussian/monolevel/uflx.mon.mean.nc", sub_folder="surface_gauss"
    )
    downloader.download_file(
        "gaussian/monolevel/vflx.mon.mean.nc", sub_folder="surface_gauss"
    )
