"""Script to download NCEP-DOE-R2."""

import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite
):
    """
    Download dataset.

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

    os.makedirs(downloader.local_folder, exist_ok=True)

    url = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/"

    downloader.download_file(
        url + "pressure/omega.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(
        url + "pressure/rhum.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(url + "pressure/air.mon.mean.nc", wget_options=[])
    downloader.download_file(
        url + "gaussian_grid/tcdc.eatm.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(
        url + "surface/pr_wtr.eatm.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(
        url + "gaussian_grid/prate.sfc.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(
        url + "gaussian_grid/uflx.sfc.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(
        url + "gaussian_grid/vflx.sfc.mon.mean.nc", wget_options=[]
    )
    downloader.download_file(
        url + "gaussian_grid/skt.sfc.mon.mean.nc", wget_options=[]
    )
