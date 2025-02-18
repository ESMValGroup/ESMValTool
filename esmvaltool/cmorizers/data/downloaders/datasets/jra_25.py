"""Script to download JRA-25 from ESGF."""

import logging
import os

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

    os.makedirs(downloader.local_folder, exist_ok=True)

    url = (
        "https://esgf.nccs.nasa.gov/thredds/fileServer/CREATE-IP/"
        "reanalysis/JMA/JRA-25/JRA-25/mon/atmos/"
    )

    downloader.download_file(
        url + "clt/clt_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
    downloader.download_file(
        url + "hus/hus_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
    downloader.download_file(
        url + "prw/prw_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
    downloader.download_file(
        url + "rlut/rlut_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
    downloader.download_file(
        url + "rlutcs/rlutcs_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
    downloader.download_file(
        url + "rsut/rsut_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
    downloader.download_file(
        url + "rsutcs/rsutcs_Amon_reanalysis_JRA-25_197901-201312.nc",
        wget_options=[],
    )
