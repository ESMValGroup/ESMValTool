"""Script to download JRA-55 from RDA."""

import logging
import os
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
    """Download dataset.

    Parameters
    ----------
    original_data_dir : Path
        Directory where original data will be stored.
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
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    os.makedirs(downloader.local_folder, exist_ok=True)

    if start_date is None:
        start_date = datetime(1958, 1, 1)
    if end_date is None:
        end_date = datetime(2023, 12, 31)
    loop_date = start_date

    # download files

    url = "https://osdf-director.osg-htc.org/ncar/gdex/d628001"
    download_options = []

    # define variables to download

    var = [
        ["011_tmp", "anl_p125"],
        ["011_tmp", "anl_surf125"],
        ["039_vvel", "anl_p125"],
        ["071_tcdc", "fcst_surf125"],
        ["054_pwat", "fcst_column125"],
        ["058_cice", "fcst_column125"],
        ["160_csusf", "fcst_phy2m125"],
        ["162_csulf", "fcst_phy2m125"],
        ["211_uswrf", "fcst_phy2m125"],
        ["212_ulwrf", "fcst_phy2m125"],
        ["227_cw", "fcst_column125"],
        ["228_clwc", "fcst_p125"],
        ["229_ciwc", "fcst_p125"],
    ]

    # download data

    while loop_date <= end_date:
        year = loop_date.year

        for item in var:
            varname = item[0]
            channel = item[1]
            fname = f"{channel}.{varname}.{year}01_{year}12"
            # download file
            downloader.download_file(
                url + f"/{channel}/{year}/" + fname,
                download_options,
            )
            # add file extension ".grb"
            os.rename(
                downloader.local_folder + "/" + fname,
                downloader.local_folder + "/" + fname + ".grb",
            )

        loop_date += relativedelta.relativedelta(years=1)
