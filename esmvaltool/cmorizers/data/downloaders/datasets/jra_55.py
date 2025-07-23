"""Script to download JRA-55 from RDA."""

import logging
import os
from datetime import datetime

from dateutil import relativedelta

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

    user = os.environ.get("rda-user")
    if user is None:
        user = str(input("RDA user name? "))
        if user == "":
            errmsg = (
                "A RDA account is required to download JRA-55 data."
                " Please visit https://rda.ucar.edu/login/register/"
                " to create an account at the Research Data Archive"
                " (RDA) if needed."
            )
            logger.error(errmsg)
            raise ValueError

    passwd = os.environ.get("rda-passwd")
    if passwd is None:
        passwd = str(input("RDA password? "))

    if start_date is None:
        start_date = datetime(1958, 1, 1)
    if end_date is None:
        end_date = datetime(2022, 12, 31)
    loop_date = start_date

    options = [
        "-O",
        "Authentication.log",
        "--save-cookies=auth.rda_ucar_edu",
        f'--post-data="email={user}&passwd={passwd}&action=login"',
    ]

    # login to Research Data Archive (RDA)

    downloader.login("https://rda.ucar.edu/cgi-bin/login", options)

    # download files

    url = "https://data.rda.ucar.edu/ds628.1"
    download_options = ["--load-cookies=auth.rda_ucar_edu"]

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
                url + f"/{channel}/{year}/" + fname, download_options
            )
            # add file extension ".grb"
            os.rename(
                downloader.local_folder + "/" + fname,
                downloader.local_folder + "/" + fname + ".grb",
            )

        loop_date += relativedelta.relativedelta(years=1)

    # clean up temporary files

    if os.path.exists("Authentication.log"):
        os.remove("Authentication.log")
    if os.path.exists("auth.rda_ucar_edu"):
        os.remove("auth.rda_ucar_edu")
