"""Script to download CFSv2 from RDA."""
import logging
import os

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


logger = logging.getLogger(__name__)


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
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
            errmsg = ("A RDA account is required to download JRA-55 data."
                      " Please visit https://rda.ucar.edu/login/register/"
                      " to create an account at the Research Data Archive"
                      " (RDA) if needed.")
            logger.error(errmsg)
            raise ValueError

    passwd = os.environ.get("rda-passwd")
    if passwd is None:
        passwd = str(input("RDA password? "))

    if start_date is None:
        start_date = datetime(2011, 1, 1)
    if end_date is None:
        #end_date = datetime(2022, 12, 31)
        end_date = datetime(2011, 12, 31)
    loop_date = start_date
    print(loop_date)

    options = ["-O", "Authentication.log", "--save-cookies=auth.rda_ucar_edu",
               f"--post-data=\"email={user}&passwd={passwd}&action=login\""]

    # login to Research Data Archive (RDA)

    downloader.login("https://rda.ucar.edu/cgi-bin/login", options)

    url = "https://rda.ucar.edu/data/ds094.2/regular/"
    download_options = ["--load-cookies=auth.rda_ucar_edu"]

    # download data
    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.month
        
        fname = f'flxl.gdas.{year}{month:02}.tar'
        downloader.download_file(url + fname, download_options)
        fname = f'pgbh.gdas.{year}{month:02}.tar'
        downloader.download_file(url + fname, download_options)

        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)

    if os.path.exists("Authentication.log"):
        os.remove("Authentication.log")
    if os.path.exists("auth.rda_ucar_edu"):
        os.remove("auth.rda_ucar_edu")
