"""Script to download ESACCI-SST."""
# Import required python modules
import logging
import os

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader

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
    if start_date is None:
        start_date = datetime(2004, 1, 1)
    if end_date is None:
        end_date = datetime(2007, 12, 31)

    loop_date = start_date

    user = os.environ.get("ceda-user")
    if user is None:
        user = str(input("CEDA user name? "))
        if user == "":
            errmsg = ("A CEDA account is required to download CCI SST data."
                      " Please visit https://services.ceda.ac.uk/cedasite/"
                      "register/info/ to create an account at CEDA if needed.")
            logger.error(errmsg)
            raise ValueError

    passwd = os.environ.get("ceda-passwd")
    if passwd is None:
        passwd = str(input("CEDA-password? "))

    downloader = FTPDownloader(
        config=config,
        server='ftp3.ceda.ac.uk',
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
        user=user,
        passwd=passwd,
    )

    downloader.connect()
    downloader.set_cwd('neodc/eocis/data/global_and_regional/'
                       'sea_surface_temperature/CDR_v3/Analysis/L4/v3.0.1/')

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.strftime("%m")
        day = loop_date.strftime("%d")
        downloader.download_folder(f'./{year}/{month}/{day}/')
        loop_date += relativedelta.relativedelta(days=1)
