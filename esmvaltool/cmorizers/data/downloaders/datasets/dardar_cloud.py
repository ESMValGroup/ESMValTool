"""Script to download DARDAR-CLOUD from ICARE data archive."""
import logging
import os

import datetime

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
    downloader = FTPDownloader(
        config=config,
        server='ftp.icare.univ-lille1.fr',
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    user = os.environ.get("ftp-user")
    if user is None:
        user = str(input("ICARE ftp user name? "))
        if user == "":
            errmsg = ("An ICARE account is required to download DARDAR-CLOUD"
                      " data. Please visit"
                      " https://www.icare.univ-lille.fr/asd-content/register"
                      " to create an account at the ICARE Data and Services "
                      " Center if needed.")
            logger.error(errmsg)
            raise ValueError

    passwd = os.environ.get("ftp-passwd")
    if passwd is None:
        passwd = str(input("ICARE ftp password? "))

    downloader.connect(user, passwd)

    if not start_date:
        start_date = datetime.datetime(2007, 1, 1)
    if not end_date:
        # end_date = datetime.datetime(2016, 12, 31)
        end_date = datetime.datetime(2007, 12, 31)

    loop_date = start_date
    while loop_date <= end_date:
        year = loop_date.year
        month = f'{loop_date.month:0>2}'
        day = f'{loop_date.day:0>2}'
        # download DARDAR files
        try:
            downloader.set_cwd("/SPACEBORNE/CLOUDSAT/DARDAR-CLOUD.v3.10")
            downloader.download_folder(
                f'{year}/{year}_{month}_{day}/',
                sub_folder=f'{year}_DARDAR',
                filter_files='DARDAR-CLOUD_.+nc')
        except Exception:
            logger.info('No DARDAR files found for date %d-%s-%s',
                        year, month, day)
        # download CLOUDSAT precipitation files
        try:
            downloader.set_cwd("/SPACEBORNE/CLOUDSAT/2C-PRECIP-COLUMN.C5")
            downloader.download_folder(
                f'{year}/{year}_{month}_{day}/',
                sub_folder=f'{year}_PRECIP-test',
                filter_files='')
        except Exception:
            logger.info('No PRECIP files found for date %d-%s-%s',
                        year, month, day)
        loop_date += relativedelta.relativedelta(days=1)
