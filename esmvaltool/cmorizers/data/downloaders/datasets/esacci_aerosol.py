"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader

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
    if start_date is None:
        start_date = datetime(1997, 1, 1)
    if end_date is None:
        end_date = datetime(2011, 12, 31)
    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    version = 'v4.3'
    algorithm = 'SU'

    # download monthly data

    loop_date = start_date
    while loop_date <= end_date:
        if loop_date.year < 2003:
            instrument = 'ATSR2'
        else:
            instrument = 'AATSR'
        rel_base_dir = f'{instrument}_{algorithm}/L3/{version}/MONTHLY'
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.download_folder(f'{loop_date.year}',
                                       f'{algorithm}-{version}-monthly')
        else:
            logger.info('%d: no data found', loop_date.year)
        loop_date += relativedelta.relativedelta(years=1)

    # download daily data

    loop_date = start_date
    while loop_date <= end_date:
        if loop_date.year < 2003:
            instrument = 'ATSR2'
        else:
            instrument = 'AATSR'
        rel_base_dir = f'{instrument}_{algorithm}/L3/{version}/DAILY'
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}')
            if downloader.exists(f"{loop_date.month:02}"):
                downloader.download_folder(f'{loop_date.month:02}',
                                           f'{algorithm}-{version}-daily')
            else:
                logger.info('%d/%d: no data found', loop_date.year,
                            loop_date.month)
        else:
            logger.info('%d: no data found', loop_date.year)
        loop_date += relativedelta.relativedelta(months=1)
