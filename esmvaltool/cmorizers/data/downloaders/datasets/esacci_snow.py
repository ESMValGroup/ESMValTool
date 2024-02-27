"""Script to download ESACCI-SNOW."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader

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
        start_date_scfg = datetime(1982, 1, 1)
        start_date_swe = datetime(1979, 1, 1)
    else:
        start_date_scfg = start_date
        start_date_swe = start_date
    if end_date is None:
        end_date_scfg = datetime(2018, 12, 31)
        end_date_swe = datetime(2019, 12, 31)
    else:
        start_date_scfg = start_date
        start_date_swe = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    version = 'v2.0'

    # download snow cover on ground (scfg)
    loop_date = start_date_scfg
    rel_base_dir = f'scfg/AVHRR_MERGED/{version}'
    while loop_date <= end_date_scfg:
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}')
            if downloader.exists(f"{loop_date.month:02}"):
                downloader.download_folder(f'{loop_date.month:02}', 'scfg')
        else:
            logger.info(f'{loop_date.year}/{loop_date.month:02}:'
                        f' no data for scfg version {version}')
        loop_date += relativedelta.relativedelta(months=1)

    # download snow water equivalent (swe)
    loop_date = start_date_swe
    rel_base_dir = f'swe/MERGED/{version}'
    while loop_date <= end_date_swe:
        downloader.set_cwd(rel_base_dir)
        if downloader.exists(f'{loop_date.year}'):
            downloader.set_cwd(f'{rel_base_dir}/{loop_date.year}')
            if downloader.exists(f"{loop_date.month:02}"):
                downloader.download_folder(f'{loop_date.month:02}', 'swe')
        else:
            logger.info(f'{loop_date.year}/{loop_date.month:02}:'
                        f' no data for swe version {version}')
        loop_date += relativedelta.relativedelta(months=1)
