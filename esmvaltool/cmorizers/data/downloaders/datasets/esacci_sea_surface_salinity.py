"""Script to download ESACCI-SOS."""

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
        start_date = datetime(2010, 1, 1)
    if end_date is None:
        end_date = datetime(2019, 1, 1)
    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()
    for version in ['v01.8', 'v02.31']:
        downloader.set_cwd(f'{version}/30days')
        loop_date = start_date
        while loop_date <= end_date:
            if downloader.exists(str(loop_date.year)):
                downloader.download_year(loop_date.year)
            else:
                logger.info('Year %s not available for version %s',
                            loop_date.year, version)
            loop_date += relativedelta.relativedelta(years=1)
