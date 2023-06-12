"""Script to download GPM-IMERG data from NASA GESDISC DATA ARCHIVE."""
import logging
import os

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import NASADownloader


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
    downloader = NASADownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    os.makedirs(downloader.local_folder, exist_ok=True)

    if start_date is None:
        start_date = datetime(2001, 1, 1)
    if end_date is None:
        end_date = datetime(2001, 12, 31)
#        end_date = datetime(2020, 12, 31)
    loop_date = start_date

    # download files

    rootpath = ('https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/'
                'GPM_3IMERGHH.06')

    loop_date = start_date
    while loop_date <= end_date:
        year = loop_date.year
        month = f'{loop_date.month:0>2}'
        day = f'{loop_date.day:0>2}'
        dayofyear = f'{loop_date.timetuple().tm_yday:0>3}'
        folder = f'{downloader.local_folder}/{year}/'
        try:
            downloader.download_folder(
                f'{rootpath}/{year}/{dayofyear}/',
                out_folder=folder)
        except Exception as exc:
            logger.info(
                f'Failed to download data for {year}-{month}-{day}. '
                f'Exception: {str(exc)}.')
        loop_date += relativedelta.relativedelta(days=1)
