"""Script to download ISCCP-H (cloud parameters, 3-hourly)."""

import logging

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

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
        start_date = datetime(1984, 1, 1, hour=0)
    if end_date is None:
        end_date = datetime(2016, 12, 31, hour=21)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    rootpath = ("https://www.ncei.noaa.gov/data/international-satellite-cloud-"
                "climate-project-isccp-h-series-data/access/isccp-basic/hgg")

    loop_date = start_date
    while loop_date <= end_date:
        year = loop_date.year
        month = f'{loop_date.month:0>2}'
        day = f'{loop_date.day:0>2}'
        hour = f'{loop_date.hour:0>2}'
        folder = f'{downloader.local_folder}/{year}/'
        try:
            downloader.download_file(
                f'{rootpath}/{year}{month}/ISCCP-Basic.HGG.v01r00.GLOBAL.'
                f'{year}.{month}.{day}.{hour}00.GPC.10KM.CS00.EA1.00.nc',
                wget_options=[],
                out_folder=folder)
        except:
            logger.info(f'No data found for {year}-{month}-{day}, {hour}h')
        loop_date += relativedelta.relativedelta(hours=3)
