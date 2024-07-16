"""Script to download ESACCI-CLOUD."""
import logging
import os

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
        start_date = datetime(2000, 1, 1)
    if end_date is None:
        end_date = datetime(2007, 12, 31)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    base_path = 'https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/L3U/AVHRR-PM/'

    wget_options = ['-r', '-nH', '-e', 'robots=off', '--cut-dirs=9', '--no-parent', '--reject="index.html*"']

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.month
        date = f'{year}{month:02}'
        if int(date) in range(200001, 200104):
            sat = 'AVHRR_NOAA-14/'
        elif int(date) in range(200104, 200509):
            sat = 'AVHRR_NOAA-16/'
        elif int(date) in range(200509, 200906):
            sat = 'AVHRR_NOAA-18/'
        else:
            logger.error("Number of instrument is not defined for date %s", date)
        folder = base_path + sat + f'{year}/{month:02}'
        logger.info("Download folder %s", folder) 
        downloader.download_folder(folder, wget_options)
        os.remove(os.path.join(downloader.local_folder, f'{month:02}'))
        loop_date += relativedelta.relativedelta(months=1)
