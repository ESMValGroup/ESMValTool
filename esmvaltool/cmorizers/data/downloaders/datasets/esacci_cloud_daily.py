"""Script to download ESACCI-CLOUD."""
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
        start_date = datetime(2000, 1, 1)
    if end_date is None:
        end_date = datetime(2000, 12, 31)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    base_path = 'https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/L3U/AVHRR-PM/AVHRR_NOAA-14/'

    wget_options = ['-r', '-nH', '-e', 'robots=off', '--cut-dirs=9', '--no-parent', '--reject="index.html*"']

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.month
        folder = base_path + f'{year}/{month:02}'
        print(folder)
        downloader.download_folder(folder, wget_options)
        loop_date += relativedelta.relativedelta(months=1)
