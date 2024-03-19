"""Script to download ESACCI-SOILMOISTURE."""
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader


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
        start_date = datetime(2003, 1, 1)
    if end_date is None:
        end_date = datetime(2007,12, 31)

    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = 'sst'
    downloader.connect()
    downloader.set_cwd('CDR_v2/Analysis/L4/v2.1/')
    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.strftime("%m")
        day = loop_date.strftime("%d")
        downloader.download_folder(f'./{year}/{month}/{day}/')
        loop_date += relativedelta.relativedelta(days=1)
