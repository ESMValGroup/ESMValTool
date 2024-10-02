"""Script to download PERSIANN-CDR."""

import os
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader


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
        start_date = datetime(1983, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 1, 1)
    loop_date = start_date

    base_path = (
        "https://www.ncei.noaa.gov/data/precipitation-persiann/access/"
        "{year}/")
    while loop_date <= end_date:
        print(base_path.format(year=loop_date.year))
        print(base_path)
        downloader = WGetDownloader(
            config=config,
            dataset=dataset,
            dataset_info=dataset_info,
            overwrite=overwrite,
        )
        downloader.download_folder(base_path.format(year=loop_date.year), [])
        os.remove(os.path.join(downloader.local_folder, 'index.html'))
        loop_date += relativedelta.relativedelta(years=1)
