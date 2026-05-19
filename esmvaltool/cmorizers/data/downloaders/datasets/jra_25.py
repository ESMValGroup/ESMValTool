"""Script to download JRA-25 from RDA."""

import logging
import os
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
    """Download dataset.

    Parameters
    ----------
    original_data_dir : Path
        Directory where original data will be stored.
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
    downloader = WGetDownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    os.makedirs(downloader.local_folder, exist_ok=True)

    if start_date is None:
        start_date = datetime(1979, 1, 1)
    if end_date is None:
        end_date = datetime(2007, 12, 31)
    loop_date = start_date

    # download files

    url = "https://data.rda.ucar.edu/ds625.1"
    download_options = ["--no-check-certificate"]

    # define files to download

    files = ["fcst_phy2m", "fcst_phy3m", "anl_chipsi", "anl_p"]

    # download data

    while loop_date <= end_date:
        year = loop_date.year
        month = f"{loop_date.month:0>2}"

        for basename in files:
            fname = f"{basename}.{year}{month}.nc"
            # download file
            downloader.download_file(
                url + "/" + basename + "/" + fname,
                download_options,
            )

        loop_date += relativedelta.relativedelta(months=1)
