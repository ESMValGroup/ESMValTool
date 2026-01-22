# pylint: disable=too-many-arguments
# pylint: disable=R0917
# pylint: disable=too-many-locals
"""Script to download IAP datasets."""

import logging
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
    if start_date is None:
        start_date = datetime(year=1940, month=1, day=1)
    if end_date is None:
        end_date = datetime(year=2024, month=12, day=31)

    loop_date = start_date

    downloader = WGetDownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        print(loop_date)
        downloader.download_file(
            "http://www.ocean.iap.ac.cn/ftp/cheng/"
            "IAPv4.2_IAP_Temperature_gridded_1month_netcdf/Monthly/"
            f"IAPv4_Temp_monthly_1_6000m_year_{loop_date.year}"
            f"_month_{loop_date.month:02d}.nc",
            wget_options=[],
        )
        loop_date += relativedelta.relativedelta(months=1)
