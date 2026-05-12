"""Script to download CALIPSO-ICECLOUD from its webpage."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import NASADownloader

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
    if not start_date:
        start_date = datetime(2007, 1, 1)
    if not end_date:
        end_date = datetime(2022, 12, 31)
    loop_date = start_date

    downloader = NASADownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        year = loop_date.year
        for month in range(1, 13):
            try:
                downloader.download_file(
                    "https://data.asdc.earthdata.nasa.gov/asdc-prod-protected/"
                    f"CALIPSO/CAL_LID_L3_Ice_Cloud-Standard-V2-00_V2-00/{year}/"
                    f"CAL_LID_L3_Ice_Cloud-Standard-V2-00.{year}-{month:02d}A.hdf",
                )
            except Exception:
                logger.info("no data downloaded for %d-%02d", year, month)
        loop_date += relativedelta.relativedelta(years=1)
