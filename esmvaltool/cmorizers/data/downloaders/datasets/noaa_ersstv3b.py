"""Script to download NOAA-ERSST-v3b."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite
):
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
        start_date = datetime(1854, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 1, 1)

    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    base_path = (
        "https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v3b/netcdf"
        "/ersst.{year}{month:02d}.nc"
    )

    while loop_date <= end_date:
        downloader.download_folder(
            base_path.format(year=loop_date.year, month=loop_date.month), []
        )
        loop_date += relativedelta.relativedelta(months=1)
