"""Script to download ESACCI-SST."""

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
        start_date = datetime(1980, 1, 1)
    if end_date is None:
        end_date = datetime(2021, 12, 31)

    loop_date = start_date

    downloader = WGetDownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    path = (
        "https://dap.ceda.ac.uk/neodc/eocis/data/global_and_regional/"
        "sea_surface_temperature/CDR_v3/Analysis/L4/v3.0.1/"
    )

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.strftime("%m")
        day = loop_date.strftime("%d")
        folder = path + f"{year}/{month}/{day}/"
        downloader.download_folder(
            folder,
            wget_options=["-e robots=off", "--no-parent", "--accept=nc"],
        )
        loop_date += relativedelta.relativedelta(days=1)
