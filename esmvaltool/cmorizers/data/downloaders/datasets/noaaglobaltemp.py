"""Script to download NOAAGlobalTemp from its webpage."""

import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import read_cmor_config

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

    filename = read_cmor_config(dataset)["filenames"]["gridded"]
    os.makedirs(downloader.local_folder, exist_ok=True)
    downloader.download_file(
        f"https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/"
        f"v5/access/gridded/"
        f"{filename}",
        wget_options=[],
    )
