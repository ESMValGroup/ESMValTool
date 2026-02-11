"""Script to download HadEX3 from its webpage."""

import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder

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

    web_source = "https://www.metoffice.gov.uk/hadobs/hadex3/data/"

    file_list = [
        "HadEX3_TXx_MON.nc.gz",
        "HadEX3_TXx_ANN.nc.gz",
        "HadEX3_TNn_MON.nc.gz",
        "HadEX3_TNn_ANN.nc.gz",
        "HadEX3_Rx1day_MON.nc.gz",
        "HadEX3_Rx1day_ANN.nc.gz",
        "HadEX3_Rx5day_MON.nc.gz",
        "HadEX3_Rx5day_ANN.nc.gz",
    ]

    for file in file_list:
        downloader.download_file(f"{web_source}{file}", wget_options=[])

    unpack_files_in_folder(downloader.local_folder)
