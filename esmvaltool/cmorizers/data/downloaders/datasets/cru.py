"""Script to download CRU from its webpage."""

import logging

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
    for var in ["tmp", "pre", "pet", "tmn", "tmx", "cld"]:
        downloader.download_file(
            "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.07/"
            f"cruts.2304141047.v4.07/{var}/"
            f"cru_ts4.07.1901.2022.{var}.dat.nc.gz",
            wget_options=[],
        )
    # for var in ['tmp', 'pre']:  # v TS4.02
    #     downloader.download_file(
    #         "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/"
    #         f"cruts.1811131722.v4.02/{var}/"
    #         f"cru_ts4.02.1901.2017.{var}.dat.nc.gz",
    #         wget_options=[],
    #     )
    unpack_files_in_folder(downloader.local_folder)
