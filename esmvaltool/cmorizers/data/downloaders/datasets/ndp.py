"""Script to download NDP from its webpage."""
import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder

logger = logging.getLogger(__name__)


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
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    if not os.path.isdir(downloader.local_folder):
        os.makedirs(downloader.local_folder)

    downloader.download_file(
        "https://data.ess-dive.lbl.gov/catalog/d1/mn/v2/object/"
        "ess-dive-ec4f4b7097524f6-20180621T213642471",
        ["-O", os.path.join(downloader.local_folder, "ndp017b.tar.gz")])

    unpack_files_in_folder(downloader.local_folder)
