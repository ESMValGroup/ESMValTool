"""Script to download HadCRUT3 from its webpage."""
import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


# pylint: disable=invalid-name
def download_dataset(config, dataset, dataset_info, _, __, overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
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

    os.makedirs(downloader.local_folder, exist_ok=True)
    downloader.download_file(
        "https://www.metoffice.gov.uk/hadobs/hadcrut3/data/HadCRUT3v.nc",
        wget_options=[])
