"""Script to download Duveiller2018 from its webpage."""
import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, _, __, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.tier = 2
    os.makedirs(downloader.local_folder, exist_ok=True)
    downloader.download_file(
        "https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz",
        wget_options=[]

    )
    downloader.download_file(
        "https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_ice.nc.gz",
        wget_options=[]

    )

    unpack_files_in_folder(downloader.local_folder)
