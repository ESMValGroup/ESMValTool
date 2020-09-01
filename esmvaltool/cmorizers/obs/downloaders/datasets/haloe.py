"""Script to download Duveiller2018 from its webpage."""
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, start_date, end_date, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.tier = 2
    downloader.download_file(
        "https://acp.copernicus.org/articles/5/2797/2005/"
        "acp-5-2797-2005-supplement.tar",
        wget_options=[]
    )
    unpack_files_in_folder(downloader.local_folder)
