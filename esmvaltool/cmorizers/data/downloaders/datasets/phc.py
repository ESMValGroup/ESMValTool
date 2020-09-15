"""Script to download Duveiller2018 from its webpage."""
import logging

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
    downloader.download_file(
        "http://psc.apl.washington.edu/nonwp_projects/PHC/Data3/"
        "phc3.0_annual.nc",
        wget_options=[]
    )
    unpack_files_in_folder(downloader.local_folder)
