"""Script to download CRU from its webpage."""
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
    for var in ['tmp', 'pre']:
        downloader.download_file(
            'https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/'
            f'cruts.1811131722.v4.02/{var}/'
            f'cru_ts4.02.1901.2017.{var}.dat.nc.gz',
            wget_options=[]
        )
    unpack_files_in_folder(downloader.local_folder)
