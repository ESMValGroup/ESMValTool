"""Script to download APHRO-MA from its webpage."""
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

    def download_file(path):
        downloader.download_file(
            "http://aphrodite.st.hirosaki-u.ac.jp/product/" + path, []
        )

    for grid in ('025deg', '050deg'):
        download_file(
            f"APHRO_V1808_TEMP/APHRO_MA/{grid}_nc/"
            f"APHRO_MA_TAVE_{grid}_V1808.nc.tgz")
        download_file(
            f"APHRO_V1808_TEMP/APHRO_MA/{grid}_nc/"
            f"APHRO_MA_TAVE_{grid}_V1808.nc.tgz")
        download_file(
            f"APHRO_V1101/APHRO_MA/{grid}_nc/"
            f"APHRO_MA_{grid}_V1101.1951-2007.nc.gz.tar")
        download_file(
            f"APHRO_V1101/APHRO_MA/{grid}_nc/"
            f"APHRO_MA_{grid}_V1101.1951-2007.nc.gz.tar")
        download_file(
            f"APHRO_V1101EX_R1/APHRO_MA/{grid}_nc/"
            f"APHRO_MA_{grid}_V1101_EXR1.nc.tgz")
        download_file(
            f"APHRO_V1101EX_R1/APHRO_MA/{grid}_nc/"
            f"APHRO_MA_{grid}_V1101_EXR1.nc.tgz")

    unpack_files_in_folder(downloader.local_folder)
