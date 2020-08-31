"""Script to download CowtanWay from its webpage."""
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
            "https://www-users.york.ac.uk/~kdc3/papers/coverage2013/" + path,
            []
        )

    download_file("had4_krig_v1_0_0.nc.gz")
    download_file("had4_uah_v1_0_0.nc.gz")
    download_file("had4_short_krig_v2_0_0.nc.gz")
    download_file("had4_short_uah_v2_0_0.nc.gz")
    download_file("ghcn_short_krig_v2_0_0.nc.gz")
    download_file("ghcn_short_uah_v2_0_0.nc.gz")
    download_file("had4sst4_krig_v2_0_0.nc.gz")
    download_file("had4_krig_v2_0_0.nc.gz")
    unpack_files_in_folder(downloader.local_folder)
