"""Script to download Duveiller2018 from its webpage."""
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.obs.utilities import (unpack_files_in_folder,
                                                read_cmor_config)

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, start_date, end_date, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.tier = 2
    cmor_config = read_cmor_config(dataset)
    raw_path = (
        "https://opendata.dwd.de/climate_environment/GPCC/"
        "full_data_2018/full_data_monthly_{version}.nc.gz"
    )
    for version in cmor_config['attributes']['version'].values():
        downloader.download_file(
            raw_path.format(version=version),
            wget_options=[]

        )
    unpack_files_in_folder(downloader.local_folder)
