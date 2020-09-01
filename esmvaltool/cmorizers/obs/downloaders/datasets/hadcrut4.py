"""Script to download Duveiller2018 from its webpage."""
import logging
import os

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, start_date, end_date, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.tier = 2
    os.makedirs(downloader.local_folder, exist_ok=True)
    downloader.download_file(
        "https://crudata.uea.ac.uk/cru/data/temperature/"
        "HadCRUT.4.6.0.0.median.nc",
        wget_options=[]

    )
    downloader.download_file(
        "https://crudata.uea.ac.uk/cru/data/temperature/absolute.nc",
        wget_options=[]

    )
