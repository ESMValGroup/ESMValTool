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
        "https://data.icos-cp.eu/objects/OT_YY6iORypk2yAcwjMpUQpo",
        wget_options=[
            f"--output-document={downloader.local_folder}/"
            "Global_Carbon_Budget_2018v1.0.xlsx"]
    )
