"""Script to download Duveiller2018 from its webpage."""
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, start_date, end_date, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    for var in ['tmp', 'pre']:
        downloader.download_file(
            'https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/'
            '9969496/albedo_IGBPgen.nc',
            wget_options=[]
        )
