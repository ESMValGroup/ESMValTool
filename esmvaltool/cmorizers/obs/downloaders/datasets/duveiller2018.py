"""Script to download APHRO-MA from its webpage."""
import datetime
from dateutil import relativedelta
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder

logger = logging.getLogger(__name__)

import os
import shutil


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
