"""Script to download APHRO-MA from its webpage."""
import shutil
import os
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.obs.utilities import unpack_files_in_folder

logger = logging.getLogger(__name__)

import os
import re
import gzip
import shutil


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

    # download_file("APHRO_V1808_TEMP/APHRO_MA/025deg_nc/APHRO_MA_TAVE_025deg_V1808.nc.tgz")
    # download_file("APHRO_V1808_TEMP/APHRO_MA/050deg_nc/APHRO_MA_TAVE_050deg_V1808.nc.tgz")
    # download_file("APHRO_V1101/APHRO_MA/025deg_nc/APHRO_MA_025deg_V1101.1951-2007.nc.gz.tar")
    # download_file("APHRO_V1101/APHRO_MA/050deg_nc/APHRO_MA_050deg_V1101.1951-2007.nc.gz.tar")
    # download_file("APHRO_V1101EX_R1/APHRO_MA/025deg_nc/APHRO_MA_025deg_V1101_EXR1.nc.tgz")
    # download_file("APHRO_V1101EX_R1/APHRO_MA/050deg_nc/APHRO_MA_050deg_V1101_EXR1.nc.tgz")

    shutil.register_unpack_format('gz', ['.gz', ], gunzip)

    unpack_files_in_folder(downloader.local_folder)
