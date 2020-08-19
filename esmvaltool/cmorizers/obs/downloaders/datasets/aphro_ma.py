"""Script to download APHRO-MA from its webpage."""
import shutil
import os
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)

import os
import re
import gzip
import shutil


def gunzip(file_name, work_dir):
    filename = os.path.split(file_name)[-1]
    filename = re.sub(r"\.gz$", "", filename, flags=re.IGNORECASE)

    with gzip.open(file_name, 'rb') as f_in:
        with open(os.path.join(work_dir, filename), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

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
    decompress =True
    while decompress:
        decompress = False
        for filename in os.listdir(downloader.local_folder):
            full_path = os.path.join(downloader.local_folder, filename)
            if os.path.isdir(full_path):
                logger.info('Moving files from folder %s', filename)
                folder_files = os.listdir(full_path)
                for file_path in folder_files:
                    shutil.move(
                        os.path.join(full_path, file_path),
                        downloader.local_folder
                    )
                os.rmdir(full_path)
                decompress = True
                continue
            if not filename.endswith(('.gz', '.tgz', '.tar')):
                continue
            logger.info('Unpacking %s', filename)
            shutil.unpack_archive(full_path, downloader.local_folder)
            os.remove(full_path)
            decompress = True
