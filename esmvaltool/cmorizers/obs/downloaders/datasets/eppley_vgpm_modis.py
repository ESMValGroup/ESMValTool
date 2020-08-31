"""Script to download APHRO-MA from its webpage."""
import datetime
from dateutil import relativedelta
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

    if not start_date:
        start_date = datetime.datetime(2002, 1, 1)
    if not end_date:
        end_date = datetime.datetime(2020, 1, 1)

    loop_date = start_date
    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_folder(
            "http://orca.science.oregonstate.edu/data/1x2/monthly/"
            f"eppley.r2018.m.chl.m.sst/hdf/eppley.m.{year}.tar",
            wget_options=["--accept=tar"]
        )
        loop_date += relativedelta.relativedelta(years=1)
    unpack_files_in_folder(downloader.local_folder)
