"""Script to download Eppley-VGPM-MODIS."""
import datetime
import logging

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    dataset_info : dict
         Dataset information from the datasets.yml file
    start_date : datetime
        Start of the interval to download
    end_date : datetime
        End of the interval to download
    overwrite : bool
        Overwrite already downloaded files
    """
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
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
            wget_options=["--accept=tar"])
        loop_date += relativedelta.relativedelta(years=1)
    unpack_files_in_folder(downloader.local_folder)
