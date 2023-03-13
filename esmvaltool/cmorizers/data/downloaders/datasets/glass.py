"""Script to download GLASS."""
import logging
from datetime import datetime

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

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
    if start_date is None:
        start_date = datetime(1981, 1, 1)
    if end_date is None:
        end_date = datetime(2021, 12, 31)
    start_year = start_date.year
    end_year = end_date.year

    root = 'http://www.glass.umd.edu/'
    dirs = {
        "LAI/AVHRR/{year}/": (1981, 2018),
        "GPP/AVHRR/{year}/": (1982, 2018),
    }
    for year in range(start_year, end_year + 1):
        logger.info("Started download for year %i", year)
        for (dir, years) in dirs.items():
            # Skip variable if years are not supported
            if not (years[0] <= year <= years[1]):
                continue
            server_path = root + dir.format(year=year)
            downloader.download_folder(
                server_path, wget_options=['-np', '-A', '*.hdf']
            )
