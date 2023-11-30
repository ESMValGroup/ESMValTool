"""Script to download NOAA-ERSST-v3b."""
import logging

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

    downloader.download_file(
        "https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v3b/netcdf/ersst.v3b.1854-2019.nc",
        wget_options=[],
    )
    downloader.download_file(
        "https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v3b/netcdf/ersst.v3b.fldmean.1854-2019.nc",
        wget_options=[],
    )
