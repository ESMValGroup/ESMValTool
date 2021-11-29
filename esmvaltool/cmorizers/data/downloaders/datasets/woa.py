"""Script to download WOA from its webpage."""
import logging

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


# pylint: disable=invalid-name
def download_dataset(config, dataset, dataset_info, _, __, overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
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

    def download(file):
        downloader.download_file(
            "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/" + file,
            wget_options=[])

    download("temperature/netcdf/decav81B0/1.00/woa18_decav81B0_t00_01.nc")
    download("salinity/netcdf/decav81B0/1.00/woa18_decav81B0_s00_01.nc")
    download("oxygen/netcdf/all/1.00/woa18_all_o00_01.nc")
    download("nitrate/netcdf/all/1.00/woa18_all_n00_01.nc")
    download("phosphate/netcdf/all/1.00/woa18_all_p00_01.nc")
    download("silicate/netcdf/all/1.00/woa18_all_i00_01.nc")
