"""Script to download WOA from its webpage."""
import logging
import os
import shutil

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

    def download(file):
        downloader.download_file(
            "https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA/" + file,
            wget_options=[])

    data_paths = [
        "nitrate/netcdf/all/1.00/woa18_all_n00_01.nc",
        "oxygen/netcdf/all/1.00/woa18_all_o00_01.nc",
        "phosphate/netcdf/all/1.00/woa18_all_p00_01.nc",
        "salinity/netcdf/decav81B0/1.00/woa18_decav81B0_s00_01.nc",
        "silicate/netcdf/all/1.00/woa18_all_i00_01.nc",
        "temperature/netcdf/decav81B0/1.00/woa18_decav81B0_t00_01.nc"
    ]

    for source_file in data_paths:
        download(source_file)
        filename = os.path.basename(source_file)
        var = source_file.split("/", maxsplit=1)[0]
        os.makedirs(os.path.join(downloader.local_folder, var), exist_ok=True)
        filepath = os.path.join(downloader.local_folder, filename)
        shutil.move(filepath,
                    os.path.join(downloader.local_folder, var, filename))
