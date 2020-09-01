"""Script to download Duveiller2018 from its webpage."""
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
    downloader.tier = 2
    def download(file):
        downloader.download_file(
            "https://data.nodc.noaa.gov/woa/WOA13/DATAv2/" + file,
            wget_options=[]

        )

    download("temperature/netcdf/decav81B0/1.00/woa13_decav81B0_t00_01.nc")
    download("salinity/netcdf/decav81B0/1.00/woa13_decav81B0_s00_01.nc")
    download("oxygen/netcdf/all/1.00/woa13_all_o00_01.nc")
    download("nitrate/netcdf/all/1.00/woa13_all_n00_01.nc")
    download("phosphate/netcdf/all/1.00/woa13_all_p00_01.nc")
    download("silicate/netcdf/all/1.00/woa13_all_i00_01.nc")
