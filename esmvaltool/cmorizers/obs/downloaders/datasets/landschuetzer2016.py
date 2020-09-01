"""Script to download Duveiller2018 from its webpage."""
import logging

from esmvaltool.cmorizers.obs.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(config, dataset, start_date, end_date, overwrite):
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.tier = 2
    downloader.download_file(
        "https://www.nodc.noaa.gov/archive/arc0105/0160558/3.3/data/0-data/"
        "spco2_1982-2015_MPI_SOM-FFN_v2016.nc",
        wget_options=[]
    )
