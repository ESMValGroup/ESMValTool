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
    downloader.download_file(
        "https://data.iac.ethz.ch/landflux/"
        "LandFluxEVAL.merged.89-05.monthly.all.nc",
        wget_options=[]
    )
