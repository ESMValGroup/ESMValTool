"""Script to download Yang2020 from its webpage."""
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
        "https://datadocs.bco-dmo.org/"
        "data/302/Multiple_Scales_N_Cycle/810032/1/data/"
        "dn2o-mapped-Yang2020.nc",
        wget_options=[],
    )
    downloader.download_file(
        "https://datadocs.bco-dmo.org/"
        "data/302/Multiple_Scales_N_Cycle/810032/1/data/"
        "n2oFlux-Yang2020.nc",
        wget_options=[],
    )
