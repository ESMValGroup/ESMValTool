"""Download Yang2020 data."""

import logging

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    config,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
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
        "https://darchive.mblwhoilibrary.org/bitstreams/"
        "c407fd4e-0fb9-52ba-93b8-5188994d02dd/download",
        wget_options=[],
        output_filename="n2oFlux-Yang2020.nc",
    )
    downloader.download_file(
        "https://darchive.mblwhoilibrary.org/bitstreams/"
        "266ae58f-b915-536e-af17-279d276e4295/download",
        wget_options=[],
        output_filename="dn2o-mapped-Yang2020.nc",
    )
