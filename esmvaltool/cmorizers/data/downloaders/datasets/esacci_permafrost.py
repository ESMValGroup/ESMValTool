"""Script to download ESACCI-PERMAFROST."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    config: dict,
    dataset: str,
    dataset_info: dict,
    start_date: datetime,
    end_date: datetime,
    overwrite: bool,
) -> None:
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
    if start_date is None:
        start_date = datetime(1997, 1, 1)
    if end_date is None:
        end_date = datetime(2023, 12, 31)

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    version = "v05.0"

    path = "https://dap.ceda.ac.uk/neodc/esacci/permafrost/data/"

    ccivars = [
        "active_layer_thickness",
        "ground_temperature",
        "permafrost_extent",
    ]

    # download active layer thickness
    loop_date = start_date
    while loop_date <= end_date:
        for var in ccivars:
            folder = (
                path
                + f"{var}/L4/area4/pp/{version}/northern_hemisphere/{loop_date.year}/"
            )
            downloader.download_folder(
                folder,
                wget_options=["-e robots=off", "--no-parent", "--accept=nc"],
            )
        loop_date += relativedelta.relativedelta(years=1)
