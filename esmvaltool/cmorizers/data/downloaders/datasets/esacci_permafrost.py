"""Script to download ESACCI-PERMAFROST."""

import datetime
import logging

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
) -> None:
    """Download dataset.

    Parameters
    ----------
    original_data_dir : Path
        Directory where original data will be stored.
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
        start_date = datetime.datetime(1997, 1, 1, tzinfo=datetime.UTC)
    if end_date is None:
        end_date = datetime.datetime(2023, 12, 31, tzinfo=datetime.UTC)

    downloader = WGetDownloader(
        original_data_dir=original_data_dir,
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

    # download cci variables active layer thickness, ground temperature and
    # permafrost extent
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
