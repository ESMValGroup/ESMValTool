"""Script to download CAMS greenhouse gas data from the Climate Data Store."""

import datetime as dt

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
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
        start_date = dt.datetime(year=1979, month=1, day=1)
    if end_date is None:
        end_date = dt.datetime(year=2023, month=12, day=31)

    downloader = CDSDownloader(
        product_name="cams-global-greenhouse-gas-inversion",
        request_dictionary={
            "variable": "carbon_dioxide",
            "quantity": "surface_flux",
            "input_observations": "surface",
            "time_aggregation": "monthly_mean",
            "version": "v23r1",
        },
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
        cds_url="https://ads.atmosphere.copernicus.eu/api",
    )

    for year in range(start_date.year, end_date.year + 1, 1):
        downloader.download_year(year)

    unpack_files_in_folder(downloader.local_folder)
