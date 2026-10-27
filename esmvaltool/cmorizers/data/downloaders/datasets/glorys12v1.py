"""Downloader for the GLORYS12V1 dataset from the Copernicus Marine Service."""

import datetime
import logging

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.cms import CMSDownloader

logger = logging.getLogger(__name__)


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
    # Download variables
    downloader = CMSDownloader(
        original_data_dir=original_data_dir,
        product_name="cmems_mod_glo_phy_my_0.083deg_P1M-m",
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
        no_directories=True,
        product_filename="mercatorglorys12v1_gl12_mean",
    )
    if start_date is None:
        start_date = datetime.datetime(1993, 1, 1, tzinfo=datetime.UTC)
    if end_date is None:
        end_date = datetime.datetime(2026, 4, 30, tzinfo=datetime.UTC)

    # Download the whole dataset
    if (start_date == datetime.datetime(1993, 1, 1, tzinfo=datetime.UTC)) and (
        end_date == datetime.datetime(2026, 4, 30, tzinfo=datetime.UTC)
    ):
        downloader.download()
    # Loop over years to download
    else:
        loop_date = start_date
        while loop_date <= end_date:
            downloader.download_year(year=loop_date.year)
            loop_date += relativedelta.relativedelta(years=1)
