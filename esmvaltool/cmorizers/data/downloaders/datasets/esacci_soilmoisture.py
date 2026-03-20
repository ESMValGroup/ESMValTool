"""Script to download ESACCI-SOILMOISTURE."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader


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
        start_date = datetime(1978, 11, 1)
    if end_date is None:
        end_date = datetime(2022, 12, 31)

    loop_date = start_date

    downloader = CCIDownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = "soil_moisture"
    downloader.connect()
    downloader.set_cwd("ancillary/v08.1/")
    downloader.download_folder(".")
    downloader.set_cwd("daily_files/COMBINED/v08.1/")
    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_year(f"{year}")
        loop_date += relativedelta.relativedelta(years=1)
