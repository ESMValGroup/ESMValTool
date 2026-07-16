"""Script to download OSI-450-sh."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader


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
        start_date = datetime(1979, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 12, 1)
    downloader = FTPDownloader(
        original_data_dir=original_data_dir,
        server="osisaf.met.no",
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.connect()

    loop_date = start_date

    downloader.set_cwd("/reprocessed/ice/conc/v3p0")

    while loop_date <= end_date:
        year = loop_date.year
        folder = f"{year}/{loop_date.month:02}"
        downloader.download_folder(
            folder,
            sub_folder=folder,
            filter_files=".*_sh.*[.]nc",
        )
        loop_date += relativedelta.relativedelta(months=1)
