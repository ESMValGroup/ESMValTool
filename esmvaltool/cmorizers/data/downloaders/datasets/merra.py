"""Script to download MERRA."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import NASADownloader


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
    if not start_date:
        start_date = datetime(1979, 1, 1)
    if not end_date:
        end_date = datetime(2015, 12, 31)
    loop_date = start_date

    downloader = NASADownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_folder(
            "https://goldsmr3.gesdisc.eosdis.nasa.gov/data/MERRA_MONTHLY/"
            f"MAIMNXINT.5.2.0/{year}/",
        )
        downloader.download_folder(
            "https://goldsmr3.gesdisc.eosdis.nasa.gov/data/MERRA_MONTHLY/"
            f"MAIMCPASM.5.2.0/{year}/",
        )
        downloader.download_folder(
            "https://goldsmr3.gesdisc.eosdis.nasa.gov/data/MERRA_MONTHLY/"
            f"MATMNXRAD.5.2.0/{year}/",
        )
        downloader.download_folder(
            "https://goldsmr3.gesdisc.eosdis.nasa.gov/data/MERRA_MONTHLY/"
            f"MATMFXCHM.5.2.0/{year}/",
        )

        loop_date += relativedelta.relativedelta(years=1)
