"""Script to download ESACCI-OC from the Climate Data Store(CDS)."""

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
    if not start_date:
        start_date = datetime(1997, 9, 1)
    if not end_date:
        end_date = datetime(2022, 12, 1)

    loop_date = start_date

    # See https://climate.esa.int/en/projects/ocean-colour/data/
    downloader = FTPDownloader(
        original_data_dir=original_data_dir,
        server="oceancolour.org",
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
        user="oc-cci-data",
        passwd="ELaiWai8ae",  # noqa: S106
    )
    downloader.connect()

    downloader.set_cwd("occci-v5.0/geographic/netcdf/monthly/chlor_a")
    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_folder(str(year))
        loop_date += relativedelta.relativedelta(years=1)
