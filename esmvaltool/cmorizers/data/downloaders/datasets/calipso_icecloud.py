"""Script to download CALIPSO-ICECLOUD from its webpage."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import NASADownloader


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
    if not start_date:
        start_date = datetime(2007, 1, 1)
    if not end_date:
        end_date = datetime(2015, 12, 31)
    loop_date = start_date

    downloader = NASADownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        year = loop_date.year
        for month in range(1, 13):

            downloader.download_file(
                "https://asdc.larc.nasa.gov/data/CALIPSO/"
                f"LID_L3_Ice_Cloud-Standard-V1-00/{year}/"
                f"CAL_LID_L3_Ice_Cloud-Standard-V1-00.{year}-{month:02}A.hdf")
        loop_date += relativedelta.relativedelta(years=1)
