"""Script to download PERSIANN-CDR."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    start_date : datetime
        Start of the interval to download
    end_date : datetime
        End of the interval to download
    overwrite : bool
        Overwrite already downloaded files
    """
    if start_date is None:
        start_date = datetime(1983, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 1, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    base_path = (
        "https://www.ncei.noaa.gov/data/precipitation-persiann/access/"
        "{year}/")
    while loop_date <= end_date:
        downloader.download_folder(base_path.format(year=loop_date.year), [])
        loop_date += relativedelta.relativedelta(years=1)
