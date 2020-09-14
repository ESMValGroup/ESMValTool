"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

from datetime import datetime
from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset ESACCI-AEROSOL."""
    if start_date is None:
        start_date = datetime(1982, 1, 1)
    if end_date is None:
        end_date = datetime(2016, 1, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    base_path = (
        "https://www.ncei.noaa.gov/data/"
        "avhrr-reflectance-cloud-properties-patmos-extended/access/{year}/")
    while loop_date <= end_date:
        downloader.download_folder(
            base_path.format(year=loop_date.year),
            # ["--accept='*NOAA*.nc'", "--reject='*preliminary*'"]
            []
        )
        loop_date += relativedelta.relativedelta(years=1)
