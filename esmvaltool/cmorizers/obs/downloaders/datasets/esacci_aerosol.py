"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

from datetime import datetime
from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.ftp import CCIDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset ESACCI-AEROSOL."""
    if start_date is None:
        start_date = datetime(1997, 1, 1)
    if end_date is None:
        end_date = datetime(2011, 1, 1)
    loop_date = start_date

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.connect()



    while loop_date <= end_date:
        year = loop_date.year
        if year < 2002:
            downloader.set_cwd('ATSR2_SU/L3/v4.21/MONTHLY')
        else:
            downloader.set_cwd('AATSR_SU/L3/v4.21/MONTHLY')

        downloader.download_year(loop_date.year)
        loop_date += relativedelta.relativedelta(years=1)
