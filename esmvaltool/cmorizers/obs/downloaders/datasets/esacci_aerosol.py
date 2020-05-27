"""Script to download ESACCI-AEROSOL from CCI CEDA ftp."""

from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.ftp import CCIDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset ESACCI-AEROSOL."""
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
