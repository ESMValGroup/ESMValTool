"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta

from esmvaltool.cmorizers.obs.downloaders.wget import NASADownloader


def download_dataset(config, dataset, start_date, end_date):
    """Download dataset cds-satellite-albedo."""
    loop_date = start_date

    downloader = NASADownloader(
        config=config,
        dataset=dataset
    )

    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_folder(
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2TMNXLND.5.12.4/{year}/"
        )
        loop_date += relativedelta.relativedelta(years=1)
