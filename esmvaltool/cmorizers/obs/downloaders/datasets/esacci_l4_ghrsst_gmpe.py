"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

from dateutil import relativedelta

from ..cds import CDSDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """Download dataset cds-satellite-albedo."""
    loop_date = start_date

    downloader = CDSDownloader(
        product_name='satellite-sea-surface-temperature-ensemble-product',
        request_dictionary={
            'variable': 'all',
            'format': 'tar',
        },
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        downloader.download(loop_date.year, loop_date.month, loop_date.day)
        loop_date += relativedelta.relativedelta(days=1)
