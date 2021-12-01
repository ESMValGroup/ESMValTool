"""Script to download cds-satellite-albedo from the CDS."""

from dateutil import relativedelta

from ..cds import CDSDownloader


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
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
    loop_date = start_date

    downloader = CDSDownloader(
        product_name='satellite-sea-surface-temperature-ensemble-product',
        request_dictionary={
            'variable': 'all',
            'format': 'tar',
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        downloader.download(loop_date.year, loop_date.month, loop_date.day)
        loop_date += relativedelta.relativedelta(days=1)
