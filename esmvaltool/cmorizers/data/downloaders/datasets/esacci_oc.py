"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)."""

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader


def download_dataset(config, dataset, start_date, end_date, overwrite):
    """
    Download dataset.

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

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )
    downloader.ftp_name = 'ocean_colour'
    downloader.connect()

    downloader.set_cwd('v3.1-release/geographic/netcdf/chlor_a/monthly/v3.1/')
    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_year(f'{year}')
        loop_date += relativedelta.relativedelta(years=1)
