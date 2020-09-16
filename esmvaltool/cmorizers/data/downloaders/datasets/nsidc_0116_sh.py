"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)."""

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import NASADownloader


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

    downloader = NASADownloader(
        config=config,
        dataset=dataset,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_file(
            "https://daacdata.apps.nsidc.org/pub/DATASETS/"
            "nsidc0116_icemotion_vectors_v4/south/daily/"
            f"icemotion_daily_sh_25km_{year}0101_{year}1231_v4.1.nc"
        )

        loop_date += relativedelta.relativedelta(years=1)
