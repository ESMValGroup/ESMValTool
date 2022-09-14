"""Script to download MERRA2."""

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
        start_date = datetime(1980, 1, 1)
    if not end_date:
        end_date = datetime(2021, 1, 1)
    loop_date = start_date

    downloader = NASADownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    while loop_date <= end_date:
        year = loop_date.year
        downloader.download_folder(
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2TMNXLND.5.12.4/{year}/")
        downloader.download_folder(
            "https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2IMNPANA.5.12.4/{year}/")
        downloader.download_folder(
            "https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2IMNPASM.5.12.4/{year}/")
        downloader.download_folder(
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2TMNXRAD.5.12.4/{year}/")
        downloader.download_folder(
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2TMNXSLV.5.12.4/{year}/")
        downloader.download_folder(
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2SMNXSLV.5.12.4/{year}/")
        downloader.download_folder(
            "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2_MONTHLY/"
            f"M2TMNXFLX.5.12.4/{year}/")
        loop_date += relativedelta.relativedelta(years=1)
