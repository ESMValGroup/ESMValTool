"""Script to download ESACCI-LANDCOVER pft data from the CEDA."""

from datetime import datetime

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader


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
    # Default start and end dates if not provided
    if not start_date:
        start_date = datetime(1992, 1, 1)
    if not end_date:
        end_date = datetime(2020, 12, 31)

    # Initialize the downloader
    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = 'land_cover'
    downloader.connect()

    # Set current working directory to the main directory with the files
    downloader.set_cwd('/pft/v2.0.8/')

    # Create a regex pattern to match any .nc files
    year_range = '|'.join(str(year) for year in range(start_date.year,
                                                      end_date.year + 1))
    pattern = rf".*-(?:{year_range}).*\.nc$"

    # Download all .nc files in the directory
    downloader.download_folder('.', filter_files=pattern)
