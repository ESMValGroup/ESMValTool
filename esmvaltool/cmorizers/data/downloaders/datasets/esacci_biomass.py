"""Script to download ESACCI-BIOMASS agb data from the CEDA."""

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader


def download_dataset(
    config,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
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
    # Initialize the downloader
    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = "biomass"
    downloader.connect()

    # Set current working directory to the main directory with the files
    downloader.set_cwd("/agb/maps/v6.0/netcdf")

    # Download 10 km file
    downloader.download_file("ESACCI-BIOMASS-L4-AGB-MERGED-10000m-fv6.0.nc")
