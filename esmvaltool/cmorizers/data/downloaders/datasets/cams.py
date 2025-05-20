"""Script to download CAMS greenhouse gas data from the Climate Data Store."""

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


def download_dataset(
    config, dataset, dataset_info, overwrite, start_year=1979, end_year=2023
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
    overwrite : bool
        Overwrite already downloaded files
    start_year : datetime
        Start of the interval to download
    end_year : datetime
        End of the interval to download
    """

    downloader = CDSDownloader(
        product_name="cams-global-greenhouse-gas-inversion",
        request_dictionary={
            "variable": "carbon_dioxide",
            "quantity": "surface_flux",
            "input_observations": "surface",
            "time_aggregation": "monthly_mean",
            "version": "v23r1",
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
        cds_url="https://ads.atmosphere.copernicus.eu/api",
    )

    for year in range(start_year, end_year + 1, 1):
        downloader.download_year(year)

    unpack_files_in_folder(downloader.local_folder)
