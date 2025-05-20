"""Script to download CAMS greenhouse gas data from the Climate Data Store."""

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


def download_dataset(
    config, dataset, dataset_info, start_date, end_date, overwrite
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

    for year in range(1979, 2024, 1):
        downloader.download_year(year)

    unpack_files_in_folder(downloader.local_folder)
