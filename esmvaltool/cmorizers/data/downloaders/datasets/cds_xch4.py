"""Script to download CDS-XCH4 from the Climate Data Store (CDS)."""

from esmvaltool.cmorizers.data.downloaders.cds import CDSDownloader
from esmvaltool.cmorizers.data.utilities import unpack_files_in_folder


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
    downloader = CDSDownloader(
        product_name='satellite-methane',
        request_dictionary={
            'format': 'tgz',
            'processing_level': 'level_3',
            'variable': 'xch4',
            'sensor_and_algorithm': 'merged_obs4mips',
            'version': '4.1',
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.download_request("CDS-XCH4.tar")
    unpack_files_in_folder(downloader.local_folder)
