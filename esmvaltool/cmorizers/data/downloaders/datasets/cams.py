"""Script to download CAMS data from the Climate Data Store."""

import datetime

from dateutil import relativedelta

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
    if start_date is None:
        start_date = datetime.datetime(1979, 1, 1)
    if end_date is None:
        end_date = datetime.datetime(2020, 12, 1)
        #end_date = datetime.datetime(2022, 12, 1)

    downloader = CDSDownloader(
        product_name='cams-global-greenhouse-gas-inversion',
        request_dictionary={
            'variable': 'carbon_dioxide',
            'quantity': 'surface_flux',
            'input_observations': 'surface',
            'time_aggregation': 'monthly_mean',
            'version': 'v20r2' #v22r2
        },
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    loop_date = start_date
    while loop_date <= end_date:
        downloader.download(loop_date.year, loop_date.month, file_format='zip')
        loop_date += relativedelta.relativedelta(months=1)

    unpack_files_in_folder(downloader.local_folder)