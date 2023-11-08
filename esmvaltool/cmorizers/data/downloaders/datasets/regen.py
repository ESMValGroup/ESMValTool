"""Script to download REGEN."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader
from esmvaltool.cmorizers.data.utilities import read_cmor_config


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
        start_date = datetime(1950, 1, 1)
    if end_date is None:
        end_date = datetime(2016, 1, 1)
    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    base_path = ("http://dapds00.nci.org.au/thredds/fileServer/ks32/CLEX_Data/"
                 "REGEN_AllStns/v1-2019/REGEN_AllStns_{version}_{year}.nc")
    version = read_cmor_config(dataset)['attributes']['version']
    while loop_date <= end_date:
        downloader.download_folder(
            base_path.format(year=loop_date.year, version=version), [])
        loop_date += relativedelta.relativedelta(years=1)
