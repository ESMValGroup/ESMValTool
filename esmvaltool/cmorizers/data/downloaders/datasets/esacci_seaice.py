"""Script to download ESACCI-SEAICE."""

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader
from esmvaltool.cmorizers.data.utilities import read_cmor_config


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
    if start_date is None:
        start_date = datetime(1991, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 12, 31)

    downloader = CCIDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    downloader.ftp_name = "sea_ice"
    downloader.connect()

    atts = read_cmor_config(dataset)["variables"]
    for att in atts:
        if "regions" in atts[att]:
            regions = atts[att]["regions"]
        else:
            regions = ("NH", "SH")

    basepath = "sea_ice_concentration/L4/ssmi_ssmis/12.5km/v3.0"

    loop_date = start_date
    while loop_date <= end_date:
        for region in regions:
            path = (
                f"{basepath}/{region}/{loop_date.year}/{loop_date.month:02d}"
            )
            downloader.set_cwd(path)
            downloader.download_folder(".", sub_folder=region)
        loop_date += relativedelta.relativedelta(months=1)
