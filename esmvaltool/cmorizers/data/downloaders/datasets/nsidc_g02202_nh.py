# pylint: disable=too-many-arguments
# pylint: disable=too-many-function-args
# pylint: disable=too-many-positional-arguments
# pylint: disable=too-many-locals
"""Script to download NSIDC-G02202-nh."""
import logging
from datetime import datetime
from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


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
        start_date = datetime(1979, 1, 1)
    if end_date is None:
        end_date = datetime(2024, 6, 1)

    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    # need area file
    area_dat = ('ftp://sidads.colorado.edu/DATASETS/seaice'
                '/polar-stereo/tools/psn25area_v3.dat')
    downloader.download_folder(area_dat, [])

    anc_path = ('https://noaadata.apps.nsidc.org/NOAA/G02202_V5/'
                'ancillary/G02202-ancillary-psn25-v05r00.nc')
    downloader.download_folder(anc_path, [])

    base_path = ('https://noaadata.apps.nsidc.org/NOAA/G02202_V5/north/monthly'
                 '/sic_psn25_{year}{month:02d}_{other}_v05r00.nc')

    datels = [datetime(1978, 11, 1), datetime(1987, 7, 30),
              datetime(1991, 12, 30), datetime(1995, 9, 30),
              datetime(2007, 12, 30), end_date]
    suffls = ['n07', 'F08', 'F11', 'F13', 'F17']
    isuf = 0
    suffix = suffls[isuf]
    # initialize suffix if dates start higher than initial
    while loop_date >= datels[isuf]:
        suffix = suffls[isuf]
        isuf += 1

    while loop_date <= end_date:

        if loop_date > datels[isuf]:
            suffix = suffls[isuf]
            isuf += 1

        downloader.download_folder(
            base_path.format(year=loop_date.year, month=loop_date.month,
                             other=suffix), [])
        loop_date += relativedelta.relativedelta(months=1)
