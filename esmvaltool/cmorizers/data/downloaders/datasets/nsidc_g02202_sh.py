"""Script to download NSIDC-G02202-sh."""

import logging
from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

logger = logging.getLogger(__name__)


def download_dataset(
    original_data_dir,
    dataset,
    dataset_info,
    start_date,
    end_date,
    overwrite,
):
    """Download dataset.

    Parameters
    ----------
    original_data_dir : Path
        Directory where original data will be stored.
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
        end_date = datetime(2023, 1, 1)

    loop_date = start_date

    downloader = WGetDownloader(
        original_data_dir=original_data_dir,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    # need area file
    area_dat = (
        "ftp://sidads.colorado.edu/DATASETS/seaice"
        "/polar-stereo/tools/pss25area_v3.dat"
    )
    downloader.download_folder(area_dat, [])

    anc_path = (
        "https://noaadata.apps.nsidc.org/NOAA/G02202_V4/"
        "ancillary/G02202-cdr-ancillary-sh.nc"
    )
    downloader.download_folder(anc_path, [])

    base_path = (
        "https://noaadata.apps.nsidc.org/NOAA/G02202_V4/south/monthly"
        "/seaice_conc_monthly_sh_{year}{month:02d}_{other}_v04r00.nc"
    )

    # regex for n07 changes to f08.. file names
    # bins #{'197811':'n07','198708':'f08',
    # '199201':'f11','199510':'f13', '200801':'f17'}
    datels = [
        datetime(1978, 11, 1),
        datetime(1987, 7, 30),
        datetime(1991, 12, 30),
        datetime(1995, 9, 30),
        datetime(2007, 12, 30),
        end_date,
    ]
    suffls = ["n07", "f08", "f11", "f13", "f17"]
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
            base_path.format(
                year=loop_date.year,
                month=loop_date.month,
                other=suffix,
            ),
            [],
        )
        loop_date += relativedelta.relativedelta(months=1)
        # check loop_date is => next bin
