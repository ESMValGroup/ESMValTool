"""Script to download NSIDC-G02202-sh"""
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
        end_date = datetime(2023, 1, 1)

    loop_date = start_date

    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    # base_path = ("https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v3b/netcdf"
    #              "/ersst.{year}{month:02d}.nc")

    # !wget --ftp-user=anonymous -nc ftp://sidads.colorado.edu/DATASETS/seaice/polar-stereo/tools/pss25area_v3.dat
    anc_path = 'https://noaadata.apps.nsidc.org/NOAA/G02202_V4/ancillary/G02202-cdr-ancillary-sh.nc'
    downloader.download_folder(anc_path, [])

    # https://polarwatch.noaa.gov/erddap/files/nsidcG02202v4shmday/seaice_conc_monthly_sh_197901_n07_v04r00.nc
    base_path = ('https://noaadata.apps.nsidc.org/NOAA/G02202_V4/south/monthly'
                 '/seaice_conc_monthly_sh_{year}{month:02d}_{other}_v04r00.nc') # regex for n07 changes to f08..? need rules

    # {'197811':'n07','198708':'f08', '199201':'f11','199510':'f13', '200801':'f17'} #bins
    dateLs = [datetime(1978,11,1),datetime(1987,7,30),datetime(1991,12,30),datetime(1995,9,30),datetime(2007,12,30),end_date]
    suffLs = ['n07','f08','f11','f13','f17']
    isuf = 0
    suffix = suffLs[isuf]
    # initialize suffix if dates start higher than initial
    while loop_date >= dateLs[isuf]:
        suffix = suffLs[isuf]
        isuf += 1

    while loop_date <= end_date:
        
        if loop_date > dateLs[isuf]:
            logging.info('... index {}..'.format(isuf))
            suffix = suffLs[isuf]
            isuf += 1

        downloader.download_folder(
            base_path.format(year=loop_date.year, month=loop_date.month, other=suffix), [])
        loop_date += relativedelta.relativedelta(months=1)
        # check loop_date is => next bin
