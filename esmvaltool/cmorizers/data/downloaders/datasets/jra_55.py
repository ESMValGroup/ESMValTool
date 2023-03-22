"""Script to download JRA-55 from RDA."""
import logging
import os

from esmvaltool.cmorizers.data.downloaders.wget import WGetDownloader

from datetime import datetime

from dateutil import relativedelta


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
    downloader = WGetDownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    os.makedirs(downloader.local_folder, exist_ok=True)

#    user = os.environ.get("rda-user")
#    if (user is None):
#        user = str(input("RDA user name? "))
#        if (user is ""):
#            print("A RDA account is required to download JRA-55 data.")
#            print("Please visit https://rda.ucar.edu/login/register/")
#            print("to create an account at the Research Data Archive (RDA)")
#            print("if needed.")
#            exit()
#
#    passwd = os.environ.get("rda-passwd")
#    if (passwd is None):
#        passwd = str(input("RDA password? "))

    if start_date is None:
#        start_date = datetime(1958, 1, 1)
        start_date = datetime(2022, 12, 1)
    if end_date is None:
        end_date = datetime(2022, 12, 31)
    loop_date = start_date

    options = ["-O", "Authentication.log", "--save-cookies=auth.rda_ucar_edu",
        f"--post-data=\"email={user}&passwd={passwd}&action=login\""]

    # login to Research Data Archive (RDA)

#    downloader.login("https://rda.ucar.edu/cgi-bin/login", options)

    # download files

    url = "https://rda.ucar.edu/data/ds628.1"
    download_options = ["--load-cookies=auth.rda_ucar_edu"]
    path = downloader.local_folder

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.month

        fname = f"anl_p125.039_vvel.{year}01_{year}12"
        print(fname)
        downloader.download_file(url + f"/anl_p125/{year}/" +
            fname, download_options)
        os.rename(downloader.local_folder + "/" + fname,
            downloader.local_folder + "/" + fname + ".grb") 

        

        loop_date += relativedelta.relativedelta(months=1)

    # add extension ".grb" to downloaded files

#    files = os.listdir(path)
#
#    for index, file in enumerate(files):
#        os.rename(os.path.join(path, file),
#            os.path.join(path, ''.join([file, '.grb'])))
