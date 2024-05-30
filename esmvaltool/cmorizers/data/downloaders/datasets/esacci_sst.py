"""Script to download ESACCI-SOILMOISTURE."""
# Import required python modules
import ftplib
import os

from datetime import datetime

from dateutil import relativedelta

from esmvaltool.cmorizers.data.downloaders.ftp import CCIDownloader
from esmvaltool.cmorizers.data.downloaders.ftp import FTPDownloader


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
        start_date = datetime(2004, 1, 1)
    if end_date is None:
        end_date = datetime(2007,12, 31)

    loop_date = start_date

    ## Define the local directory name to put data in
    #Ddir="C:\\datadir"

    ## If directory doesn't exist make it
    #If not os.path.isdir(ddir):
    #    os.mkdir(ddir)

    ## Change the local directory to where you want to put the data
    #Os.chdir(ddir)

    #user = os.environ.get("ceda-user")
    #if user is None:
    #    user = str(input("CEDA user name? "))
    #    if user == "":
    #        errmsg = ("A CEDA account is required to download CCI SST data."
    #                  " Please visit https://rda.ucar.edu/login/register/"
    #                  " to create an account at the Research Data Archive"
    #                  " (RDA) if needed.")
    #        logger.error(errmsg)
    #        raise ValueError

    #passwd = os.environ.get("ceda-passwd")
    #if passwd is None:
    #    passwd = str(input("CEDA-password? "))

    # login to FTP
    #f=ftplib.FTP("ftp3.ceda.ac.uk", user=user, passwd=passwd)
    #f=ftplib.FTP("ftp3.ceda.ac.uk", user="lbock", passwd="~?s`Cb`@E8Vf")
    #f=ftplib.FTP("ftp3.ceda.ac.uk", user="lbock", passwd="Pasw4ceda!")

    downloader = FTPDownloader(
        config=config,
        server='ftp3.ceda.ac.uk',
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )

    ##os.makedirs(downloader.local_folder, exist_ok=True)


    ##options = ["-O", "Authentication.log", "--save-cookies=auth.rda_ucar_edu",
    ##           f"--post-data=\"email={user}&passwd={passwd}&action=login\""]

    ### login to Research Data Archive (RDA)

    ##downloader.login("https://www.ceda.com/login", options)

    #downloader.ftp_name = 'eocis'
    ##downloader.ftp_name = 'sst'
    downloader.connect()
    ##downloader.set_cwd('CDR_v2/Analysis/L4/v2.1/')
    downloader.set_cwd('neodc/eocis/data/global_and_regional/sea_surface_temperature/CDR_v3/Analysis/L4/v3.0.1/')
    #downloader.set_cwd('')
    
    #f.cwd = 'neodc/eocis/data/global_and_regional/sea_surface_temperature/CDR_v3/Analysis/L4/v3.0.1/'

    while loop_date <= end_date:
        year = loop_date.year
        month = loop_date.strftime("%m")
        day = loop_date.strftime("%d")
        downloader.download_folder(f'./{year}/{month}/{day}/')
        loop_date += relativedelta.relativedelta(days=1)
