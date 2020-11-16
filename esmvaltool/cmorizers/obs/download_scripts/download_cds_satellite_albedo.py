"""Script to download cds-satellite-albedo from the Climate Data Store(CDS)"""

import argparse
import os
from datetime import datetime

import cdsapi
import yaml
from dateutil import relativedelta


def download_cds_satellite_albedo():
    """Download dataset cds-satellite-albedo."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config_file',
                        '-c',
                        default=os.path.join(os.path.dirname(__file__),
                                             'config-user.yml'),
                        help='Config file')
    parser.add_argument('--startdate',
                        '-s',
                        default='199901',
                        help='starting date as YYYYMM')
    parser.add_argument('--enddate',
                        '-e',
                        default='201405',
                        help='end date as YYYYMM')

    args = parser.parse_args()

    # get and read config file
    config_file_name = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    with open(config_file_name, 'r') as config_file:
        config = yaml.safe_load(config_file)

    rawobs_dir = os.path.abspath(
        os.path.expandvars(os.path.expanduser(config['rootpath']['RAWOBS'])))
    cds_satellite_albedo_dir = f'{rawobs_dir}/Tier3/CDS-SATELLITE-ALBEDO/'
    os.makedirs(cds_satellite_albedo_dir, exist_ok=True)

    client = cdsapi.Client()

    loopdate = datetime.strptime(args.startdate, '%Y%m')
    enddate = datetime.strptime(args.enddate, '%Y%m')

    while loopdate <= enddate:

        month = loopdate.month
        year = loopdate.year
        savename = os.path.join(
            cds_satellite_albedo_dir,
            f"cds-satellite-albedo_{year}{month:02d}.tar.gz")
        request_dictionary = {
            'format': 'tgz',
            'satellite': 'spot',
            'sensor': 'vgt',
            'product_version': 'V1',
            'horizontal_resolution': '1km',
            'variable': [
                'albb_bh',
                'albb_dh',
            ],
            'year': f"{year}",
            'month': f"{month:02d}",
            'nominal_day': '20',
        }
        client.retrieve('satellite-albedo', request_dictionary, savename)
        loopdate += relativedelta.relativedelta(months=1)


if __name__ == "__main__":
    download_cds_satellite_albedo()
