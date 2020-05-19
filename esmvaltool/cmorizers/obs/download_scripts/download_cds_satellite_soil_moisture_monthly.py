"""
Download cds-satellite-soil-moisture.

Script to download cds-satellite-soil-moisture
from the Climate Data Store(CDS).
"""

import argparse
import os
import subprocess

import yaml

import cdsapi


def download_cds_satellite_soil_moisture():
    """Download dataset cds-satellite-albedo."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config_file',
                        '-c',
                        default=os.path.join(os.path.dirname(__file__),
                                             'config-user.yml'),
                        help='Config file')
    parser.add_argument('--frequency',
                        '-f',
                        default=None,
                        help="Either monthly or daily")

    args = parser.parse_args()

    # get and read config file
    config_file_name = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    with open(config_file_name, 'r') as config_file:
        config = yaml.safe_load(config_file)

    rawobs_dir = os.path.abspath(
        os.path.expandvars(os.path.expanduser(config['rootpath']['RAWOBS'])))
    savedir = f'{rawobs_dir}/Tier3/CDS-SATELLITE-SOIL-MOISTURE_tmp/'
    os.makedirs(savedir, exist_ok=True)

    client = cdsapi.Client()

    # Download per year
    # CDR from 1978-2018 (including 2018)
    for year in range(1978, 2019):
        savename = os.path.join(
            savedir,
            # Although tgz is requested, the file
            # actually is a tar (issue has been opened on CDS)
            f"cds-satellite-soil-moisture_{year}.tar")
        request_dictionary = {
            'format':
            'tgz',
            'variable':
            'volumetric_surface_soil_moisture',
            'type_of_sensor':
            'combined_passive_and_active',
            'time_aggregation':
            'month_average',
            'year': [
                f'{year}',
            ],
            'month': [
                '01',
                '02',
                '03',
                '04',
                '05',
                '06',
                '07',
                '08',
                '09',
                '10',
                '11',
                '12',
            ],
            'day':
            '01',
            'type_of_record':
            'cdr',
            'version':
            'v201812.0.0',
        }
        client.retrieve('satellite-soil-moisture', request_dictionary,
                        savename)
    # Unpack the file
    subprocess.check_call(["tar", "-xvf", savename, '--directory', savedir])
    # Remove the tar file since it has been extracted
    subprocess.check_call(["rm", savename])

    # Now ICDR v201812.0.1 (identical to v201812.0.0 according
    # to data provider except for dates covered)
    savename = os.path.join(savedir,
                            f"cds-satellite-soil-moisture_v201812.0.1.tar.gz")
    client.retrieve(
        'satellite-soil-moisture', {
            'format': 'tgz',
            'variable': 'volumetric_surface_soil_moisture',
            'type_of_sensor': 'combined_passive_and_active',
            'time_aggregation': 'month_average',
            'year': '2019',
            'type_of_record': 'icdr',
            'month': [
                '01',
                '02',
                '03',
                '04',
                '05',
                '06',
                '07',
                '08',
                '09',
            ],
            'version': 'v201812.0.1',
            'day': '01',
        }, savename)
    # Unpack the file
    subprocess.check_call(["tar", "-xvf", savename, '--directory', savedir])
    # Remove the tar file since it has been extracted
    subprocess.check_call(["rm", savename])

    # Now ICDR v201812.0.0 (identical to v201812.0.1 according
    # to data provider except for dates covered)
    savename = os.path.join(savedir,
                            f"cds-satellite-soil-moisture_v201812.0.0.tar")
    client.retrieve(
        'satellite-soil-moisture', {
            'format': 'tgz',
            'variable': 'volumetric_surface_soil_moisture',
            'type_of_sensor': 'combined_passive_and_active',
            'time_aggregation': 'month_average',
            'year': [
                '2019',
                '2020',
            ],
            'type_of_record': 'icdr',
            'month': [
                '01',
                '02',
                '03',
                '04',
                '10',
                '11',
                '12',
            ],
            'version': 'v201812.0.0',
            'day': '01',
        }, savename)
    # Unpack the file
    subprocess.check_call(["tar", "-xvf", savename, '--directory', savedir])
    # Remove the tar file since it has been extracted
    subprocess.check_call(["rm", savename])


if __name__ == "__main__":
    download_cds_satellite_soil_moisture()
