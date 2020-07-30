"""
Download cds-satellite-soil-moisture.

Script to download cds-satellite-soil-moisture
from the Climate Data Store(CDS).
"""

import argparse
import os
import tarfile

import yaml

import cdsapi


def download_cds_satellite_sm():
    """Download dataset cds-satellite-soil-moisture."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config_file',
                        '-c',
                        help='config file',
                        required=True)
    parser.add_argument('--frequency',
                        '-f',
                        help='either monthly or daily',
                        required=True)

    args = parser.parse_args()

    # get and read config file
    config_file_name = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    with open(config_file_name, 'r') as config_file:
        config = yaml.safe_load(config_file)

    rawobs_dir = os.path.abspath(
        os.path.expandvars(os.path.expanduser(config['rootpath']['RAWOBS'])))
    savedir = f'{rawobs_dir}/Tier3/CDS-SATELLITE-SOIL-MOISTURE/'
    os.makedirs(savedir, exist_ok=True)

    # This dictionary contains the data request part that is specific to
    # the frequency (either daily/monthly) requested.
    freq_specific_kwargs = {
        'monthly': {
            'time_aggregation': 'month_average',
            'day': ['01']
        },
        'daily': {
            'time_aggregation': 'day_average',
            'day': [f'{i:02d}' for i in range(1, 32)],
        }
    }

    client = cdsapi.Client()

    # Download per year
    # CDR from 1978-2019 (including 2019)
    for year in range(1978, 2020):
        savename = os.path.join(
            savedir,
            # Although tgz is requested, the file
            # actually is a tar (issue has been opened on CDS)
            f"cds-satellite-soil-moisture_cdr_{args.frequency}_{year}.tar")
        request_dictionary = {
            'format':
            'tgz',
            'variable':
            'volumetric_surface_soil_moisture',
            'type_of_sensor':
            'combined_passive_and_active',
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
            'v201912.0.0',
            **freq_specific_kwargs[args.frequency]
        }
        client.retrieve('satellite-soil-moisture', request_dictionary,
                        savename)
        # Unpack the file
        tar = tarfile.open(savename)
        tar.extractall(path=os.path.dirname(savename))
        tar.close()
        # Remove the tar file since it has been extracted
        os.remove(savename)

    # Now ICDR v201812.0.0 for 2020
    savename = os.path.join(
        savedir,
        f"cds-satellite-soil-moisture_icdr_{args.frequency}_v201812.0.0.tar")

    client.retrieve(
        'satellite-soil-moisture', {
            'format': 'tgz',
            'variable': 'volumetric_surface_soil_moisture',
            'type_of_sensor': 'combined_passive_and_active',
            'year': [
                '2020',
            ],
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
                '10',
                '11',
                '12',
            ],
            'version': 'v201812.0.0',
            **freq_specific_kwargs[args.frequency]
        }, savename)
    # Unpack the file
    tar = tarfile.open(savename)
    tar.extractall(path=os.path.dirname(savename))
    tar.close()
    # Remove the tar file since it has been extracted
    os.remove(savename)


if __name__ == "__main__":
    download_cds_satellite_sm()
