"""Script to download era-interim data.

Before running the script:
1. Install the dependency i.e. ECMWFDataServer.
For this, run pip install ecmwf-api-client

2. Create an account at https://www.ecmwf.int/

3. Follow the instruction at:
https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch

4. Copy/paste the text in https://api.ecmwf.int/v1/key/ into a blank text file
and save it as $HOME/.ecmwfapirc

5. Use ESMValCore/esmvalcore/config-user.yml as an template
and set the rootpath of the output directory in RAWOBS

6. Check the description of the variables at
https://apps.ecmwf.int/codes/grib/param-db

7. Check the invariant variables at
https://apps.ecmwf.int/datasets/data/interim-full-invariant

```bash
python download_era_interim.py --config_file config-user.yml --start_year 2000
--end_year 2000
```

This will download and save the data in the RAWOBS directory,
under Tier3/ERA-Interim.

"""
import argparse
import os

import yaml
from ecmwfapi import ECMWFDataServer


DAY_TIMESTEPS = {
    'fc': {
        'step': '3/6/9/12',
        'time': '00:00:00/12:00:00',
        'type': 'fc',
        'levtype': 'sfc',
    },
    'accu': {
        'step': '12',
        'time': '00:00:00/12:00:00',
        'type': 'fc',
        'levtype': 'sfc',
    },
    'an': {
        'type': 'an',
        'time': '00:00:00/06:00:00/12:00:00/18:00:00',
        'step': '0',
        'levtype': 'sfc',
    },
    '3d': {
        'type': 'an',
        'time': '00:00:00/06:00:00/12:00:00/18:00:00',
        'step': '0',
        'levelist': '10/50/100/250/500/700/850/1000',  # CMIP6 day table, plev8
        'levtype': 'pl',
    }
}


DAY_PARAMS = [
    ('167.128', 't2m', 'an'),  # 2 metre temperature
    ('228.128', 'tp', 'accu'),  # Total precipitation
    ('182.128', 'e', 'accu'),  # Evaporation
    ('201.128', 'mx2t', 'fc'),  # Max. temp at 2m since previous post-proc
    ('202.128', 'mn2t', 'fc'),  # Min. temp at 2m since previous post-proc
    ('235.128', 'skt', 'an'),  # Skin temperature
    ('165.128', 'u10', 'an'),  # 10 metre U wind component
    ('166.128', 'v10', 'an'),  # 10 metre V wind component
    ('168.128', 'd2m', 'an'),  # 2 metre dewpoint temperature
    ('151.128', 'msl', 'an'),  # Mean sea level pressure
    ('134.128', 'sp', 'an'),  # Surface pressure
    ('144.128', 'sf', 'accu'),  # Snowfall
    ('176.128', 'ssr', 'accu'),  # Surface net solar radiation
    ('169.128', 'ssrd', 'accu'),  # Surface solar radiation downwards
    ('175.128', 'strd', 'accu'),  # Surface thermal radiation downwards
    ('205.128', 'ro', 'accu'),  # Runoff
    ('238.128', 'tsn', 'an'),  # Temperature of snow layer
    ('212.128', 'tisr', 'accu'),  # TOA incident solar radiation
    ('164.128', 'tcc', 'an'),  # Total cloud cover
    ('129.128', 'z', '3d'),  # Geopotential
]


MONTH_TIMESTEPS = {
    'accu': {
        'levtype': 'sfc',
        'stream': 'mdfa',
        'type': 'fc',
        'step': '0-12'
    },
    'an': {
        'levtype': 'sfc',
        'stream': 'moda',
        'type': 'an'
    },
    'fc': {
        'levtype': 'sfc',
        'stream': 'moda',
        'type': 'fc'
    },
    '3d': {
        'levtype': 'pl',
        'stream': 'moda',
        'type': 'an',
        'levelist': '1/5/10/20/30/50/70/100/150/200/' +
                    '250/300/400/500/600/700/850' +
                    '/925/1000'  # CMIP6 Amon table, plev19
    }
}


MONTH_PARAMS = [
    ('167.128', 't2m', 'an'),  # 2 metre temperature
    ('228.128', 'tp', 'accu'),  # Total precipitation
    ('182.128', 'e', 'accu'),  # Evaporation
    ('235.128', 'skt', 'an'),  # Skin temperature
    ('165.128', 'u10', 'an'),  # 10 metre U wind component
    ('166.128', 'v10', 'an'),  # 10 metre V wind component
    ('168.128', 'd2m', 'an'),  # 2 metre dewpoint temperature
    ('151.128', 'msl', 'an'),  # Mean sea level pressure
    ('144.128', 'sf', 'accu'),  # Snowfall
    ('176.128', 'ssr', 'accu'),  # Surface net solar radiation
    ('169.128', 'ssrd', 'accu'),  # Surface solar radiation downwards
    ('205.128', 'ro', 'accu'),  # Runoff
    ('238.128', 'tsn', 'an'),  # Temperature of snow layer
    ('212.128', 'tisr', 'accu'),  # TOA incident solar radiation
    ('164.128', 'tcc', 'an'),  # Total cloud cover
    ('56.162', 'p56.162', 'an'),  # Vertical integral of cloud liquid water
    ('57.162', 'p57.162', 'an'),  # Vertical integral of cloud frozen water
    ('137.128', 'tcwv', 'an'),  # Total column water vapour
    ('134.128', 'sp', 'an'),  # Surface pressure
    ('229.128', 'iews', 'fc'),  # Inst. eastward turbulent surface stress
    ('230.128', 'inss', 'fc'),  # Inst. northward turbulent surface stress
    ('34.128', 'sst', 'an'),  # Sea surface temperature
    ('177.128', 'str', 'accu'),  # Surface net thermal radiation
    ('147.128', 'slhf', 'accu'),  # Surface latent heat flux
    ('146.128', 'sshf', 'accu'),  # Surface sensible heat flux
    ('157.128', 'r', '3d'),  # Relative humidity
    ('130.128', 't', '3d'),  # Temperature
    ('131.128', 'u', '3d'),  # U component of wind
    ('132.128', 'v', '3d'),  # V component of wind
    ('135.128', 'w', '3d'),  # Vertical velocity
    ('133.128', 'q', '3d'),  # Specific humidity
    ('129.128', 'z', '3d'),  # Geopotential
]


INVARIANT_PARAMS = [
    ('172.128', 'lsm'),  # Land-sea mask
    ('129.128', 'z'),  # Geopotential (invariant at surface)
]


def _get_daily_data(params, timesteps, years, server, era_interim_dir):
    for param_id, symbol, timestep in params:
        frequency = 'daily'
        for year in years:
            server.retrieve({
                'class': 'ei',
                'dataset': 'interim',
                'date': f'{year}-01-01/to/{year}-12-31',
                'expver': '1',
                'grid': '0.75/0.75',
                'param': param_id,
                'stream': 'oper',
                'format': 'netcdf',
                'target': f'{era_interim_dir}/ERA-Interim_{symbol}'
                          f'_{frequency}_{year}.nc',
                **timesteps[timestep]
            })


def _get_monthly_data(params, timesteps, years, server, era_interim_dir):
    for param_id, symbol, timestep in params:
        frequency = 'monthly'
        for year in years:
            server.retrieve({
                'class': 'ei',
                'dataset': 'interim',
                # All months of a year eg. 19900101/.../19901101/19901201
                'date': '/'.join([f'{year}{m:02}01' for m in range(1, 13)]),
                'expver': '1',
                'grid': '0.75/0.75',
                'param': param_id,
                'format': 'netcdf',
                'target': f'{era_interim_dir}/ERA-Interim_{symbol}'
                          f'_{frequency}_{year}.nc',
                **timesteps[timestep]
            })


def _get_invariant_data(params, server, era_interim_dir):
    for param_id, symbol in params:
        server.retrieve({
            'class': 'ei',
            'dataset': 'interim',
            'date': '1989-01-01',
            'expver': '1',
            'grid': '0.75/0.75',
            'levtype': 'sfc',
            'param': param_id,
            'step': '0',
            'stream': 'oper',
            'time': '12:00:00',
            'type': 'an',
            'format': 'netcdf',
            'target': f'{era_interim_dir}/ERA-Interim_{symbol}.nc',
        })


def cli():
    """Download ERA-Interim variables from ECMWF data server."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config_file', '-c',
                        default=os.path.join(os.path.dirname(__file__),
                                             'config-user.yml'),
                        help='Config file')
    parser.add_argument('--start_year', type=int,
                        default=1979, help='Start year')
    parser.add_argument('--end_year', type=int, default=2019, help='End year')
    args = parser.parse_args()

    # get and read config file
    config_file_name = os.path.abspath(
        os.path.expandvars(os.path.expanduser(args.config_file)))

    with open(config_file_name, 'r') as config_file:
        config = yaml.safe_load(config_file)

    rawobs_dir = os.path.abspath(
        os.path.expandvars(os.path.expanduser(config['rootpath']['RAWOBS'])))
    era_interim_dir = f'{rawobs_dir}/Tier3/ERA-Interim'
    os.makedirs(era_interim_dir, exist_ok=True)

    years = range(args.start_year, args.end_year + 1)
    server = ECMWFDataServer()

    _get_daily_data(DAY_PARAMS, DAY_TIMESTEPS, years, server, era_interim_dir)
    _get_monthly_data(MONTH_PARAMS, MONTH_TIMESTEPS,
                      years, server, era_interim_dir)
    _get_invariant_data(INVARIANT_PARAMS, server, era_interim_dir)


if __name__ == "__main__":
    cli()
