#!/usr/bin/env python3
import argparse
import os

import yaml
from ecmwfapi import ECMWFDataServer

parser = argparse.ArgumentParser(description='''
Script to download era-interim data, so it can be cmorized and used by esmvaltool recipes.

Before running the script:
1. Install the dependency i.e. ECMWFDataServer. For this, run pip install ecmwf-api-client
2. Create an account at https://www.ecmwf.int/ 
3. Follow the instruction https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch
4. Copy/paste the text in https://api.ecmwf.int/v1/key/ into a blank text file and save it as $HOME/.ecmwfapirc
5. Use esmvaltool/config-user.yml as an template and set the rootpath of the output directory in RAWOBS

```bash
python download_era-interim.py -c config-user.yml
cmorize_obs -c config-user.yml
```
''')
parser.add_argument('--config_file', '-c',
                    default=os.path.join(os.path.dirname(__file__),
                                         'config-user.yml'),
                    help='Config file')
parser.add_argument('--start_year', type=int, default=1979, help='Start year')
parser.add_argument('--end_year', type=int, default=2019, help='End year')
args = parser.parse_args()

# get and read config file
config_file = os.path.abspath(
    os.path.expandvars(os.path.expanduser(args.config_file)))

with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

rawobs_dir = os.path.abspath(
    os.path.expandvars(os.path.expanduser(config['rootpath']['RAWOBS'])))
erainterim_dir = f'{rawobs_dir}/Tier3/ERA-Interim'
os.makedirs(erainterim_dir, exist_ok=True)

years = range(args.start_year, args.end_year + 1)
server = ECMWFDataServer()

day_timesteps = {
    'fc': {
        'step': '3/6/9/12',
        'time': '00:00:00/12:00:00',
        'type': 'fc',
    },
    'accu': {
        'step': '12',
        'time': '00:00:00/12:00:00',
        'type': 'fc',
    },
    'an': {
        'type': 'an',
        'time': '00:00:00/06:00:00/12:00:00/18:00:00',
        'step': '0',
    }
}

day_params = [
    ('165.128', 'u10', 'an'),  #10 metre U wind component
    ('166.128', 'v10', 'an'),  #10 metre V wind component
    ('167.128', 't2m', 'an'),  #2 metre temperature
    ('168.128', 'd2m', 'an'),  #2 metre dewpoint temperature
    ('169.128', 'ssrd', 'accu'),  #Surface solar radiation downwards
    ('201.128', 'mx2t', 'fc'),  #Max. temp at 2m since previous post-processing
    ('202.128', 'mn2t', 'fc'),  #Min. temp at 2m since previous post-processing
    ('228.128', 'tp', 'accu'),  #Total precipitation
    ('151.128', 'msl', 'an'),  #Mean sea level pressure
    ('182.128', 'e', 'accu'),  #Evaporation
    # TODO not found in era-interim yet
    # ('', 'pev'),  # Potential evaporation
]
for no, symbol, timestep in day_params:
    fr = 'daily'
    for year in years:
        server.retrieve({
            'class': 'ei',
            'dataset': 'interim',
            'date': f'{year}-01-01/to/{year}-12-31',
            'expver': '1',
            'grid': '0.75/0.75',
            'levtype': 'sfc',
            'param': no,
            'stream': 'oper',
            'format': 'netcdf',
            'target': f'{erainterim_dir}/ERA-Interim_{symbol}_{fr}_{year}.nc',
            **day_timesteps[timestep]
        })

month_timesteps = {
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
    }
}

month_params = [
    ('228.128', 'tp', 'accu'),  #Total precipitation
    ('164.128', 'tcc', 'an'),  #Total cloud cover
    ('56.162', 'p56.162', 'an'),  #Vertical integral of cloud liquid water
    ('57.162', 'p57.162', 'an'),  #Vertical integral of cloud frozen water
    ('137.128', 'tcwv', 'an'),  #Total column water vapour
    ('134.128', 'sp', 'an'),  #Surface pressure
    ('151.128', 'msl', 'an'),  #Mean sea level pressure
    ('167.128', 't2m', 'an'),  #2 metre temperature
    ('229.128', 'iews', 'fc'),  #Inst. eastward turbulent surface stress
    ('230.128', 'inss', 'fc'),  #Inst. northward turbulent surface stress
    ('34.128', 'sst', 'an'),  #Sea surface temperature
    ('235.128', 'skt', 'an'),  #Skin temperature
    ('176.128', 'ssr', 'accu'),  #Surface net solar radiation
    ('177.128', 'str', 'accu'),  #Surface net thermal radiation
    ('147.128', 'slhf', 'accu'),  #Surface latent heat flux
    ('146.128', 'sshf', 'accu'),  #Surface sensible heat flux
    # TODO not found in era-interim yet
    # ('157.128', 'r', ''),  #Relative humidity
    # ('130.128', 't', ''),  #
    # ('131.128', 'u', ''),  #
    # ('132.128', 'v', ''),  #
    # ('135.128', 'w', ''),  #
    # ('133.128', 'q', ''),  #
]

for no, symbol, timestep in month_params:
    fr = 'monthly'
    for year in years:
        server.retrieve({
            'class': 'ei',
            'dataset': 'interim',
            #All months of a year eg. 19900101/19900201/.../19901101/19901201
            'date': '/'.join([f'{year}{m:02}01' for m in range(1, 13)]),
            'expver': '1',
            'grid': '0.75/0.75',
            'param': no,
            'format': 'netcdf',
            'target': f'{erainterim_dir}/ERA-Interim_{symbol}_{fr}_{year}.nc',
            **month_timesteps[timestep]
        })

invariant_params = [
    ('172.128', 'lsm'),  # Land-sea mask
    ('129.128', 'z'),  # Geopotential
]

for no, symbol in invariant_params:
    server.retrieve({
        'class': 'ei',
        'dataset': 'interim',
        'date': '1989-01-01',
        'expver': '1',
        'grid': '0.75/0.75',
        'levtype': 'sfc',
        'param': no,
        'step': '0',
        'stream': 'oper',
        'time': '12:00:00',
        'type': 'an',
        'format': 'netcdf',
        'target': f'{erainterim_dir}/ERA-Interim_{symbol}.nc',
    })
