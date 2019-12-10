import cdsapi
from datetime import datetime
from dateutil import relativedelta
import argparse
import os
import yaml

def download_cds_satellite_albedo():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config_file', '-c',
                        default=os.path.join(os.path.dirname(__file__),
                                             'config-user.yml'),
                        help='Config file')
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


    c = cdsapi.Client()
    
    start = '199901'
    end = '201401'
    loopdate = datetime.strptime(start,'%Y%m')
    enddate = datetime.strptime(end,'%Y%m')
    
    while loopdate < enddate:
        
        loopdate += relativedelta.relativedelta(months=1)
        month = loopdate.month
        year = loopdate.year
        savename = os.path.join(cds_satellite_albedo_dir,
                                f"cds-satellite-albedo_{year}{month:02d}.tar.gz")
        c.retrieve(
            'satellite-albedo',
            {
                'format': 'tgz',
                'satellite': 'spot',
                'sensor': 'vgt',
                'product_version': 'V1',
                'horizontal_resolution': '1km',
                'variable': [
                    'albb_bh', 'albb_dh',
                ],
                'year': f"{year}",
                'month': [f"{month:02d}"],
                'nominal_day': [
                    '10', '20', '28',
                    '29', '30', '31',
                ],
            },
            savename)


if __name__ == "__main__":
    download_cds_satellite_albedo()
