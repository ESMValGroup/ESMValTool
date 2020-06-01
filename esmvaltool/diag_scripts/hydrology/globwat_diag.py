"""globwat diagnostic."""
import logging
from pathlib import Path

import iris
# import subprocess

from esmvalcore import preprocessor as preproc
from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def create_provenance_record():
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the GlobWat hydrological model.",
        'domains': ['global'],
        'authors': [
            'abdollahi_banafsheh',
            'alidoost_sarah',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [],
    }
    return record


def get_input_cubes(metadata):
    """Return a dict with all (preprocessed) input files."""
    provenance = create_provenance_record()
    all_vars = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        if short_name in all_vars:
            raise ValueError(
                f"Multiple input files found for variable '{short_name}'.")
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        all_vars[short_name] = iris.load_cube(filename)
        provenance['ancestors'].append(filename)
    return all_vars, provenance


def _get_extra_info(cube):
    """Get start/end time and origin of cube.

    Get the start and end times as an array with length 6
    and get latitude and longitude as an array with length 2
    """
    coord = cube.coord('time')
    time_start_end = []
    for index in 0, -1:
        time_val = coord.cell(index).point
        time_val = [
            time_val.year,
            time_val.month,
            time_val.day,
            time_val.hour,
            time_val.minute,
            time_val.second,
        ]
        time_val = [float(time) for time in time_val]
        time_start_end.append(time_val)

    # Add data_origin
    lat_lon = [
        cube.coord(name).points[0] for name in ('latitude', 'longitude')
    ]
    return time_start_end, lat_lon


def _shift_era5_time_coordinate(cube):
    """Shift instantaneous variables 30 minutes forward in time.

    After this shift, as an example:
    time format [1990, 1, 1, 11, 30, 0] will be [1990, 1, 1, 12, 0, 0].
    For aggregated variables, already time format is [1990, 1, 1, 12, 0, 0].
    """
    time = cube.coord(axis='T')
    time.points = time.points + 30 / (24 * 60)
    time.bounds = None
    time.guess_bounds()
    return cube


def main(cfg):
    """Process data for use as input to the GlobWat hydrological model.

    These variables are needed in all_vars:
    tas (air_temperature)
    pr (precipitation_flux)
    psl (air_pressure_at_mean_sea_level)
    rsds (surface_downwelling_shortwave_flux_in_air)
    rsdt (toa_incoming_shortwave_flux)
    """
    input_metadata = cfg['input_data'].values()
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)

        
        # Processing variables and unit conversion
        # Unit of the fluxes in GlobWat should be in mm
        logger.info("Processing variable PET")
        pet = debruin_pet(
            psl=all_vars['psl'],
            rsds=all_vars['rsds'],
            rsdt=all_vars['rsdt'],
            tas=all_vars['tas'],
        )
        # Unit conversion 'kg m-3 s-1' to 'mm' precip (multiplied by second)
        for short_name == "pr":
            if mip == "Amon":
                cube.data = cube.core_data() * 2592000
            if mip == "day"
                cube.data = cube.core_data() * 86400

        # Unit conversion 'kg m-3 s-1' to 'mm' evspsblpot (multiplied by second)
        for short_name == "evspsblpot":
            if dataset == 'ERA5':
        # to make potential evaporation by ERA5 positive,multiplied by -1
            cube.data = cube.core_data() * -2592000

        # Unit conversion 'kg m-3 s-1' to 'mm' pet calculated &
        # by debruin from ERA-Interim(multiplied by second)
            else:
            cube.data = cube.core_data() * 2592000                  

        # get start year and end year 
        coord_time = pet.coord('time')
        start_year = coord_time.cell(0).point.year
        end_year = coord_time.cell(-1).point.year     

        # Save to netcdf file 
        basename = '_'.join([
            'globwat_input',
            dataset,
            str(start_year),  
            str(end_year), 
        ])

        # output_data = {
        #     'forcing': {
        #         'precip': precip.data,
        #         'temp': temp.data,
        #         'pet': pet.data,
        #         'delta_t_days': float(1),
        #         'time_unit': 'day',
        #     },
        #     'time_start': time_start_end[0],
        #     'time_end': time_start_end[1],
        #     'data_origin': lat_lon,
        # }        
        
        output_name = get_diagnostic_filename(basename, cfg, extension='nc')
        iris.save(pet, output_name, fill_value=-9999)

        #  # slice data and save in ascii format
        # if short_name == "pr":
        #     for sub_data in data.slices_over('time'):
        #     subprocess.call('gdal_translate -of AAIGrid NETCDF:"sub_data":pr sub_data.asc')
        # if short_name == "evspsblpot":
        #     for sub_data in data.slices_over('time'):
        #     subprocess.call('gdal_translate -of AAIGrid NETCDF:"sub_data":evspsblpot sub_data.asc')

        # Store provenance
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(output_name, provenance)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
