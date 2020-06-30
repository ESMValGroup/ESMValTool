"""globwat diagnostic."""
import logging
from pathlib import Path

import iris
from iris.pandas import as_data_frame 
import dask.array as da
import numpy
import pandas as pd
import calendar

import subprocess
from esmvalcore import preprocessor as preproc
from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic,
                                            select_metadata)


logger = logging.getLogger(Path('__file__').name) #for test in python use '__file__'


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

def change_data_type(cube):
    """GlobWat input data types are float32.
    """
    # fix dtype
    cube.data = cube.core_data().astype('float32')

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float32')
        coord.bounds = None
        coord.guess_bounds()
    
    return cube

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
        cube = iris.load_cube(filename)
        #set filling values to -9999
        cube.data.set_fill_value(-9999) 
        #change cube dtype to float32
        all_vars[short_name] = change_data_type(cube)
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

def convert_units_monthly(all_vars):
    """convert pr and pet units from'kg m-2 s-1' to 'mm month-1' ."""
    # pr = all_vars['pr']
    # pr.units = pr.units / 'kg m-3'
    # pr.data = pr.core_data() / 1000.
    # pr.convert_units('mm month-1')

    if 'evspsblpot' in all_vars:
            pet = all_vars['evspsblpot']
    else:
            logger.info("Potential evapotransporation not available, deriving")
            pet = debruin_pet(
                psl=all_vars['psl'],
                rsds=all_vars['rsds'],
                rsdt=all_vars['rsdt'],
                tas=all_vars['tas'],
            )
    pet.units = pet.units / 'kg m-3'
    pet.data = pet.core_data() /1000.
    pet.convert_units('mm month-1')

    # return pr , pet
    return pet


def convert_units_daily(all_vars):
    """convert pr and pet units from'kg m-2 s-1' to 'mm day-1' ."""
    # pr = all_vars['pr']
    # pr.units = pr.units / 'kg m-3'
    # pr.data = pr.core_data() / 1000.
    # pr.convert_units('mm day-1')

    if 'evspsblpot' in all_vars:
            pet = all_vars['evspsblpot']
    else:
            logger.info("Potential evapotransporation not available, deriving")
            pet = debruin_pet(
                psl=all_vars['psl'],
                rsds=all_vars['rsds'],
                rsdt=all_vars['rsdt'],
                tas=all_vars['tas'],
            )
    pet.units = pet.units / 'kg m-3'
    pet.data = pet.core_data() /1000.
    pet.convert_units('mm day-1')

    # return pr , pet
    return pet

def get_month_number_for_output_name(cube): 
    """get month number for output name."""
    # get start year and end year  
    coord_time = cube.coord('time') 
    start_year = coord_time.cell(0).point.year 
    end_year = coord_time.cell(-1).point.year 
    year = start_year - 1   
    #looping over time dimention, get day and months 
    for i in range(0, len(cube.coord('time').points)): 
        n_month = str(coord_time.cell(i).point.month).zfill(2) 
        yield n_month 

# def get_month_day_number_for_output_name(cube): 
#     """get month and day number for output name."""
#     # get start year and end year  
#     coord_time = cube.coord('time') 
#     start_year = coord_time.cell(0).point.year 
#     end_year = coord_time.cell(-1).point.year 
#     year = start_year - 1   
#     #looping over time dimention, get day and months 
#     for i in range(0, len(cube.coord('time').points)): 
#         n_month = str(coord_time.cell(i).point.month).zfill(2) 
#         nday = calendar.monthrange(year,int(n_month))[1] 
#         for daynumber in range(1, nday+1):  
#             n_day_month = str(n_month) + str(daynumber).zfill(2) 
#             yield n_day_month #or it is better to write yield n_month here as well????

def get_output_stem_monthly(metadata,cube):
    """Get output file stem, specific to HYPE."""
    short_to_stem = dict(pr="prc",
                         evspsblpot="eto"
                         )
# TODO: add the dataset and year to the directory name.
# Outputs name are the same for all years. 
# So they have to be saved in different folders based on years and dataset.  
    # for nyear in range (start_year , end_year+1):
    #     if dataset == "ERA5":
    #         dir_name = 'ERA5_' + str(nyear) 
    for i in range(0, len(metadata)): 
        shortname = metadata[i]['short_name']
        for n_month in get_month_number_for_output_name(cube):
            if shortname in short_to_stem:
                stem = Path(short_to_stem[shortname] + str(n_month) + 'wb')
            else:
                stem = Path(metadata[i]['filename']).stem + '_gloabwat'
    print(stem)

        # stem = metadata[i]['alias'] / stem

    return stem

# def get_output_stem_daily(metadata):
#     """Get output file stem, specific to HYPE."""
#     short_to_stem = dict(pr="prc",
#                          evspsblpot="eto"
#                          )
#     for i in range(0, len(metadata)): 
#         shortname = metadata[i]['short_name']
#         for n_day_month in get_month_day_number_for_output_name():
#             if shortname in short_to_stem:
#                 stem = Path(short_to_stem[shortname] + str(n_day_month) + 'wb')
#             else:
#                 stem = Path(metadata[i]['filename']).stem + '_gloabwat'

#         # stem = metadata[i]['alias'] / stem

#     return stem

def main(cfg):
    """Process data for use as input to the GlobWat hydrological model.

    These variables are needed in all_vars:
    evspsblpot (potential_evapotranspiration_flux)
    pr (precipitation_flux)
    psl (air_pressure_at_mean_sea_level)
    rsds (surface_downwelling_shortwave_flux_in_air)
    rsdt (toa_incoming_shortwave_flux)
    tas (air_temperature)
    """
    input_metadata = cfg['input_data'].values()
    # for checking the code in ipython I added print(input_metadata).
    # Run the script and use print(input_metadata) in the log.txt as input_metadata
    print(input_metadata)
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)

        # In ERA5 potential evaporation dataset, positive values represent condensation
        # and negative values represent evaporation. So, ERA5 is multiplied by -1 to 
        # have positive values for evaporation 
        if dataset == 'ERA5':
            all_vars['evspsblpot'].data *= -1  
        
        # extract pr cube from all_vars
        pr = all_vars['pr']
        cube = all_vars['pr']  


        # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm month-1'
        for mip in 'Amon', 'day':
            if mip == 'Amon':
                convert_units_monthly(all_vars)
            # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm day-1'
            elif mip == 'day':
                convert_units_daily(all_vars)

        output_file = get_diagnostic_filename(get_output_stem_monthly(metadata,cube),
                                                  cfg, 'asc')
        Path(output_file).parent.mkdir(exist_ok=True)
        
        for key in all_vars:
            key_cube = all_vars[key]
            for sub_cube in key_cube.slices_over('time'): 
                frame = iris.pandas.as_data_frame(sub_cube, copy=True) 
                frame.to_csv(output_file,
                            sep=' ',
                            na_rep='-9999',
                            index_label="DATE", 
                            float_format='%.1f'
                            ) 
            # # Store provenance
            # with ProvenanceLogger(cfg) as provenance_logger:
            #     provenance_logger.log(output_file, provenance)



if __name__ == '__main__':
   with run_diagnostic() as config:
       main(config)