"""globwat diagnostic."""
import logging
from pathlib import Path

import iris
import calendar
import dask.array as da
import numpy
import pandas as pd 
import os
import os.path

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

def dtype_change(cube):
    """GlobWat input types are float32.
    """
    # fix dtype
    cube.data = cube.core_data().astype('float32')

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float32')
        coord.bounds = None
        coord.guess_bounds()
    
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
    # for checking the code in ipython I added print(input_metadata).
    # Run the script and use print(input_metadata) in the log.txt as input_metadata
    print(input_metadata)
    # for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
    #     all_vars, provenance = get_input_cubes(metadata)

    for dataset, metadata in group_metadata(input_metadata,
                                            'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)
        for short_name in "pr", "tas":
            logger.info("Processing variable %s for dataset %s", short_name,
                        dataset)

            # Load preprocessed cubes for normal data and climatology
            var = select_metadata(metadata, variable_group=short_name)[0]
            cube = iris.load_cube(var['filename'])

        #change cube dtype to float32
            cube = dtype_change(cube)

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
        # get start year and end year 
        coord_time = cube.coord('time')


        # coord_time = pet.coord('time')
        start_year = coord_time.cell(0).point.year
        end_year = coord_time.cell(-1).point.year 
        print(start_year)
        print(end_year)
        
        # TODO: check the unit conversion by the preprocessor 
        # TODO: if needed make a function for unit conversion 
        year = start_year - 1
        for j in range (end_year - start_year + 1):
            year = year + 1
            for i in range (1, 12):
                nday = calendar.monthrange(year,i)[1]
                for short_name in "pr", "evspsblpot":
                    logger.info("Processing variable %s for dataset %s", short_name,
                                dataset)
                    for mip in "Amon", "day":
                        logger.info("Processing variable %s for dataset %s", mip,
                                    dataset)

                        # Load preprocessed cubes for normal data 
                        # print(len(select_metadata(metadata, variable_group=short_name)))
                        # var = select_metadata(metadata, variable_group=short_name)[0]
                        # cube = iris.load_cube(var['filename'])
                        # cube = dtype_change(cube)

                        # print('cube.name', cube.name)
                        # print("cube banafsheh=",cube)
                        if short_name == "pr" and mip == "Amon":
                            cube.units = cube.units / 'mm'
                            cube.data = cube.core_data() * nday * 86400
                        elif short_name == "evspsblpot":
                            # to make potential evaporation by ERA5 positive,multiplied by -1
                            if dataset == "ERA5":
                                cube.units = cube.units / 'mm'
                                cube.data = cube.core_data() * nday * -86400
                            else:
                                cube.units = cube.units / 'mm'
                                cube.data = cube.core_data() * nday * 86400
            
                        if short_name == "pr" and mip == "day":
                                cube.units = cube.units / 'mm'
                                cube.data = cube.core_data() * 86400
                            
        # Load preprocessed cubes for normal data 
        # var = select_metadata(metadata, variable_group=short_name)[0]
        # cube = iris.load_cube(var['filename'])
        print("cube shape=",cube.shape)
        print("cube ndim=",cube.ndim)


        #looping over time dimention, get time and data for each time
        for i in range(0, len(cube.coord('time').points)):
            n_month = str(coord_time.cell(i).point.month).zfill(2)
            nday = calendar.monthrange(year,int(n_month))[1]
            for daynumber in range(1, nday+1): 
                n_day_month = str(n_month) + str(daynumber).zfill(2) 

            # df = pd.DataFrame(cube.data[i])
            # data_slice = df.replace(numpy.nan, -9999)

        for sub_cube in cube.slices_over('time'):
            for i in range(0, len(cube.coord('time').points)):
                data_slice = sub_cube.data[i]
                df = pd.DataFrame(data_slice)
                data_ch_nan = df.replace(numpy.nan, -9999)

        #     df = pd.DataFrame(sub_cube)
        #     data_slice = df.replace(numpy.nan, -9999)
            # make data structure
                for nyear in range (start_year , end_year+1):
                    if dataset == "ERA5":
                        dir_name = 'ERA5_' + str(nyear) 
                        os.mkdir(dir_name)
                        if mip == 'Amon':
                            for var in range (0, len(data_slice)):
                                if data_slice[var] == 'pr':
                                    output = pr.data_slice
                                    basename = '_'.join([
                                        'prc',
                                        n_month,
                                        'wb' 
                                    ])
                                    output_name = get_diagnostic_filename(basename, cfg, extension='asc')
                                    file_path = os.path.join(dir_name, output_name)
                                    output.to_csv(basename)

                                    if all_vars[var] == 'evspsblpot':
                                        output = evspsblpot.data_slice
                                        basename = '_'.join([
                                            'eto',
                                            n_month,
                                            'wb' 
                                        ])
                                        output_name = get_diagnostic_filename(basename, cfg, extension='asc')
                                        file_path = os.path.join(dir_name, output_name)
                                        output.to_csv(basename)
                            else:
                                for var in range (0, len(data_slice)):
                                    if data_slice[var] == 'pr':
                                        output = pr.data_slice
                                        basename = '_'.join([
                                            'prc',
                                            n_day_month,
                                            'wb' 
                                        ])
                                        output_name = get_diagnostic_filename(basename, cfg, extension='asc')
                                        file_path = os.path.join(dir_name, output_name)
                                        output.to_csv(basename)

                        if dataset == 'ERA-Interim':
                            dir_name = 'ERA-Interim_' + str(nyear) 
                            os.mkdir(dir_name)
                            if mip == 'Amon':
                                for var in range (0, len(data_slice)):
                                    if data_slice[var] == 'pr':
                                        output = pr.data_slice
                                        basename = '_'.join([
                                            'prc',
                                            n_month,
                                            'wb' 
                                        ])
                                        output_name = get_diagnostic_filename(basename, cfg, extension='asc')
                                        file_path = os.path.join(dir_name, output_name)
                                        output.to_csv(basename)
                                    else:
                                        output = pet
                                        basename = '_'.join([
                                            'eto',
                                            n_month,
                                            'wb' 
                                        ])
                                        output_name = get_diagnostic_filename(basename, cfg, extension='asc')
                                        file_path = os.path.join(dir_name, output_name)
                                        output.to_csv(basename)
                            else:
                                for var in range (0, len(data_slice)):
                                    if data_slice[var] == 'pr':
                                        output = pr.data_slice
                                        basename = '_'.join([
                                            'prc',
                                            n_day_month,
                                            'wb' 
                                        ])
                                        output_name = get_diagnostic_filename(basename, cfg, extension='asc')
                                        file_path = os.path.join(dir_name, output_name)
                                        output.to_csv(basename)

        # # Store provenance
        # with ProvenanceLogger(cfg) as provenance_logger:
        #     provenance_logger.log(output_name, provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)