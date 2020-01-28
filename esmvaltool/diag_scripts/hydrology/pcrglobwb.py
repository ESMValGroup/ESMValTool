"""PCR-GLOBWB diagnostic."""
import logging
from pathlib import Path

import dask.array as da
import iris
import numpy as np
import xarray as xr

from esmvalcore.cmor.table import CMOR_TABLES
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            select_metadata,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_file):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the PCR-GLOBWB hydrological model.",
        'domains': ['global'],
        'authors': [
            'aerts_jerom',
            'andela_bouwe',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': [ancestor_file],
    }
    return record

def add_spinup_year(cube, cube_climatology, varname):

    # Remove leap year day from climatology
    cube_climatology = cube_climatology.extract(iris.Constraint(day_of_year=lambda cell: cell<366))
    
    # Set climatology year in front of regular startyear
    points = cube.coord('time').points[0] - 366 + cube_climatology.coord('day_of_year').points
    time = cube.coord('time').copy(points)

    # Drop dimension day_of_year
    iris.util.demote_dim_coord_to_aux_coord(cube_climatology, 'day_of_year')
    cube_climatology.remove_coord('day_of_year')

    # Add dimension time
    cube_climatology.add_dim_coord(time, 0)

    # Round times to integer number of days
    time_coord = cube_climatology.coord('time')
    time_coord.points = da.floor(time_coord.core_points())
    time_coord.bounds = None
    
    # Set cube cell_methods to None
    cube.cell_methods = ()
    cube_climatology.cell_methods = ()

    # Set dtype
    cube.data = cube.core_data().astype('float32')
    cube_climatology.data = cube_climatology.core_data().astype('float32')
   

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float32')
        coord.bounds = None
        coord.guess_bounds()
        coord_climatology = cube_climatology.coord(coord_name)
        coord_climatology.points = coord_climatology.core_points().astype('float32')
        coord_climatology.bounds = None
        coord_climatology.guess_bounds()


    # Create CubeList and concatenate
    cube_list = iris.cube.CubeList([cube, cube_climatology])
    new_cube = iris.cube.CubeList(cube_list).concatenate_cube()
    
    return new_cube

def main(cfg):
    """Process data for use as input to the PCR-GLOBWB hydrological model."""
    input_data = cfg['input_data'].values()

    # Loop over variables
    for short_name in "pr", "tas":

        # Select and load in cube regular variable timeseries
        metadata = select_metadata(input_data, variable_group=short_name)[0]
        input_file = metadata['filename']
        cube = iris.load_cube(input_file)

        # Round times to integer number of days
        time_coord = cube.coord('time')
        time_coord.points = da.floor(time_coord.core_points())
        time_coord.bounds = None
        time_coord.guess_bounds()

        # Select and load in cube climatology variable timeseries
        metadata_climatology = select_metadata(input_data, variable_group=short_name+'_climatology')[0]
        input_climatology_file = metadata_climatology['filename']
        cube_climatology = iris.load_cube(input_climatology_file)
        
        # Run function to add spinup year to regular variable timeseries
        cube = add_spinup_year(cube, cube_climatology, short_name)

        # Set lat from highest to lowest value
        cube = cube[:, ::-1, ...]

        # Unit conversion 'kg m-3 day-1' to 'm' precip
        if short_name == "pr":
            cube.units = cube.units / 'kg m-3 day-1'
            cube.data = cube.core_data() / 1000

        # Save data
        output_file = get_diagnostic_filename(
            Path(input_file).stem + '_pcrglobwb', cfg)
        iris.save(cube, output_file, fill_value=1.e20)


        print(cube)
        # Store provenance
        provenance_record = get_provenance_record([input_file,input_climatology_file])
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)