"""PCR-GLOBWB diagnostic."""
import logging
from pathlib import Path

import dask.array as da
import iris

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic,
                                            select_metadata)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption':
        "Forcings for the PCR-GLOBWB hydrological model.",
        'domains': ['global'],
        'authors': [
            'aerts_jerom',
            'andela_bouwe',
            'alidoost_sarah',
            'kalverla_peter',
        ],
        'projects': [
            'ewatercycle',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def add_spinup_year(cube, cube_climatology):
    """Prepend the climatology to the cube.

    To reach the equilibrium, the model was spun up using
    the average climatological forcing over each year.
    """
    # Remove leap year day from climatology
    cube_climatology = cube_climatology.extract(
        iris.Constraint(day_of_year=lambda cell: cell < 366))

    # Set climatology year in front of regular startyear
    points = cube.coord('time').points[0] - 366
    points += cube_climatology.coord('day_of_year').points
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

    # fix dtype
    cube.data = cube.core_data().astype('float32')
    cube_climatology.data = cube_climatology.core_data().astype('float32')

    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float32')
        coord.bounds = None
        coord.guess_bounds()
        coord_climatology = cube_climatology.coord(coord_name)
        coord_climatology.points = (
            coord_climatology.core_points().astype('float32'))
        coord_climatology.bounds = None
        coord_climatology.guess_bounds()

    # Create CubeList and concatenate
    cube_list = iris.cube.CubeList([cube, cube_climatology])
    new_cube = iris.cube.CubeList(cube_list).concatenate_cube()

    return new_cube


def main(cfg):
    """Process data for use as input to the PCR-GLOBWB hydrological model."""
    for dataset, metadata in group_metadata(cfg['input_data'].values(),
                                            'dataset').items():
        for short_name in "pr", "tas":
            logger.info("Processing variable %s for dataset %s", short_name,
                        dataset)

            # Load preprocessed cubes for normal data and climatology
            var = select_metadata(metadata, variable_group=short_name)[0]
            cube = iris.load_cube(var['filename'])
            var_climatology = select_metadata(
                metadata,
                variable_group=short_name + '_climatology',
            )[0]
            cube_climatology = iris.load_cube(var_climatology['filename'])

            # Create a spin-up year for pcrglob based on the climatology data
            cube = add_spinup_year(cube, cube_climatology)

            # Round times to integer number of days
            time_coord = cube.coord('time')
            time_coord.points = da.floor(time_coord.core_points())
            time_coord.bounds = None
            time_coord.guess_bounds()

            # Set lat from highest to lowest value
            cube = cube[:, ::-1, ...]

            # Workaround for bug in PCRGlob
            # (see https://github.com/UU-Hydro/PCR-GLOBWB_model/pull/13)
            for coord_name in ['latitude', 'longitude']:
                coord = cube.coord(coord_name)
                coord.points = coord.points + 0.001

            # Unit conversion 'kg m-3 day-1' to 'm' precip (divide by density)
            if short_name == "pr":
                cube.units = cube.units / 'kg m-3 day-1'
                cube.data = cube.core_data() / 1000

            # Save data
            basename = '_'.join([
                'pcrglobwb',
                Path(var['filename']).stem,
                cfg['basin'],
            ])
            output_file = get_diagnostic_filename(basename, cfg)
            iris.save(cube, output_file, fill_value=1.e20)

            # Store provenance
            provenance_record = get_provenance_record(
                [var['filename'], var_climatology['filename']])
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
