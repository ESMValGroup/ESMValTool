"""LISFLOOD diagnostic."""
import logging
from pathlib import Path

# import dask.array as da
import iris

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata, run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def get_provenance_record(ancestor_files):
    """Create a provenance record."""
    record = {
        'caption': "Forcings for the LISFLOOD hydrological model.",
        'domains': ['global'],
        'authors': [
            'verhoeven_stefan',
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


def get_input_cubes(metadata):
    """Create a dict with all (preprocessed) input files."""
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
        cube.attributes.clear()
        all_vars[short_name] = cube
        provenance['ancestors'].append(filename)

    return all_vars, provenance


def shift_era5_time_coordinate(cube, shift=30):
    """Shift instantaneous variables (default = 30 minutes forward in time).

    After this shift, as an example:
    time format [1990, 1, 1, 11, 30, 0] will be [1990, 1, 1, 12, 0, 0].
    For aggregated variables, already time format is [1990, 1, 1, 12, 0, 0].
    """
    time = cube.coord(axis='T')
    time.points = time.points + shift / (24 * 60)
    time.bounds = None
    time.guess_bounds()
    return cube


def compute_vapour_pressure(tdps):
    """Compute vapour pressure using tetens formula."""
    e0 = 6.11
    tdps_c = d2m - 273.15
    e = e0 * np.exp(7.5*(tdps_c) / (237.3 + (tdps_c)))
    return e


def windspeed(u, v):
    return (u**2+v**2)**.5


def main(cfg):
    """Process data for use as input to the LISFLOOD hydrological model """
    input_metadata = cfg['input_data'].values()
    logger.info(input_metadata)

    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)

        if dataset == 'ERA5':
            shift_era5_time_coordinate(all_vars['tas'])
            shift_era5_time_coordinate(all_vars['tdps'])
            shift_era5_time_coordinate(all_vars['uas'])
            shift_era5_time_coordinate(all_vars['vas'])

        pr = all_vars['pr']
        tas = all_vars['tas']
        tasmax = all_vars['tasmax']
        tasmin = all_vars['tasmin']
        rsds = all_vars['rsds']

        # Compute additional variables as input for lisvap
        tdps = all_vars['tdps']
        uas = all_vars['uas']
        vas = all_vars['vas']
        windspeed_inputs = ['uas', 'vas']
        e_inputs = ['tdps']
        e = compute_vapour_pressure(tdps)
        wspd = compute_windspeed(uas, vas)

        for outvar in (pr, tas, tasmax, tasmin, rsds, e, wspd):
            # Target output file
            # In the cdo example, the output format was:
            # 2m_temperature_final_1990-1996.nc

            # Save data
            time_coord = outvar.coord('time')
            start_year = time_coord.cell(0).point.year
            end_year = time_coord.cell(-1).point.year
            var_name = outvar.var_name
            basename = '_'.join([
                'lisflood',
                var_name,
                cfg['catchment'],
                str(start_year),
                str(end_year),
            ])
            output_file = get_diagnostic_filename(basename, cfg)
            iris.save(outvar, output_file, fill_value=1.e20)

            # Store provenance
            provenance_record = get_provenance_record(
                [outvar['filename']]
            )
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(output_file, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
