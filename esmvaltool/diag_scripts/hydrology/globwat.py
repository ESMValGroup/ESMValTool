"""Globwat diagnostic."""
import logging
from pathlib import Path

import numpy as np
import xarray as xr
import pandas as pd
import dask.array as da
import iris

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.hydrology.compute_chunks import compute_chunks
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            group_metadata,
                                            run_diagnostic)


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
            'debruin16ams',
            'hoogeveen15hess',
            'langbein1949usgs',
        ],
        'ancestors': [],
    }
    return record


def rechunk_and_regrid(src, tgt, scheme):
    """Rechunk cube src and regrid it onto the grid of cube tgt."""
    src_chunks = compute_chunks(src, tgt)
    src.data = src.lazy_data().rechunk(src_chunks)
    return regrid(src, tgt, scheme)


def change_data_type(cube):
    """Change data type to float32."""
    cube.data = cube.core_data().astype('float32')
    for coord_name in 'latitude', 'longitude', 'time':
        coord = cube.coord(coord_name)
        coord.points = coord.core_points().astype('float32')
        coord.bounds = None
        coord.guess_bounds()
    return cube


def _convert_units(cube):
    """Convert unit of cube, used only for water variables.

    From kg m-2 s-1 to kg m-2 month-1 or kg m-2 day-1.
    Note that the unit kg m-2 s-1 is equivalent to mm s-1.
    """
    mip = cube.attributes['mip']

    if mip == 'Amon':
        cube.convert_units('kg m-2 month-1')  # equivalent to mm/month
    elif mip == 'day':
        cube.convert_units('kg m-2 day-1')  # equivalent to mm/day
    return cube


def _fix_negative_values(cube):
    """Change negative values to zero."""
    cube.data = da.where(cube.core_data() < 0, 0, cube.core_data())
    return cube


def get_input_cubes(metadata):
    """Return a dictionary with all (preprocessed) input files."""
    provenance = create_provenance_record()
    all_vars = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        all_vars[short_name] = change_data_type(cube)
        cube.attributes['mip'] = attributes['mip']
        provenance['ancestors'].append(filename)
    return all_vars, provenance


def load_target(cfg):
    """Load target grid."""
    filename = Path(cfg['auxiliary_data_dir']) / cfg['target_grid_file']
    cube = iris.load_cube(str(filename))
    for coord in 'longitude', 'latitude':
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    return cube


def langbein_pet(tas):
    """Calculate potential ET using Langbein method.

    The Langbein curve represents an empirical relationship between temperature
    and potential ET (pet). Where T is the annual average temperature in degree
    Celcius, pet is in mm per year, and a, b, c are unitless empirical
    constants.
    Reference: https://doi.org/10.3133/cir52 page 8, figure 1.
    An example of using Langbein method can be found at:
    https://doi.org/10.1080/02626667.2017.1332416 page 1472, equation 7.
    """
    tas.convert_units('degC')
    constant_a = iris.coords.AuxCoord(np.float32(325),
                                      long_name='first constant', units=None)
    constant_b = iris.coords.AuxCoord(np.float32(21),
                                      long_name='second constant', units=None)
    constant_c = iris.coords.AuxCoord(np.float32(0.9),
                                      long_name='third constant', units=None)

    # assumption here: tas is constant over time, then the monthly/daily
    # average value is equal to the annual average.
    pet = (tas) * constant_b + (tas ** 2) * constant_c + constant_a
    pet.units = 'kg m-2 year-1'  # equivalent to mm year-1
    pet.convert_units('kg m-2 s-1')  # convert to a cmor compatible unit
    pet.var_name = 'evspsblpot'
    pet.standard_name = 'water_potential_evaporation_flux'
    pet.long_name = 'Potential Evapotranspiration'
    return pet


def get_cube_time_info(cube):
    """Return year, month and day from the cube."""
    coord_time = cube.coord('time')
    time = coord_time.cell(0).point
    time_step = time.strftime("%Y%m%d")
    return time_step


def get_cube_data_info(cube):
    """Return short_name, and mip from the cube."""
    short_name = cube.var_name
    mip = cube.attributes['mip']
    return short_name, mip


def _swap_western_hemisphere(cube):
    """Set longitude values in range -180, 180.

    Western hemisphere longitudes should be negative.
    """
    array = xr.DataArray.from_iris(cube)

    # Set longitude values in range -180, 180.
    array['lon'] = (array['lon'] + 180) % 360 - 180

    # Re-index data along longitude values
    west = array.where(array.lon < 0, drop=True)
    east = array.where(array.lon >= 0, drop=True)
    return west.combine_first(east)


def _flip_latitudes(array):
    """Flip latitudes for writing as ascii.

    Latitudes order should be in range 90, -90.
    """
    flipped = array[::-1, ...]
    flipped['lat'] = array['lat'] * -1
    return flipped


def save_to_ascii(cube, file_name):
    """Save data to an ascii file.

    Data with index [0,0] should be in -180, 90 lon/lat.
    """
    # Re-index data
    array = _swap_western_hemisphere(cube)
    array = _flip_latitudes(array)

    # Set nodata values
    array = array.fillna(-9999)

    xmin = array['lon'].min().values
    ymin = array['lat'].min().values
    xres = array['lon'].values[1] - array['lon'].values[0]
    output = open(file_name, "w")
    output.write(f"ncols {array.shape[1]}\n")
    output.write(f"nrows {array.shape[0]}\n")
    output.write(f"xllcorner     {xmin}\n")
    output.write(f"yllcorner     {ymin}\n")
    output.write(f"cellsize      {xres}\n")
    output.write(f"NODATA_value  {np.int32(-9999)}\n")
    output.close()

    data_frame = pd.DataFrame(array.values, dtype=array.dtype)
    data_frame.to_csv(file_name, sep=' ', na_rep='-9999', float_format=None,
                      header=False, index=False, mode='a')


def make_filename(dataset_name, cfg, cube, extension='asc'):
    """Return a valid path for saving a diagnostic data file.

    filenames are specific to Globwat.
    """
    time_stamp = get_cube_time_info(cube)
    short_name, mip = get_cube_data_info(cube)
    if cfg['evaporation_method'] == 'langbein':
        pet_method_name = 'langbein_'
    else:
        pet_method_name = 'debruin_'

    if short_name == 'pet':
        pet_method = pet_method_name
    else:
        pet_method = ''

    base_name = (f"globwat_{dataset_name}_{mip}_{short_name}_{pet_method}"
                 f"{time_stamp}")
    filename = get_diagnostic_filename(base_name, cfg, extension=extension)
    return filename


def _shift_era5_time_coordinate(cube):
    """Shift instantaneous variables 30 minutes forward in time.

    After this shift, as an example:
    time format [1990, 1, 1, 11, 30, 0] will be [1990, 1, 1, 12, 0, 0].
    For aggregated variables, already time format is [1990, 1, 1, 12, 0, 0].
    """
    if not cube.attributes['mip'] == 'Amon':
        time = cube.coord(axis='T')
        time.points = time.points + 30 / (24 * 60)
        time.bounds = None
        time.guess_bounds()
    return cube


def main(cfg):
    """Process data for GlobWat hydrological model.

    These variables are needed in all_vars:
    pr (precipitation_flux)
    psl (air_pressure_at_mean_sea_level)
    rsds (surface_downwelling_shortwave_flux_in_air)
    rsdt (toa_incoming_shortwave_flux)
    tas (air_temperature)
    """
    # Load target grid to be used in re-gridding
    target_cube = load_target(cfg)

    input_metadata = cfg['input_data'].values()
    for dataset_name, metadata in group_metadata(input_metadata,
                                                 'dataset').items():
        all_vars, provenance = get_input_cubes(metadata)

        # Fix time coordinate of ERA5 instantaneous variables
        if dataset_name == 'ERA5':
            _shift_era5_time_coordinate(all_vars['tas'])

        if cfg['evaporation_method'] == 'langbein':
            logger.info("Calculation PET uisng arora method")
            all_vars.update(pet=langbein_pet(all_vars['tas']))
        else:
            logger.info("Calculation PET uisng debruin method")
            # Fix time coordinate of ERA5 instantaneous variables
            if dataset_name == 'ERA5':
                _shift_era5_time_coordinate(all_vars['psl'])
            all_vars.update(pet=debruin_pet(
                psl=all_vars['psl'],
                rsds=all_vars['rsds'],
                rsdt=all_vars['rsdt'],
                tas=all_vars['tas']))

        # Add mip to pet cube attribute
        all_vars['pet'].attributes['mip'] = all_vars['pr'].attributes['mip']
        all_vars['pet'].var_name = 'pet'

        # Change negative values for pr to zero
        _fix_negative_values(all_vars['pr'])

        for key in ['pr', 'pet']:
            cube = all_vars[key]

            # Convert unit
            _convert_units(cube)

            # Re-grid data according to the target cube
            cube = rechunk_and_regrid(cube, target_cube, cfg['regrid_scheme'])

            # Save data as an ascii file per each time step
            for sub_cube in cube.slices_over('time'):

                # Removes time dimension of length one
                new_cube = iris.util.squeeze(sub_cube)

                # Make a file name
                filename = make_filename(
                    dataset_name, cfg, new_cube, extension='asc'
                )

                # Save to ascii
                save_to_ascii(new_cube, filename)

                # Store provenance
                with ProvenanceLogger(cfg) as provenance_logger:
                    provenance_logger.log(filename, provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
