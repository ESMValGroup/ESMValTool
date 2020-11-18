"""Globwat diagnostic.

Currently, the model reads monthly reference evaporation and then multiplies
it by kc and then divide it by the number of days in each month and uses it
as daily actual evaporation in the model. We can change the model code to
read daily reference evaporation and does the calculation, but we
did not change the model code for this aim yet.
"""
import logging
from pathlib import Path

import iris
import iris.analysis
import numpy as np
import xarray as xr
import pandas as pd
import dask.array as da

from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.hydrology.lazy_regrid import lazy_regrid
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
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
        ],
        'ancestors': [],
    }
    return record


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
    """Convert unit of cube."""
    cube.units = cube.units / 'kg m-3'
    cube.data = cube.core_data() / 1000.

    mip = cube.attributes['mip']

    # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm month-1'
    if mip == 'Amon':
        cube.convert_units('mm month-1')

    # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm day-1'
    elif mip == 'day':
        cube.convert_units('mm day-1')
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
        provenance['ancestors'].append(filename)
    return all_vars, provenance


def load_target(cfg):
    """Load target grid."""
    filename = Path(cfg['auxiliary_data_dir']) / cfg['target_file']
    cube = iris.load_cube(str(filename))
    return cube


def arora_pet(tas):
    """Calculate potential ET using Arora method.

    Arora is a temperature-based method for calculating potential ET.
    In this equation, t is the monthly average temperature in degree Celcius
    pet is in mm per month (Arora, V. K., 2002).
    https://doi.org/10.1016/S0022-1694(02)00101-4
    """
    tas.convert_units('degC')
    constant_a = iris.coords.AuxCoord(np.float32(325),
                                      long_name='first constant', units=None)
    constant_b = iris.coords.AuxCoord(np.float32(21),
                                      long_name='second constant', units=None)
    constant_c = iris.coords.AuxCoord(np.float32(0.9),
                                      long_name='third constant', units=None)
    mip = tas.attributes['mip'] 
    if mip == "Amon":
        conversion = 12
        unit = 'mm month-1'
    else:
        conversion = tas.coord('time').shape[0]
        unit = 'mm day-1'
    # assumption here: tas is constant over time, then the monthly/daily    
    # average value is equal to the annual average.
    pet_annual = constant_a + constant_b * (tas) + constant_c * (tas ** 2)
    pet = pet_annual / conversion
    pet.units = unit
    pet.var_name = 'evspsblpot'
    pet.standard_name = 'water_potential_evaporation_flux'
    pet.long_name = 'Potential Evapotranspiration'
    pet.attributes['mip'] = tas.attributes['mip']
    return pet


def get_cube_time_info(cube):
    """Get year, month and day from the cube."""
    coord_time = cube.coord('time')
    year = coord_time.cell(0).point.year
    month = str(coord_time.cell(0).point.month).zfill(2)
    day = str(coord_time.cell(0).point.day).zfill(2)
    return year, month, day


def get_cube_data_info(cube):
    """Get short_name, and mip from the cube."""
    short_name = cube.var_name
    mip = cube.attributes['mip']
    return short_name, mip


def _reindex_data(cube, target):
    """Reindex data to a global coordinates set.

    Globwat works with global extent.
    """
    array = xr.DataArray.from_iris(cube).to_dataset()
    target_ds = xr.DataArray.from_iris(target).to_dataset()
    reindex_ds = array.reindex(
        {"lat": target_ds["lat"], "lon": target_ds["lon"]},
        method="nearest",
        tolerance=1e-2,
    )
    new_cube = reindex_ds.to_array().to_iris()
    new_cube.var_name = cube.var_name
    new_cube.attributes = cube.attributes
    return new_cube


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


def _flip_vertically(array):
    """Flip vertically for writing as ascii.

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
    array = _flip_vertically(array)

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
    """Get a valid path for saving a diagnostic data file.

    filenames are specific to Globwat.
    """
    names_map = {'pr': 'prc', 'evspsblpot': 'eto'}

    nyear, nmonth, nday = get_cube_time_info(cube)
    short_name, mip = get_cube_data_info(cube)

    if mip == 'Amon':
        filename = f"{names_map[short_name]}{nmonth}wb.{extension}"
        freq = 'Monthly'

    else:
        filename = f"{names_map[short_name]}{nmonth}{nday}wb.{extension}"
        freq = 'Daily'

    data_dir = Path(f"{cfg['work_dir']}/{dataset_name}/{nyear}/{freq}")
    data_dir.mkdir(parents=True, exist_ok=True)
    return str(data_dir / f"{filename}")


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
        if cfg['pet_arora']:
            logger.info("Calculation PET uisng arora method")
            all_vars.update(pet=arora_pet(all_vars['tas']))
        else:
            logger.info("Calculation PET uisng debruin method")
            all_vars.update(pet=debruin_pet(
                psl=all_vars['psl'],
                rsds=all_vars['rsds'],
                rsdt=all_vars['rsdt'],
                tas=all_vars['tas']))
            # Convert unit of pet
            _convert_units(all_vars['pet'])

        # Convert unit of pr
        _convert_units(all_vars['pr'])

        # Change negative values for pr to zero
        _fix_negative_values(all_vars['pr'])

        for key in ['pr', 'pet']:
            cube = all_vars[key]

            # Re-grid data according to the target cube
            cube = lazy_regrid(cube, target_cube, cfg['regrid_scheme'])

            # Reindex data to a global coordinates set
            cube = _reindex_data(cube, target_cube)

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
