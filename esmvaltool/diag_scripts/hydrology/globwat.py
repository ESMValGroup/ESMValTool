"""Globwat diagnostic.

Currently, the model reads monthly reference evaporation and then multiplies
it by kc and then divide it by the number of days in each month and uses it
as daily actual evaporation in the model. We can change the model code to
read daily reference evaporation and does the calculation, but we
did not change the model code for this aim yet.
"""
import logging
from pathlib import Path

import calendar
import iris
import iris.analysis
import numpy as np

from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.shared import (ProvenanceLogger, group_metadata, run_diagnostic)


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


def get_input_cubes(metadata):
    """Return a dictionary with all (preprocessed) input files."""
    provenance = create_provenance_record()
    all_vars = {}
    time_step = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        time_step['mip'] = attributes['mip']
        for key,value in time_step.items():
            if value not in time_step.values():
                time_step[key] = value
        if short_name in all_vars:
            raise ValueError(
                f"Multiple input files found for variable '{short_name}'.")
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.data.set_fill_value(-9999)
        all_vars[short_name] = change_data_type(cube)

        # convert unit of pr and pet
        for var in ['pr', 'pet']:
            _convert_units(all_vars[var], time_step)
        provenance['ancestors'].append(filename)
    return all_vars, provenance, time_step


def load_target(cfg):
    """Load target grid."""
    logger.info("Reading digital elevation model from %s", filename)
    filename = Path(cfg['auxiliary_data_dir']) / cfg['target_file']
    if filename.suffix.lower() == '.nc':
        cube = iris.load_cube(str(filename))
    elif filename.suffix.lower() == '.map':
        cube = _load_pcraster_dem(filename)
    else:
        raise ValueError(f"Unknown file format {filename}. Supported formats "
                         "are '.nc' and '.map'.")
    for coord in 'longitude', 'latitude':
        if not cube.coord(coord).has_bounds():
            logger.warning("Guessing DEM %s bounds", coord)
            cube.coord(coord).guess_bounds()
    return cube


def get_month_day_for_output_name(cube):
    """Get month and day number for output name."""
    coord_time = cube.coord('time')
    start_year = coord_time.cell(0).point.year
    end_year = coord_time.cell(-1).point.year
    year = start_year - 1
    months = []
    days = []
    for i in range(0, len(cube.coord('time').points)):
        n_month = str(coord_time.cell(i).point.month).zfill(2)
        months.append(n_month)
        nday = calendar.monthrange(year,int(n_month))[1]
        for daynumber in range(1, nday+1):
            days.append(str(n_month) + str(daynumber).zfill(2))
    return months , days


def make_output_name(cube):
    """Get output file name, specific to Globwat."""
    monthly = []
    daily = []
    output_name = {'pr':{'Amon':{}, 'day':{}},
                   'pet':{'Amon':{}, 'day':{}},
                   'pet_arora':{'Amon':{}, 'day':{}}}
    months , days = get_month_day_for_output_name(cube)

    names = ['prc', 'eto', 'eto_arora']
    for shortname, name in zip(output_name.keys(), names):
        for month in months:
            monthly.append(name + str(month) + 'wb')
        for day in days:
            daily.append(name + str(day) + 'wb')
        output_name[shortname]['Amon'] = monthly
        output_name[shortname]['day'] = daily
    return output_name


def arora_pet(tas):
    """Calculating potential ET using Arora method.

    Arora is a temperature-based method for calculating potential ET.
    In this equation, t is the monthly average temperature in degree Celcius
    pet is in mm per month.
    """
    tas.convert_units('degC')
    a = iris.coords.AuxCoord(np.float32(325),
                             long_name='first constant',
                             units=None)
    b = iris.coords.AuxCoord(np.float32(21),
                             long_name='second constant',
                             units=None)
    c = iris.coords.AuxCoord(np.float32(0.9),
                             long_name='third constant',
                             units=None)
    arora_pet = a + b * tas + c * (tas ** 2)
    arora_pet.units = 'mm month-1'
    arora_pet.rename("arora potential evapotranspiration")
    return arora_pet


def _convert_units(cube, time_step):
    """Convert unit of cube."""

    cube.units = cube.units / 'kg m-3'
    cube.data = cube.core_data() / 1000.

    # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm month-1'
    if time_step['mip'] == 'Amon':
        cube.convert_units('mm month-1')

    # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm day-1'
    elif time_step['mip'] == 'day':
        cube.convert_units('mm day-1')
    return cube

def get_cube_info(cube):
    """"""
    coord_time = cube.coord('time')
    start_year = coord_time.cell(0).point.year
    end_year = coord_time.cell(-1).point.year
    return start_year, end_year

#TODO: We must work on saving data as ascii format.
def save_ascii(cube, filename):
    """"""
    # array = sub_cube.data
    # points = sub_cube.coord('longitude').points
    # sub_cube.coord('longitude').points = (points + 180) % 360 - 180
    # lon = (sub_cube.coord('longitude').points + 180) % 360 -180
    # lat = sub_cube.coord('latitude').points
    # sub_cube = sub_cube[:, ::-1, ...]
    # nrows,ncols = np.shape(array)
    # xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
    # xres = (xmax-xmin)/float(ncols)
    # yres = (ymax-ymin)/float(nrows)
    # transform = Affine.translation(xmin + xres / 2, ymin - xres / 2) * Affine.scale(xres, xres)
    # file_name = output_name_amon[i]
    # new_dataset = rasterio.open(
    # f"{data_dir}/{file_name}.asc",
    # 'w',
    # driver='GTiff',
    # height=array.shape[1],
    # width=array.shape[0],
    # count=1,
    # dtype='float32',
    # crs='+proj=latlong',
    # transform=transform,
    # nodata= -9999,
    # )
    # new_dataset.write(array, 1)
    # new_dataset.close()


def main(cfg):
    """Process data for GlobWat hydrological model.

    These variables are needed in all_vars:
    pr (precipitation_flux)
    psl (air_pressure_at_mean_sea_level)
    rsds (surface_downwelling_shortwave_flux_in_air)
    rsdt (toa_incoming_shortwave_flux)
    tas (air_temperature)
    """
    input_metadata = cfg['input_data'].values()
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance, mip = get_input_cubes(metadata)

        # load target grid for using in regriding function
        target_cube = load_target(cfg)


        logger.info("Calculation PET uisng debruin method")
        all_vars.update(pet = debruin_pet(
            psl=all_vars['psl'],
            rsds=all_vars['rsds'],
            rsdt=all_vars['rsdt'],
            tas=all_vars['tas'],
        ))

        logger.info("Calculation PET uisng arora method")
        all_vars.update(pet_arora = arora_pet(all_vars['tas']))

        output_name = make_output_name(all_vars['pr'])
        start_year, end_year = get_cube_info(all_vars['pr'])

        if mip == 'Amon':
            freq = 'Monthly'
        else:
            freq = 'Daily'

        for nyear in range(start_year , end_year+1):
            for key in ['pr', 'pet', 'pet_arora']:
                file_names = output_name[key][mip]
                data_dir = Path(f"{dataset}/{nyear}/{freq}")
                data_dir.mkdir(parents=True, exist_ok=True)
                for file_name in file_names:
                    cube = all_vars[key]
                    for sub_cube in cube.slices_over('time'):
                        # set negative values for precipitaion to zero
                        if key == 'pr':
                            sub_cube.data[(-9999<sub_cube.data) & (sub_cube.data<0)] = 0
                        #regrid data baed on target cube
                        sub_cube_regrided = sub_cube.regrid(target_cube, iris.analysis.Linear())
                        #Save data as netcdf
                        path = data_dir / f"{file_name}.nc"
                        iris.save(sub_cube_regrided, path , fill_value=-9999)

                        # save_ascii(sub_cube, path)

            # Store provenance
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(path, provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
