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
import xarray as xr

import rasterio
from rasterio.transform import Affine

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


def _convert_units(cube, time_step):
    """Convert unit of cube."""

    cube.units = cube.units / 'kg m-3'
    cube.data = cube.core_data() / 1000.

    # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm month-1'
    if time_step == 'Amon':
        cube.convert_units('mm month-1')

    # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm day-1'
    elif time_step == 'day':
        cube.convert_units('mm day-1')
    return cube


def get_input_cubes(metadata):
    """Return a dictionary with all (preprocessed) input files."""
    provenance = create_provenance_record()
    all_vars = {}
    for attributes in metadata:
        short_name = attributes['short_name']
        time_step = attributes['mip']
        filename = attributes['filename']
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.data.set_fill_value(-9999)
        all_vars[short_name] = change_data_type(cube)
        provenance['ancestors'].append(filename)
    return all_vars, provenance, time_step


def load_target(cfg):
    """Load target grid."""
    logger.info("Reading a sample input file from %s")
    filename = Path(cfg['auxiliary_data_dir']) / cfg['target_file']
    if filename.suffix.lower() == '.nc':
        cube = iris.load_cube(str(filename))
    for coord in 'longitude', 'latitude':
        if not cube.coord(coord).has_bounds():
            logger.warning("Guessing target %s bounds", coord)
            cube.coord(coord).guess_bounds()
    return cube


def get_month_day_for_output_name(cube):
    """Get month and day number for output name."""
    coord_time = cube.coord('time')
    start_year = coord_time.cell(0).point.year
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
    monthly_pr = [] 
    daily_pr = []
    monthly_pet = [] 
    daily_pet = [] 
    monthly_pet_arora = [] 
    daily_pet_arora = [] 
    output_name = {'pr':{'Amon':{}, 'day':{}}, 
                   'pet':{'Amon':{}, 'day':{}},
                   'pet_arora':{'Amon':{}, 'day':{}}}      
    months , days = get_month_day_for_output_name(cube) 
    for mip in 'Amon', 'day': 
        if mip == 'Amon': 
            for i in range(0, len(months)):
                for shortname in ['pr', 'pet', 'pet_arora']:
                    if shortname == 'pr':
                        monthly_pr.append('prc'+ str(months[i]) + 'wb')
                    elif shortname == 'pet':
                        monthly_pet.append('eto'+ str(months[i]) + 'wb')
                    else:
                        monthly_pet_arora.append('eto_arora_'+ str(months[i]) + 'wb')                 
        elif mip == 'day':
            for i in range(0, len(days)):
                for shortname in ['pr', 'pet', 'pet_arora']:
                    if shortname == 'pr':
                        daily_pr.append('prc'+ str(days[i]) + 'wb')
                    elif shortname == 'pet':
                        daily_pet.append('eto'+ str(days[i]) + 'wb')
                    else:
                        daily_pet_arora.append('eto_arora_'+ str(days[i]) + 'wb')     
    output_name['pr']['Amon'] = monthly_pr
    output_name['pr']['day'] = daily_pr
    output_name['pet']['Amon']  = monthly_pet
    output_name['pet']['day'] = daily_pet
    output_name['pet_arora']['Amon']  = monthly_pet_arora
    output_name['pet_arora']['day']  = daily_pet_arora
    return output_name

#TODO: the function does not loop over names and it only shows prc. 
# def make_output_name(cube):
#     """Get output file name, specific to Globwat."""
#     monthly = []
#     daily = []
#     output_name = {'pr':{'Amon':{}, 'day':{}},
#                    'pet':{'Amon':{}, 'day':{}},
#                    'pet_arora':{'Amon':{}, 'day':{}}}
#     months , days = get_month_day_for_output_name(cube)

#     names = ['prc', 'eto', 'eto_arora']
#     for shortname, name in zip(output_name.keys(), names):
#         for month in months:
#             monthly.append(name + str(month) + 'wb')
#         for day in days:
#             daily.append(name + str(day) + 'wb')
#         output_name[shortname]['Amon'] = monthly
#         output_name[shortname]['day'] = daily
#     return output_name


def monthly_arora_pet(tas):
    """Calculating potential ET using Arora method.

    Arora is a temperature-based method for calculating potential ET.
    In this equation, t is the monthly average temperature in degree Celcius
    pet is in mm per month (Arora, V. K., 2002).
    https://doi.org/10.1016/S0022-1694(02)00101-4
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
    monthly_arora_pet = (a + b * tas + c * (tas ** 2))/12
    monthly_arora_pet.units = 'mm month-1'
    monthly_arora_pet.rename("arora potential evapotranspiration")
    return monthly_arora_pet


def get_cube_info(cube):
    """"""
    coord_time = cube.coord('time')
    year = coord_time.cell(0).point.year
    return year

#TODO: We must work on saving data as ascii format.
# def save_ascii(cube, filename):
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


# def affine_rotate(transform, degrees=15.0):
#     parameters = np.array(transform.GetParameters())
#     new_transform = sitk.AffineTransform(transform)
#     matrix = np.array(transform.GetMatrix()).reshape((dimension,dimension))
#     radians = -np.pi * degrees / 180.
#     rotation = np.array([[np.cos(radians), -np.sin(radians)],[np.sin(radians), np.cos(radians)]])
#     new_matrix = np.dot(rotation, matrix)
#     new_transform.SetMatrix(new_matrix.ravel())
#     resampled = resample(grid, new_transform)
#     print(new_matrix)
#     myshow(resampled, 'Rotated')
#     return new_transform
    
#TODO: the maps are saving upside down in ascii output. in addition 
# a part of Africa is cloase to USA which must be corrected
def save_to_ascii(cube, file_name):
    """Save array as ascii grid file."""
    dataset = xr.DataArray.from_iris(cube)
    array = dataset.fillna(-9999)

    xmin = array['lon'].min().values
    ymax = array['lat'].max().values
    # ymin = array['lat'].min().values

    xres = array['lon'].values[1] - array['lon'].values[0]
    yres = array['lat'].values[0] - array['lat'].values[1]

    Affine_translation = Affine.translation(xmin - xres / 2, ymax - yres / 2)
    Affine_scale = Affine.scale(xres, yres)
    transform = Affine_translation * Affine_scale
    
    new_dataset = rasterio.open(
        file_name,
        'w',
        driver='AAIGrid',
        height=array.shape[0],  # lat
        width=array.shape[1],  # lon
        count=1,
        dtype ='float32',
        # dtype=array.dtype,
        crs='+proj=latlong',
        transform=transform,
        nodata= -9999,
        )
    new_dataset.write(array, 1)
    new_dataset.close()


def _make_drs(dataset, nyear, mip):
    """making the structure of saving directory."""
    if mip == 'Amon':
        freq = 'Monthly'
    else:
        freq = 'Daily'

    data_dir = Path(f"{dataset}/{nyear}/{freq}")
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def _make_path(data_dir, file_name, extention="nc"):
    """making the saving the directory."""
    path = data_dir / f"{file_name}.{extention}"
    return path


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

        logger.info("Calculation PET uisng debruin method")
        all_vars.update(pet = debruin_pet(
            psl=all_vars['psl'],
            rsds=all_vars['rsds'],
            rsdt=all_vars['rsdt'],
            tas=all_vars['tas'],
        ))

        logger.info("Calculation PET uisng arora method")
        all_vars.update(pet_arora = monthly_arora_pet(all_vars['tas']))

        # convert unit of pr and pet
        for var in ['pr', 'pet']:
            _convert_units(all_vars[var], mip)

        output_name = make_output_name(all_vars['pr'])
        for key in ['pr', 'pet', 'pet_arora']:
            cube = all_vars[key]
            #TODO: we need a loop heir to loop over months. otherwise the code
            # stop running because output_name stuck in the first name. I just
            #added the for below which might be a better solusion for it.
            #In the previous code I looped over the lengh of the output name.
            for i in range(0,12):
                # for time in cube.coord('time'):
                # for sub_cube in cube.slices_over('time'):
                sub_cube = cube[i]
                print('sub_cub',sub_cube)
                nyear = get_cube_info(sub_cube)
                print ('sub_cube*********', sub_cube)

                # set negative values for precipitaion to zero
                if key == 'pr':
                    sub_cube.data[(-9999<sub_cube.data) & (sub_cube.data<0)] = 0
                #TODO: here the code faces memory error:
                #"MemoryError: Unable to allocate 285. MiB for an array 
                # with shape (37324800,) and data type int64
                
                # # load target grid for using in regriding function"
                # target_cube = load_target(cfg)
                # #regrid data baed on target cube
                # sub_cube_regrided = sub_cube.regrid(target_cube, iris.analysis.Linear())

                output_name_save = output_name[key][mip]
                file_name = output_name_save[i]
                data_dir = _make_drs(dataset, nyear, mip)

                #Save data as netcdf
                path = _make_path(data_dir, file_name, extention="nc")
                iris.save(sub_cube, path , fill_value=-9999)

                # save_ascii(sub_cube_regrided, path)
                path = _make_path(data_dir, file_name, extention="asc")
                save_to_ascii(sub_cube, path)

            # # Store provenance
            # with ProvenanceLogger(cfg) as provenance_logger:
            #     provenance_logger.log(path, provenance)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
