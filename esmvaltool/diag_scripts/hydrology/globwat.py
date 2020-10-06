"""Globwat diagnostic.

Currently, the model reads monthly reference evaporation and then multiplies it by kc and then
divide it by the number of days in each month and uses it as daily actual evaporation in the model.
We can change the model code to read daily reference evaporation and does the calculation, but we
did not change the model code for this aim yet.
"""
import logging
from pathlib import Path

import calendar
import iris
import iris.analysis
import numpy as np

from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.shared import (get_diagnostic_filename,
                                            group_metadata, run_diagnostic,
                                            select_metadata)


logger = logging.getLogger(Path(__file__).name) #for test in python use '__file__'


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
        #set filling values to -9999
        cube.data.set_fill_value(-9999)
        #change cube dtype to float32
        all_vars[short_name] = change_data_type(cube)
        # Since the code faces memory error escaped change to floaat 32
        # all_vars[short_name] = cube
        provenance['ancestors'].append(filename)
    return all_vars, provenance, time_step

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

def load_target(filename):
    """Load DEM into iris cube."""
    logger.info("Reading digital elevation model from %s", filename)
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
    """get month and day number for output name."""
    # get start year and end year
    coord_time = cube.coord('time')
    start_year = coord_time.cell(0).point.year
    end_year = coord_time.cell(-1).point.year
    year = start_year - 1
    #looping over time dimention, get day and months
    months = []
    days = []
    for i in range(0, len(cube.coord('time').points)):
        n_month = str(coord_time.cell(i).point.month).zfill(2)
        months.append(n_month)
        nday = calendar.monthrange(year,int(n_month))[1]
        for daynumber in range(1, nday+1):
            days.append(str(n_month) + str(daynumber).zfill(2))
    return months ,days

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
                        monthly_pet_arora.append('eto_arora'+ str(months[i]) + 'wb')
        elif mip == 'day':
            for i in range(0, len(days)):
                for shortname in ['pr', 'pet', 'pet_arora']:
                    if shortname == 'pr':
                        daily_pr.append('prc'+ str(days[i]) + 'wb')
                    elif shortname == 'pet':
                        daily_pet.append('eto'+ str(days[i]) + 'wb')
                    else:
                        daily_pet_arora.append('eto_arora'+ str(days[i]) + 'wb')
    output_name['pr']['Amon'] = monthly_pr
    output_name['pr']['day'] = daily_pr
    output_name['pet']['Amon']  = monthly_pet
    output_name['pet']['day'] = daily_pet
    output_name['pet_arora']['Amon']  = monthly_pet_arora
    output_name['pet_arora']['day']  = daily_pet_arora
    return output_name

def arora_pet(tas):
    """Arora is a temperature-based method for calculating potential ET.
    In this equation, t is the monthly average temperature in degree Celcius
    pet is in mm per month."""
    # convert temperature unit to degC
    tas.convert_units('degC')
    # define constants of formula
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

def main(cfg):
    """Process data for use as input to the GlobWat hydrological model.
    These variables are needed in all_vars:
    pr (precipitation_flux)
    psl (air_pressure_at_mean_sea_level)
    rsds (surface_downwelling_shortwave_flux_in_air)
    rsdt (toa_incoming_shortwave_flux)
    tas (air_temperature)
    """
    input_metadata = cfg['input_data'].values()
    # for checking the code in ipython I added print(input_metadata).
    # Run the script and use print(input_metadata) in the log.txt as input_metadata
    # print('input_metadata', input_metadata)
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance, time_step = get_input_cubes(metadata)
        #load target grid for using in regriding function
        target_cube_path = Path(cfg['auxiliary_data_dir']) / cfg['target_file']
        target_cube = load_target(target_cube_path)

        cube = all_vars['pr']
        logger.info("Potential evapotransporation not available, deriving and adding to all_vars dictionary")
        all_vars.update(pet = debruin_pet(
            psl=all_vars['psl'],
            rsds=all_vars['rsds'],
            rsdt=all_vars['rsdt'],
            tas=all_vars['tas'],
        ))
        #Take tas values to calculate pet_arora
        tas = all_vars['tas']
        all_vars.update(pet_arora = arora_pet(tas))

        logger.info("Converting units")
        pr = all_vars['pr']
        pet = all_vars['pet']
        pr.units = pr.units / 'kg m-3'
        pr.data = pr.core_data() / 1000.
        pet.units = pet.units / 'kg m-3'
        pet.data = pet.core_data() /1000.
        # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm month-1'
        if time_step['mip'] == 'Amon':
            pr.convert_units('mm month-1')
            pet.convert_units('mm month-1')
        # Unit conversion of pr and pet from 'kg m-2 s-1' to 'mm day-1'
        elif time_step['mip'] == 'day':
            pr.convert_units('mm day-1')
            pet.convert_units('mm day-1')

        coord_time = cube.coord('time')
        output_name = make_output_name(cube)
        start_year = coord_time.cell(0).point.year
        end_year = coord_time.cell(-1).point.year
        for nyear in range (start_year , end_year+1):
            for key in ['pr', 'pet', 'pet_arora']:
                if time_step['mip'] == 'Amon':
                    output_name_amon = output_name[key]['Amon']
                    data_dir = Path(f"{dataset}/{nyear}/{'Monthly'}")
                    data_dir.mkdir(parents=True, exist_ok=True)
                    for i in range(len(output_name_amon)):
                        key_cube = all_vars[key]
                        for sub_cube in key_cube.slices_over('time'):
                            # set negative values for precipitaion to zero
                            if key == 'pr':
                                sub_cube.data[(-9999<sub_cube.data) & (sub_cube.data<0)] = 0
                            #regrid data baed on target cube
                            sub_cube_regrided = sub_cube.regrid(target_cube, iris.analysis.Linear())
                            #Save data as netcdf
                            file_name = output_name_amon[i]
                            iris.save(sub_cube_regrided, data_dir / f"{file_name}.nc", fill_value=-9999)

                            #TODO: We must work on saving data as ascii format.
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
                else:
                    #getting output names
                    output_name_amon = output_name[key]['Amon']
                    output_name_day = output_name[key]['day']
                    #making output directory
                    data_dir = Path(f"{dataset}/{nyear}/{'Daily'}")
                    data_dir.mkdir(parents=True, exist_ok=True)
                    #saving output
                    for i in range(len(output_name_day)):
                        key_cube = all_vars[key]
                        for sub_cube in key_cube.slices_over('time'):
                            # set negative values for precipitaion to zero
                            if key == 'pr':
                                sub_cube.data[(-9999<sub_cube.data) & (sub_cube.data<0)] = 0

                            sub_cube_regrided = sub_cube.regrid(target_cube, iris.analysis.Linear())
                            file_name = output_name_day[i]
                            iris.save(sub_cube, data_dir / f"{file_name}.nc", fill_value=-9999)

            # # Store provenance
            # with ProvenanceLogger(cfg) as provenance_logger:
            #     provenance_logger.log(output_file, provenance)



if __name__ == '__main__':
   with run_diagnostic() as config:
       main(config)
