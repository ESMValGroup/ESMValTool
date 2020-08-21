"""globwat diagnostic."""
import logging
from pathlib import Path

import calendar
import os.path 
import iris
from iris.pandas import as_data_frame 
import iris.analysis    
import pandas as pd
import numpy as np
# from osgeo import gdal
# from osgeo import gdal_array
# from osgeo import osr
import matplotlib.pylab as plt
# from osgeo import gdal, osr
import os
# import shlex
# import shutil
import rasterio
from rasterio.transform import Affine


import subprocess
from esmvaltool.diag_scripts.hydrology.lazy_regrid import lazy_regrid
from esmvalcore import preprocessor as preproc
from esmvaltool.diag_scripts.hydrology.derive_evspsblpot import debruin_pet
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
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

# def get_input_cubes(metadata):
#     """Return a dict with all (preprocessed) input files."""
#     provenance = create_provenance_record()
#     all_vars = {}
#     for attributes in metadata:
#         short_name = attributes['short_name']
#         if short_name in all_vars:
#             raise ValueError(
#                 f"Multiple input files found for variable '{short_name}'.")
#         filename = attributes['filename']
#         logger.info("Loading variable %s", short_name)
#         cube = iris.load_cube(filename)
#         #set filling values to -9999
#         cube.data.set_fill_value(-9999)
#         #change cube dtype to float32
#         all_vars[short_name] = change_data_type(cube)
#         # Since the code faces memory error escaped change to floaat 32
#         # all_vars[short_name] = cube
#         provenance['ancestors'].append(filename)
#     return all_vars, provenance

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

# def regrid_temperature(src_temp, src_height, target_height, scheme):
#     """Convert temperature to target grid with lapse rate correction."""
#     # Convert 2m temperature to sea-level temperature (slt)
#     src_dtemp = lapse_rate_correction(src_height)
#     src_slt = src_temp.copy(data=src_temp.core_data() + src_dtemp.core_data())

#     # Interpolate sea-level temperature to target grid
#     target_slt = lazy_regrid(src_slt, target_height, scheme)

#     # Convert sea-level temperature to new target elevation
#     target_dtemp = lapse_rate_correction(target_height)
#     target_temp = target_slt
#     target_temp.data = target_slt.core_data() - target_dtemp.core_data()

#     return target_temp

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

# def _load_pcraster_dem(filename):
#     """Load DEM from a PCRASTER .map file."""
#     dataset = gdal.Open(str(filename))
#     lon_offset, lon_step, _, lat_offset, _, lat_step = dataset.GetGeoTransform(
#     )
#     lon_size, lat_size = dataset.RasterXSize, dataset.RasterYSize
#     data = dataset.ReadAsArray()
#     data = np.ma.masked_less(data, -1e8)
#     dataset = None

#     lons = lon_offset + lon_step * (np.arange(lon_size) + 0.5)
#     lats = lat_offset + lat_step * (np.arange(lat_size) + 0.5)

#     lon_coord = iris.coords.DimCoord(
#         lons,
#         var_name='lon',
#         standard_name='longitude',
#         units='degrees',
#     )
#     lat_coord = iris.coords.DimCoord(
#         lats,
#         var_name='lat',
#         standard_name='latitude',
#         units='degrees',
#     )

#     cube = iris.cube.Cube(
#         data,
#         var_name='height',
#         units='m',
#         dim_coords_and_dims=[
#             (lat_coord, 0),
#             (lon_coord, 1),
#         ],
#     )
#     return cube


# def check_dem(dem, cube):
#     """Check that the DEM and extract_region parameters match."""
#     for coord in ('longitude', 'latitude'):
#         start_dem_coord = dem.coord(coord).cell(0).point
#         end_dem_coord = dem.coord(coord).cell(-1).point
#         start_cube_coord = cube.coord(coord).cell(0).point
#         end_cube_coord = cube.coord(coord).cell(-1).point
#         if start_dem_coord < start_cube_coord:
#             logger.warning(
#                 "Insufficient data available, input data starts at %s "
#                 "degrees %s, but should be at least one grid "
#                 "cell larger than the DEM start at %s degrees %s.",
#                 start_cube_coord, coord, start_dem_coord, coord)
#         if end_dem_coord > end_cube_coord:
#             logger.warning(
#                 "Insufficient data available, input data ends at %s "
#                 "degrees %s, but should be at least one grid "
#                 "cell larger than the DEM end at %s degrees %s.",
#                 end_cube_coord, coord, end_dem_coord, coord)
    
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
    output_name = {'pr':{'Amon':{}, 'day':{}}, 'pet':{'Amon':{}, 'day':{}}}      
    months , days = get_month_day_for_output_name(cube) 
    for mip in 'Amon', 'day': 
        if mip == 'Amon': 
            for i in range(0, len(months)):
                for shortname in ['pr', 'pet']:
                    if shortname == 'pr':
                        monthly_pr.append('prc'+ str(months[i]) + 'wb')
                    else:
                        monthly_pet.append('eto'+ str(months[i]) + 'wb')           
        elif mip == 'day':
            for i in range(0, len(days)):
                for shortname in ['pr', 'pet']:
                    if shortname == 'pr':
                        daily_pr.append('prc'+ str(days[i]) + 'wb')
                    else:
                        daily_pet.append('eto'+ str(days[i]) + 'wb')
    output_name['pr']['Amon'] = monthly_pr
    output_name['pr']['day'] = daily_pr
    output_name['pet']['Amon']  = monthly_pet
    output_name['pet']['day'] = daily_pet
    return output_name

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
    print('input_metadata', input_metadata)
    for dataset, metadata in group_metadata(input_metadata, 'dataset').items():
        all_vars, provenance, time_step = get_input_cubes(metadata)

        target_cube_path = Path(cfg['auxiliary_data_dir']) / cfg['dem_file']
        target_cube = load_target(target_cube_path)
        # check_dem(dem, all_vars['pr'])

        # logger.info("Processing variable precipitation_flux")
        # scheme = cfg['regrid']  #'linear'
        # pr_dem = lazy_regrid(all_vars['pr'], dem, scheme)

        # logger.info("Processing variable temperature")
        # tas_dem = regrid_temperature(
        #     all_vars['tas'],
        #     all_vars['orog'],
        #     dem,
        #     scheme,
        # )

        # logger.info("Processing variable potential evapotranspiration")
        # if 'evspsblpot' in all_vars:
        #     pet = all_vars['evspsblpot']
        #     pet_dem = lazy_regrid(pet, dem, scheme)
        # else:
        #     logger.info("Potential evapotransporation not available, deriving")
        #     psl_dem = lazy_regrid(all_vars['psl'], dem, scheme)
        #     rsds_dem = lazy_regrid(all_vars['rsds'], dem, scheme)
        #     rsdt_dem = lazy_regrid(all_vars['rsdt'], dem, scheme)
        #     pet_dem = debruin_pet(
        #         tas=tas_dem,
        #         psl=psl_dem,
        #         rsds=rsds_dem,
        #         rsdt=rsdt_dem,
        #     )
        # pet_dem.var_name = 'pet'

        
        # # Interpolating variables onto the dem grid
        # # Read the target cube, which contains target grid and target elevation
        # dem_path = Path(cfg['auxiliary_data_dir']) / cfg['dem_file']
        # dem = load_dem(dem_path)
        # # # check_dem(dem, all_vars['pr'])

        # logger.info("Processing variable precipitation_flux")
        # scheme = cfg['regrid']
        # all_vars.update(pr = lazy_regrid(all_vars['pr'], dem, scheme)
        # # print('I could regrid prc')

        # # logger.info("Potential evapotransporation not available, deriving")
        # # psl_dem = lazy_regrid(all_vars['psl'], dem, scheme)
        # # rsds_dem = lazy_regrid(all_vars['rsds'], dem, scheme)
        # # rsdt_dem = lazy_regrid(all_vars['rsdt'], dem, scheme)
        # # tas_dem = lazy_regrid(all_vars['tas'], dem, scheme)
        # # all_vars.update(pet = debruin_pet(
        # #             psl= psl_dem,
        # #             rsds= rsds_dem,
        # #             rsdt= rsdt_dem,
        # #             tas= tas_dem,
        # #         ))

        # # # debruin_pet(tas=lazy_regrid(all_vars['tas'], dem, scheme),
        # # #                                   psl=lazy_regrid(all_vars['psl'], dem, scheme),
        # # #                                   rsds=lazy_regrid(all_vars['rsds'], dem, scheme),
        # # #                                   rsdt=lazy_regrid(all_vars['rsdt'], dem, scheme),
        # # #                                  )
            
        cube = all_vars['pr']  
        logger.info("Potential evapotransporation not available, deriving and adding to all_vars dictionary")
        all_vars.update(pet = debruin_pet(
            psl=all_vars['psl'],
            rsds=all_vars['rsds'],
            rsdt=all_vars['rsdt'],
            tas=all_vars['tas'],
        ))

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
            # data_dir = Path(f"{dataset}/{nyear}")
            # data_dir.mkdir(parents=True, exist_ok=True)
            for key in ['pr', 'pet']:
                for mip in ['Amon', 'day']:
                    if mip == 'Amon':
                        output_name_amon = output_name[key]['Amon']
                        data_dir = Path(f"{dataset}/{nyear}/{'Monthly'}")
                        data_dir.mkdir(parents=True, exist_ok=True)
                        for i in range(len(output_name_amon)):
                            key_cube = all_vars[key]
                            for sub_cube in key_cube.slices_over('time'): 
                                #precipitation data has negative values except than fill_value
                                #which results in gaps in model results so here we set those negative 
                                # values to zero
                                sub_cube.data[(-9999<sub_cube.data) & (sub_cube.data<0)] = 0 
                                sub_cube_regrided = sub_cube.regrid(target_cube, iris.analysis.Linear())  

                                #The preprocessed data are saved as netcdf files which then go through  
                                #a set of modifications including regriding, reprojection and netcdf to
                                #AAIGrid conversion by running a bash file in the directory where the netcdf
                                #files have been saved. 

                                #Save data as netcdf 
                                file_name = output_name_amon[i] 
                                iris.save(sub_cube_regrided, data_dir / f"{file_name}.nc", fill_value=-9999)

                                #TODO: We must work on saving data as ascii format. 
                                #We tried different saving methods including save to csv, gdal
                                # and, rasterio. However, in all methods we must use AAIGrid driver
                                #which rotate the maps 180 degree (maps are Southup insteas of Northup)
                                #TODO: The maps are Pacific-Centered while they must be Earth-centered(WGS84).
                                # We tried to apply such a change by changing latitude and longitude points, but
                                #this method was not sucessful.

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
                                #                             f"{data_dir}/{file_name}.asc", 
                                #                             'w', 
                                #                             driver='GTiff', 
                                #                             height=array.shape[1], 
                                #                             width=array.shape[0], 
                                #                             count=1, 
                                #                             dtype='float32', 
                                #                             crs='+proj=latlong', 
                                #                             transform=transform, 
                                #                             nodata= -9999,
                                #                             )    
                                # new_dataset.write(array, 1)
                                # new_dataset.close()
                    # TODO: the code faces memory error while running for daily data, need to fix 
                    elif mip == 'day':
                        output_name_day = output_name[key]['day']
                        output_name_amon = output_name[key]['Amon']
                        data_dir = Path(f"{dataset}/{nyear}/{'Daily'}")
                    for var_name in output_name: 
                        for time_step in output_name[var_name]: 
                            print(output_name[var_name][time_step]) 
                            for i in range(len(output_name_day)):
                                key_cube = all_vars[key]
                                for sub_cube in key_cube.slices_over('time'): 
                                    # frame = iris.pandas.as_data_frame(sub_cube, copy=True)
                                    sub_cube.data[(-9999<sub_cube.data) & (sub_cube.data<0)] = 0 
                                    file_name = output_name_day[i]
                                    iris.save(sub_cube, data_dir / f"{file_name}.nc", fill_value=-9999)
                                # print(frame) 
                                # final_path = os.path.join(path, output_name[key][mip][i])
                                # Path(output_file[i]).parent.mkdir(exist_ok=True)
                                # os.makedirs(final_path)
                                # frame.to_csv(output_name_day[i],
                                #             sep=' ',
                                #             na_rep='-9999',
                                #             float_format='%.1f'
                                #             ) 
            # # Store provenance
            # with ProvenanceLogger(cfg) as provenance_logger:
            #     provenance_logger.log(output_file, provenance)



if __name__ == '__main__':
   with run_diagnostic() as config:
       main(config)