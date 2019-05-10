# pylint: disable=invalid-name
"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
   http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid

Last access
   20190510

Download and processing instructions
   Download the desired files from:
   https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/

"""

import logging
import os
import glob

import numpy as np
from cf_units import Unit
import iris
from iris.coords import AuxCoord, DimCoord

from .utilities import (constant_metadata, fix_coords,
                        fix_var_metadata, read_cmor_config, save_variable,
                        set_global_atts)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = read_cmor_config('PIOMAS.yml')


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in CFG['variables'].items():
        grid_path = os.path.join(in_dir, vals['grid_file'])
        file_expression = os.path.join(in_dir, '{0}.H????'.format(vals['raw']))
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        for file_path in glob.glob(file_expression):
            read_binary_file(file_path, grid_path, var_info, glob_attrs, out_dir)


def read_binary_file(data_path, grid_path, var_info, attrs, out_dir):
    nx = 360
    ny = 120

    grids = np.loadtxt(grid_path)
    grids = grids.reshape(2, ny, nx)  #2010

    lon_coord = AuxCoord(
        grids[0, ...],
        standard_name='longitude',
        var_name='lon',
        units='degrees_east'
    )

    lat_coord = AuxCoord(
        grids[1, ...],
        standard_name='latitude',
        var_name='lat',
        units='degrees_north'
    )

    fd_data = open(data_path, 'rb')
    data = np.fromfile(fd_data, dtype=np.dtype('f'), count=-1)
    days = data.shape[0] // nx // ny
    year = int(data_path[-4:])

    time_coord = DimCoord(
        np.arange(0, days),
        standard_name='time',
        var_name='time',
        units=Unit('days since {}-01-01'.format(year), calendar='noleap'),
    )
    data = data.reshape(days, ny, nx)

    cube = iris.cube.Cube(
        data,
        standard_name=var_info.standard_name,
        var_name=var_info.short_name,
        units='m',
    )
    cube.add_dim_coord(time_coord, 0)
    cube.add_aux_coord(lon_coord, (1, 2))
    cube.add_aux_coord(lat_coord, (1, 2))

    logger.info(cube)

    # fix_var_metadata(cube, var_info)
    # fix_coords(cube)
    # _fix_data(cube, var)
    set_global_atts(cube, attrs)
    save_variable(cube, var_info.short_name, out_dir, attrs)

    # for month in range(1, 13):
    #     sicvar.units = '1'
    #     sicvar.missing_value = '-9.e+33'
    #     timevar_area.units = 'days since '+str(y)+'-'+str(month).zfill(2)+'-01'
    #     areavar_area.units = 'm2'
    #     lonvar_area[:,:] = grids[1,:,:]
    #     latvar_area[:,:] = grids[0,:,:]
    #     #  areavar_area[:,:] = 0.5*np.abs(np.sin(np.deg2rad(90-grids[6,:,:])))*(grids[2,:,:]*grids[3,:,:]+grids[4,:,:]*grids[5,:,:])
    #     areavar_area[:,:] = 0.5*(grids[2,:,:]*grids[3,:,:]+grids[4,:,:]*grids[5,:,:])
    #     timevar_area[:] = [15]
    #     sicvar[0,:,:] = area[month-1,:,:]
    #     # sicvar[sicvar == 0.] = 9999
    #     # sicvar = np.where(sicvar==0.,9999.,sicvar)
    #     ncfile_area.close
    #     #heff
    #     ncfilename_heff = basedir+'monthly_mean/sithick/sithick_'+str(y)+str(month).zfill(2)+'.nc'
    #     ncfile_heff = nc.Dataset(ncfilename_heff,"w",format='NETCDF4')
    #     ncfile_heff.createDimension('lat',ny)
    #     ncfile_heff.createDimension('lon',nx)
    #     ncfile_heff.createDimension('time',None)
    #     lonvar_heff = ncfile_heff.createVariable('lon','f',('lat','lon'),zlib=True)
    #     latvar_heff = ncfile_heff.createVariable('lat','f',('lat','lon'),zlib=True)
    #     timevar_heff = ncfile_heff.createVariable('time','i',('time'),zlib=True)
    #     areavar_heff = ncfile_heff.createVariable('area','f',('lat','lon'),zlib=True)
    #     sithickvar = ncfile_heff.createVariable('sithick','f',('time','lat','lon'),zlib=True)
    #     sithickvar.units = 'm'
    #     timevar_heff.units = 'days since '+str(y)+'-'+str(month).zfill(2)+'-01'
    #     areavar_heff.units = 'm2'
    #     lonvar_heff[:,:] = grids[1,:,:]
    #     latvar_heff[:,:] = grids[0,:,:]
    #     # areavar_heff[:,:] = 0.5*np.abs(np.sin(np.deg2rad(90-grids[6,:,:])))*(grids[2,:,:]*grids[3,:,:]+grids[4,:,:]*grids[5,:,:])
    #     areavar_heff[:,:] = 0.5*(grids[2,:,:]*grids[3,:,:]+grids[4,:,:]*grids[5,:,:])
    #     timevar_heff[:] = [15]
    #     sithickvar[0,:,:] = heff[month-1,:,:]
    #     # sithickvar[sithickvar == 0.] = 9999
    #     # sithickvar = np.where(sithickvar==0,-9999.,sithickvar)
    #     ncfile_heff.close()