# -*- coding: utf-8 -*-
"""Part of the ESMValTool Arctic Ocean diagnostics.

This module contains functions for extracting the data
from netCDF files and prepearing them for plotting.
"""
import logging
import os
import ESMF
import numpy as np
from netCDF4 import Dataset, num2date

from esmvaltool.diag_scripts.arctic_ocean.regions import (hofm_regions,
                                                          transect_points)
from esmvaltool.diag_scripts.arctic_ocean.utils import (genfilename,
                                                        point_distance,
                                                        get_fx_filenames,
                                                        get_series_lenght,
                                                        get_provenance_record)
from esmvaltool.diag_scripts.shared import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


def load_meta(datapath, fxpath=None):
    """Load metadata of the netCDF file.

    Parameters
    ----------
    datapath: str
        path to the netCDF file with data
    fxpath: str
        path to the netCDF file with fx files

    Returns
    -------
    datafile: instance of netCDF4 Dataset
        points to the file
    lon2d: numpy array
        Two dimentional longitude information
    lat2d: numpy array
        Two dimentional latitude information
    lev: numpy array
        depths of the model levels
    time: numpy array
        dates are converted to datetime objects
    areacello: numpy array
        values of areacello
    """
    datafile = Dataset(datapath)

    if fxpath:
        datafile_area = Dataset(fxpath)
        areacello = datafile_area.variables['areacello'][:]
    else:
        areacello = None

    lon = datafile.variables['lon'][:]
    lat = datafile.variables['lat'][:]
    lev = datafile.variables['lev'][:]
    time = num2date(datafile.variables['time'][:],
                    datafile.variables['time'].units)
    # hack for HadGEM2-ES
    lat[lat > 90] = 90

    if lon.ndim == 2:
        lon2d, lat2d = lon, lat
    elif lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)

    metadata = {}
    metadata['datafile'] = datafile
    metadata['lon2d'] = lon2d
    metadata['lat2d'] = lat2d
    metadata['lev'] = lev
    metadata['time'] = time
    metadata['areacello'] = areacello
    return metadata


def hofm_extract_region(metadata, cmor_var, indexes, level, time=0):
    """Calculates mean over the region."""
    # fix for climatology
    if metadata['datafile'].variables[cmor_var].ndim < 4:
        level_pp = metadata['datafile'].variables[cmor_var][level, :, :]
    else:
        level_pp = metadata['datafile'].variables[cmor_var][time, level, :, :]
    if not isinstance(level_pp, np.ma.MaskedArray):
        level_pp = np.ma.masked_equal(level_pp, 0)
    data_mask = level_pp[indexes[0], indexes[1]].mask
    area_masked = np.ma.masked_where(
        data_mask, metadata['areacello'][indexes[0], indexes[1]])
    result = (area_masked *
              level_pp[indexes[0], indexes[1]]).sum() / area_masked.sum()
    return result


def hofm_save_data(cfg, data_info, oce_hofm):
    """Save data for Hovmoeller diagrams."""

    ofiles = {}
    ofiles['ofilename'] = genfilename(**data_info, data_type='hofm')
    ofiles['ofilename_levels'] = genfilename(**data_info, data_type='levels')
    ofiles['ofilename_time'] = genfilename(**data_info, data_type='time')

    np.save(ofiles['ofilename'], oce_hofm)
    provenance_record = get_provenance_record(data_info, 'hofm', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename'] + '.npy', provenance_record)

    if isinstance(data_info['levels'], np.ma.core.MaskedArray):
        np.save(ofiles['ofilename_levels'],
                data_info['levels'][0:data_info['lev_limit']].filled())
    else:
        np.save(ofiles['ofilename_levels'],
                data_info['levels'][0:data_info['lev_limit']])
    provenance_record = get_provenance_record(data_info, 'lev', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_levels'] + '.npy',
                              provenance_record)

    np.save(ofiles['ofilename_time'], data_info['time'])
    provenance_record = get_provenance_record(data_info, 'time', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_time'] + '.npy',
                              provenance_record)


def hofm_data(cfg, model_filenames, mmodel, cmor_var, region):
    """Extract data for Hovmoeller diagrams from monthly values.

    Saves the data to files in `diagworkdir`.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    mmodel: str
        model name that will be processed.
    cmor_var: str
        name of the CMOR variable
    areacello_fx: OrderedDict.
        dictionary with model names as keys and paths to fx files as values.
    max_level: float
        maximum depth level the Hovmoeller diagrams should go to.
    region: str
        name of the region predefined in `hofm_regions` function.
    diagworkdir: str
        path to work directory.

    Returns
    -------
    None
    """
    logger.info("Extract  %s data for %s, region %s", cmor_var, mmodel, region)
    areacello_fx = get_fx_filenames(cfg, 'areacello')
    metadata = load_meta(datapath=model_filenames[mmodel],
                         fxpath=areacello_fx[mmodel])

    lev_limit = metadata['lev'][
        metadata['lev'] <= cfg['hofm_depth']].shape[0] + 1

    indexes = hofm_regions(region, metadata['lon2d'], metadata['lat2d'])

    series_lenght = get_series_lenght(metadata['datafile'], cmor_var)

    oce_hofm = np.zeros((metadata['lev'][0:lev_limit].shape[0], series_lenght))
    for mon in range(series_lenght):
        for ind, _ in enumerate(metadata['lev'][0:lev_limit]):
            oce_hofm[ind, mon] = hofm_extract_region(metadata, cmor_var,
                                                     indexes, ind, mon)
    data_info = {}
    data_info['basedir'] = cfg['work_dir']
    data_info['variable'] = cmor_var
    data_info['mmodel'] = mmodel
    data_info['region'] = region
    data_info['time'] = metadata['time']
    data_info['levels'] = metadata['lev']
    data_info['lev_limit'] = lev_limit
    data_info['ori_file'] = model_filenames[mmodel]
    data_info['areacello'] = areacello_fx[mmodel]

    hofm_save_data(cfg, data_info, oce_hofm)

    metadata['datafile'].close()


def transect_level(datafile, cmor_var, level, grid, locstream):
    """Interpolation for one level of transect."""

    sourcefield = ESMF.Field(
        grid,
        staggerloc=ESMF.StaggerLoc.CENTER,
        name='MPI',
    )
    # load model data
    model_data = datafile.variables[cmor_var][0, level, :, :]

    # ESMF do not understand masked arrays, so fill them
    if isinstance(model_data, np.ma.core.MaskedArray):
        sourcefield.data[...] = model_data.filled(0).T
    else:
        sourcefield.data[...] = model_data.T

    # create a field we giong to intorpolate TO
    dstfield = ESMF.Field(locstream, name='dstfield')
    dstfield.data[:] = 0.0

    # create an object to regrid data
    # from the source to the destination field
    dst_mask_values = None
    # if domask:
    dst_mask_values = np.array([0])

    regrid = ESMF.Regrid(
        sourcefield,
        dstfield,
        regrid_method=ESMF.RegridMethod.NEAREST_STOD,
        # regrid_method=ESMF.RegridMethod.BILINEAR,
        unmapped_action=ESMF.UnmappedAction.IGNORE,
        dst_mask_values=dst_mask_values)

    # do the regridding from source to destination field
    dstfield = regrid(sourcefield, dstfield)
    return dstfield


def transect_save_data(cfg, data_info, secfield, lon_s4new, lat_s4new):
    """Save data for transects."""

    ofiles = {}
    ofiles['ofilename'] = genfilename(**data_info, data_type='transect')
    ofiles['ofilename_depth'] = genfilename(
        data_info['basedir'], 'depth', data_info['mmodel'],
        data_info['region'], 'transect_' + data_info['variable'])
    ofiles['ofilename_dist'] = genfilename(data_info['basedir'], 'distance',
                                           data_info['mmodel'],
                                           data_info['region'],
                                           'transect_' + data_info['variable'])

    np.save(ofiles['ofilename'], secfield)
    print(ofiles['ofilename'])
    provenance_record = get_provenance_record(data_info, 'transect', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename'] + '.npy', provenance_record)
    # we have to fill masked arrays before saving

    if isinstance(data_info['levels'], np.ma.core.MaskedArray):
        np.save(ofiles['ofilename_depth'], data_info['levels'].filled())
    else:
        np.save(ofiles['ofilename_depth'], data_info['levels'])
    print(ofiles['ofilename_depth'])
    provenance_record = get_provenance_record(data_info, 'levels', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_depth'] + '.npy',
                              provenance_record)

    np.save(ofiles['ofilename_dist'], point_distance(lon_s4new, lat_s4new))
    provenance_record = get_provenance_record(data_info, 'distance', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_dist'] + '.npy',
                              provenance_record)
    print(ofiles['ofilename_dist'])


def transect_data(cfg, mmodel, cmor_var, region, mult=2):
    """Extract data for transects (defined in regions.transect_points).

    Parameters
    ----------
    mmodel: str
        model name that will be processed.
    cmor_var: str
        name of the CMOR variable
    region: str
        name of the region predefined in `transect_points` function.
    diagworkdir: str
        path to work directory.
    mult: integer
        multiplicator for the number of points in the transect.
        Can be used to increase transect resolution.
    observations: str
        name of the observation dataset.
    """
    logger.info("Extract  %s transect data for %s, region %s", cmor_var,
                mmodel, region)
    # get the path to preprocessed file
    ifilename = genfilename(cfg['work_dir'],
                            cmor_var,
                            mmodel,
                            data_type='timmean',
                            extension='.nc')
    # open with netCDF4
    datafile = Dataset(ifilename)
    # open with ESMF
    grid = ESMF.Grid(filename=ifilename, filetype=ESMF.FileFormat.GRIDSPEC)

    # get depth of the levels
    lev = datafile.variables['lev'][:]

    # indexesi, indexesj = hofm_regions(region, lon2d, lat2d)
    lon_s4new, lat_s4new = transect_points(region, mult=mult)

    # masking true
    # domask = True

    # create instans of the location stream (set of points)
    locstream = ESMF.LocStream(lon_s4new.shape[0],
                               name="Atlantic Inflow Section",
                               coord_sys=ESMF.CoordSys.SPH_DEG)

    # appoint the section locations
    locstream["ESMF:Lon"] = lon_s4new
    locstream["ESMF:Lat"] = lat_s4new
    # if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_s4new.shape[0]),
                                      dtype=np.int32)
    # initialise array for the section
    secfield = np.zeros(
        (lon_s4new.shape[0], datafile.variables[cmor_var].shape[1]))

    # loop over depth levels
    for level in range(0, datafile.variables[cmor_var].shape[1]):
        secfield[:, level] = transect_level(datafile, cmor_var, level, grid,
                                            locstream).data
    data_info = {}
    data_info['basedir'] = cfg['work_dir']
    data_info['variable'] = cmor_var
    data_info['mmodel'] = mmodel
    data_info['region'] = region
    data_info['levels'] = lev
    data_info['ori_file'] = ifilename
    data_info['areacello'] = None

    transect_save_data(cfg, data_info, secfield, lon_s4new, lat_s4new)

    datafile.close()


def tsplot_extract_data(mmodel, observations, metadata_t, metadata_s, ind):
    """Extracts level data from the files for TS plots."""

    if mmodel != observations:
        level_pp = metadata_t['datafile'].variables['thetao'][0, ind, :, :]
        level_pp_s = metadata_s['datafile'].variables['so'][0, ind, :, :]
    else:
        level_pp = metadata_t['datafile'].variables['thetao'][0, ind, :, :]
        level_pp_s = metadata_s['datafile'].variables['so'][0, ind, :, :]
    # This is fix fo make models with 0 as missing values work,
    # should be fixed in fixes that do not work for now in the new backend
    if not isinstance(level_pp, np.ma.MaskedArray):
        level_pp = np.ma.masked_equal(level_pp, 0)
        level_pp_s = np.ma.masked_equal(level_pp_s, 0)
    return level_pp, level_pp_s


def tsplot_save_data(cfg, data_info, temp, salt, depth_model):
    """Save data for TS plots."""

    ofiles = {}
    data_info['variable'] = 'thetao'
    ofiles['ofilename_t'] = genfilename(**data_info, data_type='tsplot')
    np.save(ofiles['ofilename_t'], temp)
    provenance_record = get_provenance_record(data_info, 'tsplot', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_t'] + '.npy',
                              provenance_record)

    data_info['variable'] = 'so'
    ofiles['ofilename_s'] = genfilename(**data_info, data_type='tsplot')
    np.save(ofiles['ofilename_s'], salt)
    provenance_record = get_provenance_record(data_info, 'tsplot', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_s'] + '.npy',
                              provenance_record)

    data_info['variable'] = 'depth'
    ofiles['ofilename_depth'] = genfilename(**data_info, data_type='tsplot')
    np.save(ofiles['ofilename_depth'], depth_model)
    provenance_record = get_provenance_record(data_info, 'tsplot', 'npy')
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(ofiles['ofilename_depth'] + '.npy',
                              provenance_record)


def tsplot_data(cfg, mmodel, region, observations='PHC'):
    """Extract data for TS plots from one specific model.

    Parameters
    ----------

    mmodel: str
        model name
    max_level: int
        maximum level (depth) of TS data to be used
    region: str
        region as defined in `hofm_regions`
    observations: str
        name of the observations

    Returns
    -------
    None
    """
    logger.info("Extract  TS data for %s, region %s", mmodel, region)

    # generate input names for T and S. The files are generated by the
    # `timmean` function.
    ifilename_t = genfilename(cfg['work_dir'],
                              'thetao',
                              mmodel,
                              data_type='timmean',
                              extension='.nc')
    ifilename_s = genfilename(cfg['work_dir'],
                              'so',
                              mmodel,
                              data_type='timmean',
                              extension='.nc')
    # get the metadata for T and S

    metadata_t = load_meta(datapath=ifilename_t, fxpath=None)
    metadata_s = load_meta(datapath=ifilename_s, fxpath=None)

    # find index of the max_level
    lev_limit = metadata_t['lev'][
        metadata_t['lev'] <= cfg['tsdiag_depth']].shape[0] + 1
    # find indexes of data that are in the region
    indexes = hofm_regions(region, metadata_t['lon2d'], metadata_t['lat2d'])

    temp = np.array([])
    salt = np.array([])
    depth_model = np.array([])
    # loop over depths
    for ind, depth in enumerate(metadata_t['lev'][0:lev_limit]):
        level_pp, level_pp_s = tsplot_extract_data(mmodel, observations,
                                                   metadata_t, metadata_s, ind)
        # select individual points for T, S and depth
        temp = np.hstack((temp, level_pp[indexes[0], indexes[1]].compressed()))
        salt = np.hstack(
            (salt, level_pp_s[indexes[0], indexes[1]].compressed()))
        depth_temp = np.zeros_like(
            level_pp[indexes[0], indexes[1]].compressed())
        depth_temp[:] = depth
        depth_model = np.hstack((depth_model, depth_temp))

    # Saves the data to individual files
    data_info = {}
    data_info['basedir'] = cfg['work_dir']
    data_info['mmodel'] = mmodel
    data_info['region'] = region
    data_info['levels'] = metadata_t['lev']
    data_info['ori_file'] = [ifilename_t, ifilename_s]
    data_info['areacello'] = None
    tsplot_save_data(cfg, data_info, temp, salt, depth_model)

    metadata_t['datafile'].close()
    metadata_s['datafile'].close()


def aw_core(model_filenames, diagworkdir, region, cmor_var):
    """Calculate Atlantic Water (AW) core depth the region.

    The AW core is defined as water temperature maximum
    between 200 and 1000 meters. Can be in future generalised
    to find the depth of specific water masses.

    The function relies on the data for the profiles, so
    this information should be available.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    diagworkdir: str
        path to work directory.
    region: str
        one of the regions from `hofm_regions`,
        the data from the mean vertical profiles should be available.
    cmor_var: str
        name of the variable.

    Returns
    -------
    aw_core_parameters: dict
        For each model there is maximum temperature, depth level in the model,
        index of the depth level in the model.
    """
    logger.info("Calculate AW core statistics")
    aw_core_parameters = {}

    for mmodel in model_filenames:
        aw_core_parameters[mmodel] = {}
        logger.info("Plot profile %s data for %s, region %s", cmor_var, mmodel,
                    region)
        ifilename = genfilename(diagworkdir, cmor_var, mmodel, region, 'hofm',
                                '.npy')
        ifilename_levels = genfilename(diagworkdir, cmor_var, mmodel, region,
                                       'levels', '.npy')

        hofdata = np.load(ifilename, allow_pickle=True)
        lev = np.load(ifilename_levels, allow_pickle=True)

        profile = (hofdata)[:, :].mean(axis=1)
        maxvalue = np.max(profile[(lev >= 200) & (lev <= 1000)])
        maxvalue_index = np.where(profile == maxvalue)[0][0]
        maxvalue_depth = lev[maxvalue_index]

        if maxvalue > 100:
            maxvalue = maxvalue - 273.15

        aw_core_parameters[mmodel]['maxvalue'] = maxvalue
        aw_core_parameters[mmodel]['maxvalue_index'] = maxvalue_index
        aw_core_parameters[mmodel]['maxvalue_depth'] = maxvalue_depth

    return aw_core_parameters
