"""
Cloud Regime Error Metrics (CREM).

  Author: Keith Williams (Metoffice, UK)

  Project: ESA-CMUG

  Description
    Calculates the Cloud Regime Error Metric (CREM) following Williams and
    Webb (2009, Clim. Dyn.). Regridding to the 2.5x2.5 degree ISCCP grid is
    done by the ESMValTool preprocessor.

  Required diag_script_info attributes (diagnostics specific)
    none

  Optional diag_script_info attributes (diagnostic specific)
    none

  Required variable_info attributes (variable specific)
    none

  Optional variable_info attributes (variable specific)
    none

  Caveats
    none

  Modification history
    20190216-lauer_axel: outsourced regridding to preprocessor
    20190215-lauer_axel: added metadata to netcdf output and plot
    20190213-lauer_axel: made code more flexible to support CMIP6 data
    20181012-lauer_axel: extended (optional) netCDF output
    20180920-lauer_axel: code adapted for ESMValTool v2.0
    20171128-lauer_axel: added author and diagname to meta data
                         switched off "replacing of exact values"
                         in regridding function
    20170713-lauer_axel: added tagging (for reporting)
    20151117-lauer_axel: added parameters for call to "write_references"
    20151113-lauer_axel: added creation of directory for plots if needed
                         (code was crashing if directory does not exist)
    20151029-lauer_axel: added output of acknowledgements + processed files
                         to log-file
    20150903-lauer_axel: ESMValTool implementation.
    20150521-williams_keith: CREM routines written.
"""
import logging
import os
import sys
from pprint import pformat

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from esmvaltool.diag_scripts.shared import (
    group_metadata, ProvenanceLogger, run_diagnostic, select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    # get description of the preprocessed data
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'dataset')
    nummod = len(grouped_input_data)
    crems = np.empty(nummod)

    # list of variables needed for CREM calculations
    ww_vars = ('albisccp', 'pctisccp', 'cltisccp', 'rsut', 'rsutcs', 'rlut',
               'rlutcs', 'sic')
    ww_vars_plus = ('snc', 'snw')
    # alternative variable names to check if variable was not found (CMIP6)
    ww_vars_alternative = {'sic': 'siconc'}

    # for human readable output
    # regions/regimes as they come from the CREM calculation
    regions = {'tropics': ['shallow cumulus', 'congestus', 'thin cirrus',
                           'stratocumulus/cumulus transition',
                           'anvil cirrus', 'deep convection',
                           'stratocumulus'],
               'ice-free-extra-tropics': ['shallow cumulus', 'congestus',
                                          'stratocumulus/cumulus transition',
                                          'cirrus', 'stratocumulus',
                                          'frontal', 'thin cirrus'],
               'snow-ice-covered': ['shallow cumulus', 'stratocumulus',
                                    'thick mid-level', 'frontal',
                                    'thin mid-level', 'thin cirrus']}
    # regimes as we write them to netCDF
    allregimes = ['shallow cumulus', 'congestus', 'thin cirrus',
                  'stratocumulus/cumulus transition', 'anvil cirrus',
                  'deep convection', 'stratocumulus', 'cirrus',
                  'frontal', 'thick mid-level', 'thin mid-level']
    # field for (optional) netCDF output of individual regions and regimes
    r_crems = np.empty((nummod, len(regions), len(allregimes)))
    r_crems[:] = 999.9

    # provenance information
    climofiles = []

    # create list of dataset names (plot labels)
    models = []

    i = 0
    missing_vars = []

    for dataset in grouped_input_data:
        models.append(dataset)
        pointers = {}

        for var in ww_vars:
            selection = select_metadata(input_data, dataset=dataset,
                                        short_name=var)
            alt_var = None
            if not selection:
                # try alternative variable name (if defined)
                if var in ww_vars_alternative:
                    alt_var = ww_vars_alternative[var]
                    selection = select_metadata(input_data, dataset=dataset,
                                                short_name=alt_var)
            if not selection:
                missing_vars.append(var)
            else:
                key_nc = var + '_nc'
                key_var = var
                pointers[key_nc] = selection[0]['filename']
                if alt_var is None:
                    pointers[key_var] = var
                else:
                    pointers[key_var] = alt_var

        # snow variable: use 'snc' if available or alternatively use 'snw'

        missing_snow = True

        for var in ww_vars_plus:
            selection = select_metadata(input_data, dataset=dataset,
                                        short_name=var)
            key_nc = var + '_nc'
            key_var = var
            if not selection:
                logger.info("%s: no data for variable snc found, trying "
                            "variable snw instead", dataset)
                pointers[key_nc] = ""
                pointers[key_var] = ""
            else:
                pointers[key_nc] = selection[0]["filename"]
                pointers[key_var] = var
                missing_snow = False
                break

        if missing_snow:
            missing_vars.append(ww_vars_plus[0] + " or " + ww_vars_plus[1])

        for key in pointers:
            if key[-3:] == '_nc':
                climofiles.append(pointers[key])

        # check if all variables are available

        if missing_vars:
            printlist = ', '.join(missing_vars)
            logger.error("error: the following variables are not "
                         "available: %s", printlist)
            raise Exception('Variables missing (see log file for details).')

        # calculate CREM

        (crem_pd, r_crem_pd) = crem_calc(pointers)

        crems[i] = crem_pd

        # sort results into output array

        j = 0
        for region in regions:
            regime = regions[region]
            k = 0
            for reg in regime:
                idx = allregimes.index(reg)
                r_crems[i, j, idx] = r_crem_pd[j, k]
                k = k + 1
            j = j + 1

        i = i + 1

    logger.info("==================================")
    logger.info("*** Cloud Regime Error Metrics ***")
    logger.info("==================================")
    logger.info(crems)
    logger.info("==================================")

    # define diagnostic internal provenance data

    provenance_record = {
        'caption': 'Cloud Regime Error Metric (CREM) following Williams ' +
                   'and Webb (2009, Clim. Dyn.).',
        'statistics': ['other'],
        'domains': ['global'],
        'plot_type': 'bar',
        'authors': [
            'williams_keith',
            'lauer_axel',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': climofiles,
    }

    # plot results

    if cfg['write_plots']:
        plotname = os.path.join(
            cfg['plot_dir'],
            'ww09_metric_multimodel.' + cfg['output_file_type'],
        )
        logger.debug("Plotting results to %s", plotname)

        plt.figure()
        ypos = np.arange(nummod)
        plt.barh(ypos, crems, align='center')
        plt.yticks(ypos, models)
        plt.xlabel('Cloud Regime Error Metric')

        # draw observational uncertainties (dashed red line)
        plt.plot([0.96, 0.96], [-0.5, nummod - 0.5], 'r--')

        plt.savefig(plotname, bbox_inches='tight')

        provenance_record['plot_file'] = plotname

    # save results to netcdf

    oname = os.path.join(cfg['work_dir'], 'ww09_metric_multimodel.nc')
    logger.debug("Saving results to %s", oname)
    # convert strings
    modstr_out = np.array(models, dtype=object)
    regionstr_out = np.array(list(regions.keys()), dtype=object)
    regimestr_out = np.array(allregimes, dtype=object)
    # open a new netCDF file for writing
    ncfile = Dataset(oname, 'w')
    # create dimensions
    ncfile.createDimension('model', nummod)
    ncfile.createDimension('region', len(regions))
    ncfile.createDimension('regime', len(allregimes))
    # create variables
    data = ncfile.createVariable('crem', np.dtype('float32').char, ('model'))
    r_data = ncfile.createVariable('r_crem', np.dtype('float32').char,
                                   ('model', 'region', 'regime'),
                                   fill_value=999.9)
    mod = ncfile.createVariable('model', np.dtype('int32').char, ('model'))
    reg = ncfile.createVariable('region', np.dtype('int32').char, ('region'))
    rgm = ncfile.createVariable('regime', np.dtype('int32').char, ('regime'))
    mod_name = ncfile.createVariable('model_name', str, ('model'))
    reg_name = ncfile.createVariable('region_name', str, ('region'))
    rgm_name = ncfile.createVariable('regime_name', str, ('regime'))
    # write data to variable
    data[:] = crems
    r_data[:, :, :] = r_crems
    mod[:] = range(nummod)
    reg[:] = range(len(regions))
    rgm[:] = range(len(allregimes))
    mod_name[:] = modstr_out
    reg_name[:] = regionstr_out
    rgm_name[:] = regimestr_out
    # close the file
    ncfile.close()

    # add provenance data to netcdf and plot

    logger.info("Recording provenance of %s:\n%s", oname,
                pformat(provenance_record))

    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(oname, provenance_record)


def read_and_check(srcfilename, varname, lons2, lats2, time2):
    """
    Function for reading and checking for correct regridding of input data.

    Parameters
    ----------
    srcfilename : str
        filename containing input data
    varname : str
        variable name in netcdf
    lons2 : float
        longitudes of target grid (ISCCP)
    lats2 : float
        latitudes of target grid (ISCCP)
    time2: integer
        number of time steps
    """
    nlon = len(lons2)
    nlat = len(lats2)

    src_dataset = Dataset(srcfilename, 'r')

    n_time = len(src_dataset.variables['time'][:])
    logger.debug('Number of data times in file %s is %i', srcfilename, n_time)

    # check number of time steps is matching

    if n_time != time2:
        logger.error("error: number of time steps in input files are "
                     "not equal")
        raise Exception('Variables contain different number of time steps '
                        '(see log file for details).')

    grid_mismatch = False
    coord_mismatch = False

    # check longitudes

    lons = src_dataset.variables['lon'][:]
    if nlon != len(lons):
        grid_mismatch = True
    if np.amax(np.absolute(lons - lons2)) > 1.0e-3:
        coord_mismatch = True

    # check latitudes

    lats = src_dataset.variables['lat'][:]
    if nlat != len(lats):
        grid_mismatch = True
    if np.amax(np.absolute(lats - lats2)) > 1.0e-3:
        coord_mismatch = True

    if grid_mismatch:
        logger.error("error: input data are not on 2.5x2.5 deg ISCCP grid."
                     "lons = %i (required: %i), lats = %i (required: %i)",
                     len(lons), nlon, len(lats), nlat)

    if coord_mismatch:
        logger.error("error: input data are not on 2.5x2.5 deg ISCCP grid, "
                     "longitudes and/or latitudes differ from ISCCP grid by "
                     "more than 1.0e-3")

    if (grid_mismatch or coord_mismatch):
        raise Exception('Input variables are not on 2.5x2.5 deg ISCCP grid '
                        '(see log file for details).')

    # read data
    src_data = src_dataset.variables[varname]

    # create mask (missing values)
    try:
        data = np.ma.masked_equal(src_data, getattr(src_data, "_FillValue"))
        rgmasked = np.ma.masked_invalid(data)
    except AttributeError:
        rgmasked = np.ma.masked_invalid(src_data)
    np.ma.set_fill_value(rgmasked, 0.0)

    return np.ma.filled(rgmasked)


def crem_calc(pointers):
    """
    Main program for calculating Cloud Regime Error Metric.

    Following equation 4 in Williams and Webb (2009) (WW09).

    Parameters
    ----------
    pointers : dict
        Keys in dictionary are: albisccp_nc, pctisccp_nc, cltisccp_nc,
        rsut_nc, rsutcs_nc, rlut_nc, rlutcs_nc, snc_nc, sic_nc

    For CMIP5, snc is in the CMIP5 table 'day'. All other variables
    are in the CMIP5 table 'cfday'. A minimum of 2 years, and ideally 5
    years, of data are required. The observational regime characteristics
    were calculated for the period Mar 1985 - Feb 1990.

    If snc is not available then snw can be used instead. In this case
    pointers[snc_nc] should be set to None and snw_nc set.

    Returns
    -------
    crem_pd : float
        present-day cloud regime error metric of WW09.
    r_crem_pd : float
        component from each regime.
    """
    # Lookup arrays
    # Observational regime centroids for assignment of the model data.
    # These are taken from Table 3 of Williams and Webb (2009)
    # (999.9 represents missing data). The observational regime
    # characteristics were calculated for the period Mar 1985 - Feb 1990.

    obs_alb = np.array([[0.261, 0.339, 0.211, 0.338, 0.313, 0.532, 0.446],
                        [0.286, 0.457, 0.375, 0.325, 0.438, 0.581, 0.220],
                        [0.433, 0.510, 0.576, 0.505, 0.343, 0.247, 999.9]])

    obs_pct = np.array([[0.652, 0.483, 0.356, 0.784, 0.327, 0.285, 0.722],
                        [0.643, 0.607, 0.799, 0.430, 0.723, 0.393, 0.389],
                        [0.582, 0.740, 0.620, 0.458, 0.595, 0.452, 999.9]])

    obs_clt = np.array([[0.314, 0.813, 0.740, 0.640, 0.944, 0.979, 0.824],
                        [0.473, 0.932, 0.802, 0.914, 0.900, 0.978, 0.713],
                        [0.356, 0.747, 0.778, 0.884, 0.841, 0.744, 999.9]])

    # Observed regime RFO's taken from Table 3 of WW09
    obs_rfo = np.array([[0.375, 0.195, 0.119, 0.103, 0.091, 0.064, 0.052],
                        [0.354, 0.170, 0.114, 0.104, 0.091, 0.083, 0.083],
                        [0.423, 0.191, 0.139, 0.111, 0.094, 0.042, 999.9]])

    # Observed regime net cloud forcing (Figure 2f of WW09)
    obs_ncf = np.array([[-10.14, -25.45, -5.80, -27.40, -16.83, -48.45,
                         -55.84],
                        [-13.67, -58.28, -36.26, -25.34, -64.27, -56.91,
                         -11.63],
                        [-3.35, -16.66, -13.76, -8.63, -12.17, 1.45, 999.9]])

    # aw in eq 3 of WW09
    area_weights = np.array([0.342, 0.502, 0.156])
    # weighting for swcf to account for lack of ISCCP diagnostics
    # during polar night (p153 of WW09)
    solar_weights = np.array([1.000, 0.998, 0.846])

    # number of regimes in each region (Table 3 of WW09)
    nregimes = {'tropics': 7, 'extra-tropics': 7, 'snow-ice': 6}

    # -----------------------------------------------------------

    # Section to re-grid onto 2.5 degr lat long grid.
    # Note this has been tested with regular lat-long grids - other grid
    # types may need changes to the regrid subroutine.

    # target grid spec
    npts = 144
    nrows = 72
    z_x = -1.25
    d_x = 2.5
    z_y = -91.25
    d_y = 2.5

    lons2 = np.array([z_x + d_x * (i + 1.0) for i in range(npts)])
    lats2 = np.array([z_y + d_y * (j + 1.0) for j in range(nrows)])

    # Read input data
    # ---------------
    # pointers['xxx_nc'] = file name of input file
    # pointers['xxx'] = actual variable name in input file

    logger.debug('Reading albisccp')
    ntime2 = len(Dataset(pointers['albisccp_nc'], 'r').variables['time'][:])
    albisccp_data = read_and_check(pointers['albisccp_nc'],
                                   pointers['albisccp'], lons2, lats2, ntime2)
    logger.debug('Reading pctisccp')
    pctisccp_data = read_and_check(pointers['pctisccp_nc'],
                                   pointers['pctisccp'], lons2, lats2, ntime2)
    logger.debug('Reading cltisccp')
    cltisccp_data = read_and_check(pointers['cltisccp_nc'],
                                   pointers['cltisccp'], lons2, lats2, ntime2)
    logger.debug('Reading rsut')
    rsut_data = read_and_check(pointers['rsut_nc'],
                               pointers['rsut'], lons2, lats2, ntime2)
    logger.debug('Reading rsutcs')
    rsutcs_data = read_and_check(pointers['rsutcs_nc'],
                                 pointers['rsutcs'], lons2, lats2, ntime2)
    logger.debug('Reading rlut')
    rlut_data = read_and_check(pointers['rlut_nc'],
                               pointers['rlut'], lons2, lats2, ntime2)
    logger.debug('Reading rlutcs')
    rlutcs_data = read_and_check(pointers['rlutcs_nc'],
                                 pointers['rlutcs'], lons2, lats2, ntime2)
    logger.debug('Reading sic')
    sic_data = read_and_check(pointers['sic_nc'],
                              pointers['sic'], lons2, lats2, ntime2)
    if not pointers['snc_nc']:
        logger.debug('Reading snw')
        snc_data = read_and_check(pointers['snw_nc'],
                                  pointers['snw'], lons2, lats2, ntime2)
    else:
        logger.debug('Reading snc')
        snc_data = read_and_check(pointers['snc_nc'],
                                  pointers['snc'], lons2, lats2, ntime2)

    # -----------------------------------------------------------

    # Set up storage arrays
    numreg = len(nregimes)            # = 3
    numrgm = nregimes[max(nregimes)]  # = 7

    model_rfo = np.zeros((numreg, numrgm))
    model_ncf = np.zeros((numreg, numrgm))
    r_crem_pd = np.zeros((numreg, numrgm))
    model_rfo[:] = 999.9
    model_ncf[:] = 999.9
    r_crem_pd[:] = 999.9

    # Normalize data used for assignment to regimes to be in the range 0-1
    pctisccp_data = pctisccp_data / 100000.0
    cltisccp_data = cltisccp_data / 100.0

    # Calculate cloud forcing
    swcf_data = rsutcs_data - rsut_data
    lwcf_data = rlutcs_data - rlut_data

    # loop over 3 regions
    # (0 = tropics, 1 = ice-free extra-tropics, 2 = snow/ice covered)
    for idx_region, (region, regime) in enumerate(nregimes.items()):

        # Set up validity mask for region

        mask = pctisccp_data.copy()
        if region == 'tropics':
            mask[:, (lats2 < -20) | (lats2 > 20), :] = np.NAN
        elif region == 'extra-tropics':
            mask[:, (lats2 >= -20) & (lats2 <= 20), :] = np.NAN
            mask[(snc_data >= 0.1) | (sic_data >= 0.1)] = np.NAN
        elif region == 'snow-ice':
            mask[:, (lats2 >= -20) & (lats2 <= 20), :] = np.NAN
            mask[(snc_data < 0.1) & (sic_data < 0.1)] = np.NAN

        mask[cltisccp_data == 0.0] = np.NAN

        points = np.isfinite(mask)
        npoints = len(mask[points])  # Number of valid data points in region

        group = np.zeros(npoints)
        e_d = np.zeros((npoints, regime))

        swcf_data_pts = swcf_data[points]
        lwcf_data_pts = lwcf_data[points]

        # Assign model data to observed regimes

        for i in range(regime):
            e_d[:, i] = \
                ((albisccp_data[points] - obs_alb[idx_region, i]) ** 2) + \
                ((pctisccp_data[points] - obs_pct[idx_region, i]) ** 2) + \
                ((cltisccp_data[points] - obs_clt[idx_region, i]) ** 2)

        group[:] = np.argmin(e_d, axis=1)

        for i in range(regime):
            mem = (group == i)

            count = len(group[mem])

            if count > 0:

                model_rfo[idx_region, i] = float(count) / float(npoints)
                model_ncf[idx_region, i] = np.average(swcf_data_pts[mem]) \
                    * solar_weights[idx_region] +                         \
                    np.average(lwcf_data_pts[mem])
            else:
                logger.info("Model does not reproduce all observed cloud "
                            "regimes.")
                logger.info("Cannot calculate CREM. Abort.")
                sys.exit()
                model_rfo[idx_region, i] = 0.0
                model_ncf[idx_region, i] = 0.0

    # Calculation of eq 3 in WW09
    for idx_region, (region, regime) in enumerate(nregimes.items()):
        r_crem_pd[idx_region, 0:regime] = area_weights[idx_region] * \
            (((model_ncf[idx_region, 0:regime] -
               obs_ncf[idx_region, 0:regime]) *
              obs_rfo[idx_region, 0:regime]) ** 2 +
             ((model_rfo[idx_region, 0:regime] -
               obs_rfo[idx_region, 0:regime]) *
              obs_ncf[idx_region, 0:regime]) ** 2) ** 0.5

    # Calculation of eq 4 in WW09
    crem_pd = ((np.sum(r_crem_pd[0, :] ** 2) + np.sum(r_crem_pd[1, :] ** 2) +
                np.sum(r_crem_pd[2, 0:5] ** 2)) / 20.0) ** 0.5

    # A perfect crem_pd with respect to ISCCP would be 0.0
    # An estimate of observational uncertainty (obtained by calculating
    # crem_pd wrt MODIS/ERBE) is 0.96 (i.e. models with crem_pd less than
    # 0.96 may be regarded as within observational uncertainty overall,
    # although not necessarily for every regime)'.
    # Interrogation of the r_crem_pd array from this program will indicate
    # which regimes contribute most to the total crem_pd (elements ordered
    # as Table 3 of WW09)'

    return crem_pd, r_crem_pd


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
