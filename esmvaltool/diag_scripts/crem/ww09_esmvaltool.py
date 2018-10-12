"""
Cloud Regime Error Metrics (CREM).

  Author: Keith Williams (Metoffice, UK)

  Project: ESA-CMUG

  Description
    Calculates the Cloud Regime Error Metric (CREM) following Williams and
    Webb (2009, Clim. Dyn.)

  Required diag_script_info attributes (diagnostics specific)
    none

  Optional diag_script_info attributes (diagnostic specific)
    none

  Required variable_info attributes (variable specific)
    none

  Optional variable_info attributes (variable specific)
    none

  Caveats
    TO DO:
      1) add metadata to plot
      2) add metadata to netcdf output
      3) use preprocessor for regridding input data

  Modification history
    20181012-A_laue_ax: extended (optional) netCDF output
    20180920-A_laue_ax: code adapted for ESMValTool v2.0
    20171128-A_laue_ax: added author and diagname to meta data
                        switched off "replacing of exact values"
                        in regridding function
    20170713-A_laue_ax: added tagging (for reporting)
    20151117-A_laue_ax: added parameters for call to "write_references"
    20151113-A_laue_ax: added creation of directory for plots if needed
                        (code was crashing if directory does not exist)
    20151029-A_laue_ax: added output of acknowledgements + processed files
                        to log-file
    20150903-A_laue_ax: ESMValTool implementation.
    20150521-A_will_ke: CREM routines written.
"""
import sys
import logging
import os
import numpy as np
from scipy.ndimage.interpolation import map_coordinates as interp2d
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata)

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
    vartags = []
    modeltags = []

    for var in ww_vars:
        vartags.append("V_" + var)

    # create list of dataset names (plot labels)
    models = []

    i = 0
    missing_vars = []

    for dataset in grouped_input_data:
        models.append(dataset)
        modeltags.append('M_' + dataset)

        pointers = {}

        for var in ww_vars:
            selection = select_metadata(input_data, dataset=dataset,
                                        short_name=var)
            if not selection:
                missing_vars.append(var)
            else:
                key = var + '_nc'
                pointers[key] = selection[0]['filename']

        # snow variable: use 'snc' if available or alternatively use 'snw'

        missing_snow = True

        for var in ww_vars_plus:
            selection = select_metadata(input_data, dataset=dataset,
                                        short_name=var)
            key = var + '_nc'
            if not selection:
                logger.info("%s: no data for variable snc found, trying "
                            "variable snw instead", dataset)
                pointers[key] = ""
            else:
                pointers[key] = selection[0]["filename"]
                vartags.append('V_' + var)
                missing_snow = False
                break

        if missing_snow:
            missing_vars.append(ww_vars_plus[0] + " or " + ww_vars_plus[1])

        for filen in pointers:
            climofiles.append(','.join(filen))

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

        # logger.info("==================================")
        # logger.info(dataset)
        # logger.info("==================================")
        # logger.info(crem_pd)
        # logger.info("..................................")
        j = 0
        for region in regions:
            # logger.info("*** " + region + " ***")
            regime = regions[region]
            k = 0
            for reg in regime:
                idx = allregimes.index(reg)
                r_crems[i, j, idx] = r_crem_pd[j, k]
                # printstr = ": %f" % r_crem_pd[j, k]
                # logger.info("  * " + reg + printstr)
                k = k + 1
            j = j + 1
        # logger.info("==================================")

        i = i + 1

    logger.info("==================================")
    logger.info("*** Cloud Regime Error Metrics ***")
    logger.info("==================================")
    logger.info(crems)
    logger.info("==================================")

    # plot results

    if cfg['write_plots']:
        oname = os.path.join(
            cfg['plot_dir'],
            'ww09_metric_multimodel.' + cfg['output_file_type'],
        )
        logger.debug("Plotting results to %s", oname)

        plt.figure()
        ypos = np.arange(nummod)
        plt.barh(ypos, crems, align='center')
        plt.yticks(ypos, models)
        plt.xlabel('Cloud Regime Error Metric')

        # draw observational uncertainties (dashed red line)
        plt.plot([0.96, 0.96], [-0.5, nummod - 0.5], 'r--')

        plt.savefig(oname, bbox_inches='tight')

        # add meta data to plot (for reporting)

#        basetags = 'TO BE DONE'
#
#        ESMValMD("both",
#            oname,
#            basetags + ['DM_global', 'PT_bar'] + modeltags + vartags,
#            'Cloud Regime Error Metric (CREM) following Williams and Webb '
#            '(2009, Clim. Dyn.).',
#            '#ID_ww09_crem',
#            ','.join(climofiles), 'ww09_ESMValTool.py', 'A_will_ke')

    if cfg['write_netcdf']:
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
        data = ncfile.createVariable('crem', np.dtype('float32').char,
                                     ('model'))
        r_data = ncfile.createVariable('r_crem', np.dtype('float32').char,
                                       ('model', 'region', 'regime'),
                                       fill_value=999.9)
        mod = ncfile.createVariable('model', np.dtype('int32').char,
                                    ('model'))
        reg = ncfile.createVariable('region', np.dtype('int32').char,
                                    ('region'))
        rgm = ncfile.createVariable('regime', np.dtype('int32').char,
                                    ('regime'))
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


# Reading and Regridding functions (scroll down for main program)


def regrid(a_in, x_in, y_in, x_out, y_out, fixmdis=True, x_cyclic=0.0):
    """
    Function for regridding.

    Regridding data onto 2.5 degree lat-long grid as the ISCCP
    obs data used for comparison was stored on.

    a_in : float
        input data
    x_in : float
        x coordinates (longitudes) of input data
    y_in : float
        y coordinates (latitudes) of input data
    x_out : float
        x coordinates (longitudes) of target grid
    y_out : float
        y coordinates (latitudes) of target grid
    fixmdis : Bool
        post-process regridded results: replace with original values for
        any "exact" coordinate matches
    x_cyclic : float
        xxxxx
    """
    # first represent missing data as np.NAN
    # - this replicates the default "hard MDI" behaviour of IDL regrid
    # code used by WW09 (in conjunction with the post-regrid "put back
    # exact coord matches" - see last part)
    a_in = a_in.copy()  # avoid overwriting
    if isinstance(a_in, np.ma.masked_array):
        a_in[np.ma.getmaskarray(a_in)] = np.NAN

    # replicate a column to the right if we have "x-cyclic" data

    # copy inputs to avoid changing them
    x_in = x_in.copy()
    y_in = y_in.copy()

    # sort the input Xs and Ys to guarantee ascending order
    n_x = len(x_in)
    n_y = len(y_in)
    i_sort_x = np.argsort(x_in)
    x_in = np.array(x_in)[i_sort_x]
    a_in = a_in[:, i_sort_x]
    i_sort_y = np.argsort(y_in)
    y_in = np.array(y_in)[i_sort_y]
    a_in = a_in[i_sort_y, :]

    # simulate cyclic X-coords, if enabled
    if x_cyclic > 0.0:
        a_inew = list(range(n_x)) + [0]   # Python 2-->3: range-->list(range)
        n_x += 1
        a_in = a_in[:, a_inew]   # recopy one lhs column on rhs
        x_in = x_in[a_inew]      # ditto for coords
        x_in[-1] += x_cyclic    # bump last element by range

    # convert input+output coordinate specs to "fractional coordinate values"
    xinds = np.interp(x_out, x_in, range(n_x))
    yinds = np.interp(y_out, y_in, range(n_y))

    # make a full coordinate mesh
    ainds = np.meshgrid(xinds, yinds)
    ainds = np.array(ainds)
    ainds = ainds[[1, 0]]

    # do main interpolation
    result = interp2d(a_in, ainds, order=1, mode='nearest', cval=np.NAN,
                      prefilter=False)
    # 1st-order spline is just bilinear interpolation

    # post-process replacing originals for any "exact" coordinate matches
    if fixmdis:
        bx_exact = abs(xinds - np.round(xinds, 0)) < 1e-6
        i_xout_exact = np.arange(n_x)[bx_exact]
        i_xin_exact = [int(round(ix)) for ix in xinds[i_xout_exact]]

        by_exact = abs(yinds - np.round(yinds, 0)) < 1e-6
        i_yout_exact = np.arange(n_y)[by_exact]
        i_yin_exact = [int(round(iy)) for iy in yinds[i_yout_exact]]

        for (i, ix_out) in enumerate(i_xout_exact):
            for (j, iy_out) in enumerate(i_yout_exact):
                result[iy_out, ix_out] = a_in[i_yin_exact[j], i_xin_exact[i]]

    return result


def read_and_regrid(srcfilename, varname, lons2, lats2):
    """
    Function for reading and regridding cmor compliant input data.

    Parameters
    ----------
    srcfilename : str
        filename containing input data
    varname : str
        variable name in netcdf
    lons2 : float
        longitudes of target grid
    lats2 : float
        latitudes of target grid
    """
    npts = len(lons2)
    nrows = len(lats2)

    n_time = len(Dataset(srcfilename, 'r').variables['time'][:])
    data_rg = np.zeros((n_time, nrows, npts))
    logger.debug('Number of data times in file %i', n_time)

    # read data
    src_dataset = Dataset(srcfilename, 'r')
    src_data = src_dataset.variables[varname]

    # grid of input data
    lats = src_dataset.variables['lat'][:]
    lons = src_dataset.variables['lon'][:]

    # create mask (missing values)
    data = np.ma.masked_equal(src_data, getattr(src_data, "_FillValue"))

    for i_t in range(n_time):    # range over fields in the file
        data_rg[i_t, :, :] = regrid(data[i_t, :, :], lons, lats, lons2, lats2,
                                    False, x_cyclic=360.0)

    rgmasked = np.ma.masked_invalid(data_rg)
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

    usrnames = ['albisccp', 'pctisccp', 'cltisccp', 'rsut', 'rsutcs', 'rlut',
                'rlutcs', 'snc', 'sic']

    varnames = usrnames[:]    # names used for nc vars.

    if not pointers['snc_nc']:
        varnames[7] = 'snw'

    # target grid spec
    npts = 144
    nrows = 72
    z_x = -1.25
    d_x = 2.5
    z_y = -91.25
    d_y = 2.5

    lons2 = np.array([z_x + d_x * (i + 1.0) for i in range(npts)])
    lats2 = np.array([z_y + d_y * (j + 1.0) for j in range(nrows)])

    # Read in and regrid input data

    logger.debug('Reading and regridding albisccp_nc')
    albisccp_data = read_and_regrid(pointers['albisccp_nc'], varnames[0],
                                    lons2, lats2)
#    E.add_to_filelist(pointers['albisccp_nc'])
    logger.debug('Reading and regridding pctisccp_nc')
    pctisccp_data = read_and_regrid(pointers['pctisccp_nc'], varnames[1],
                                    lons2, lats2)
#    E.add_to_filelist(pointers['pctisccp_nc'])
    logger.debug('Reading and regridding cltisccp_nc')
    cltisccp_data = read_and_regrid(pointers['cltisccp_nc'], varnames[2],
                                    lons2, lats2)
#    E.add_to_filelist(pointers['cltisccp_nc'])
    logger.debug('Reading and regridding rsut_nc')
    rsut_data = read_and_regrid(pointers['rsut_nc'], varnames[3],
                                lons2, lats2)
#    E.add_to_filelist(pointers['rsut_nc'])
    logger.debug('Reading and regridding rsutcs_nc')
    rsutcs_data = read_and_regrid(pointers['rsutcs_nc'], varnames[4],
                                  lons2, lats2)
#    E.add_to_filelist(pointers['rsutcs_nc'])
    logger.debug('Reading and regridding rlut_nc')
    rlut_data = read_and_regrid(pointers['rlut_nc'], varnames[5],
                                lons2, lats2)
#    E.add_to_filelist(pointers['rlut_nc'])
    logger.debug('Reading and regridding rlutcs_nc')
    rlutcs_data = read_and_regrid(pointers['rlutcs_nc'], varnames[6],
                                  lons2, lats2)
#    E.add_to_filelist(pointers['rlutcs_nc'])
    logger.debug('Reading and regridding sic_nc')
    sic_data = read_and_regrid(pointers['sic_nc'], varnames[8],
                               lons2, lats2)
#    E.add_to_filelist(pointers['sic_nc'])
    if not pointers['snc_nc']:
        logger.debug('Reading and regridding snw_nc')
        snc_data = read_and_regrid(pointers['snw_nc'], varnames[7],
                                   lons2, lats2)
#        E.add_to_filelist(pointers['snw_nc'])
    else:
        logger.debug('Reading and regridding snc_nc')
        snc_data = read_and_regrid(pointers['snc_nc'], varnames[7],
                                   lons2, lats2)
#        E.add_to_filelist(pointers['snc_nc'])

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
