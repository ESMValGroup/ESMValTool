"""
; ############################################################################
; Cloud Regime Error Metrics (CREM)
; Author: Keith Williams (Metoffice, UK)
; ESA-CMUG project
; ############################################################################
; Description
;    Calculates the Cloud Regime Error Metric (CREM) following Williams and
;    Webb (2009, Clim. Dyn.)
;
; Required diag_script_info attributes (diagnostics specific)
;    none
;
; Optional diag_script_info attributes (diagnostic specific)
;    none
;
; Required variable_info attributes (variable specific)
;    none
;
; Optional variable_info attributes (variable specific)
;    none
;
; Caveats
;    TO DO:
;       1) add metadata to plot
;       2) add metadata to netcdf output
;       3) use preprocessor for regridding input data
;
; Modification history
;    20180920-A_laue_ax: code adapted for ESMValTool v2.0
;    20171128-A_laue_ax: added author and diagname to meta data
;                        switched off "replacing of exact values"
;                        in regridding function
;    20170713-A_laue_ax: added tagging (for reporting)
;    20151117-A_laue_ax: added parameters for call to "write_references"
;    20151113-A_laue_ax: added creation of directory for plots if needed
;                        (code was crashing if directory does not exist)
;    20151029-A_laue_ax: added output of acknowledgements + processed files
;                        to log-file
;    20150903-A_laue_ax: ESMValTool implementation.
;    20150521-A_will_ke: CREM routines written.
;
; ############################################################################
"""

import sys
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.interpolation import map_coordinates as interp2d
from netCDF4 import Dataset
from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
import matplotlib
matplotlib.use('Agg')

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
    vars = ('albisccp', 'pctisccp', 'cltisccp', 'rsut', 'rsutcs', 'rlut',
            'rlutcs', 'sic')
    vars_plus = ('snc', 'snw')

    # provenance information
    climofiles = []
    vartags = []
    modeltags = []

    for var in vars:
        vartags.append("V_" + var)

    # create list of dataset names (plot labels)
    models = []

    i = 0
    missing_vars = []

    for dataset in grouped_input_data:
        models.append(dataset)
        modeltags.append('M_' + dataset)

        pointers = {}

        for var in vars:
            selection = select_metadata(input_data, dataset=dataset,
                                        short_name=var)
            if not selection:
                missing_vars.append(var)
            else:
                key = var + '_nc'
                pointers[key] = selection[0]['filename']

        # snow variable: use 'snc' if available or alternatively use 'snw'

        missing_snow = True

        for var in vars_plus:
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
            missing_vars.append(vars_plus[0] + " or " + vars_plus[1])

        for fn in pointers:
            climofiles.append(','.join(fn))

        # check if all variables are available

        if missing_vars:
            printlist = ', '.join(missing_vars)
            logger.error("error: the following variables are not "
                         "available: %s", printlist)
            raise Exception('Variables missing (see log file for details).')

        # calculate CREM

        (CREMpd, __) = crem_calc(pointers)

        crems[i] = CREMpd
        i = i + 1

    logger.info("------------------------------------")
    logger.info(crems)
    logger.info("------------------------------------")

    # plot results

    if cfg['write_plots']:
        oname = os.path.join(
            cfg['plot_dir'],
            'ww09_metric_multimodelname.' + cfg['output_file_type'],
        )
        logger.debug("Plotting results to %s", oname)

        fig = plt.figure()
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
        oname = os.path.join(cfg['work_dir'], 'ww09_metric_multimodelname.nc')
        logger.debug("Saving results to %s", oname)
        # convert strings
        str_out = np.array(models, dtype=object)
        # open a new netCDF file for writing
        ncfile = Dataset(oname, 'w')
        # create the x dimension
        ncfile.createDimension('n', nummod)
        # create the variable
        # first argument is name of variable, second is datatype, third is
        # a tuple with the name(s) of dimension(s)
        data = ncfile.createVariable('crem', np.dtype('float32').char, ('n'))
        mods = ncfile.createVariable('model', str, ('n'))
        # write data to variable
        data[:] = crems
        mods[:] = str_out
        # close the file
        ncfile.close()


# Reading and Regridding functions (scroll down for main program)


def regrid(aIn, xIn, yIn, xOut, yOut, fixmdis=True, xCyclic=0.0):
    """
    Function for regridding onto 2.5 degree lat-long grid as the ISCCP
    obs data used for comparison was stored on.

    aIn : xx
        xxx
    xIn : xxx
        xxx
    yIn : xxx
        xxx
    xOut : xxx
        xxx
    yOut : xxx
        xxx
    fixmdis : Bool
        xxx
    xCyclic : float
        xxxxx
    """
    # first represent missing data as np.NAN
    # - this replicates the default "hard MDI" behaviour of IDL regrid
    # code used by WW09 (in conjunction with the post-regrid "put back
    # exact coord matches" - see last part)
    aIn = aIn.copy()  # avoid overwriting
    if isinstance(aIn, np.ma.masked_array):
        aIn[np.ma.getmaskarray(aIn)] = np.NAN

    # replicate a column to the right if we have "x-cyclic" data

    # copy inputs to avoid changing them
    xIn = xIn.copy()
    yIn = yIn.copy()

    # sort the input Xs and Ys to guarantee ascending order
    nx = len(xIn)
    ny = len(yIn)
    iSortX = np.argsort(xIn)
    xIn = np.array(xIn)[iSortX]
    aIn = aIn[:, iSortX]
    iSortY = np.argsort(yIn)
    yIn = np.array(yIn)[iSortY]
    aIn = aIn[iSortY, :]

    # simulate cyclic X-coords, if enabled
    if xCyclic > 0.0:
        aiNew = list(range(nx)) + [0]   # Python 2 --> 3: range --> list(range)
        nx += 1
        aIn = aIn[:, aiNew]   # recopy one lhs column on rhs
        xIn = xIn[aiNew]      # ditto for coords
        xIn[-1] += xCyclic    # bump last element by range

    # convert input+output coordinate specs to "fractional coordinate values"
    xinds = np.interp(xOut, xIn, range(nx))
    yinds = np.interp(yOut, yIn, range(ny))

    # make a full coordinate mesh
    ainds = np.meshgrid(xinds, yinds)
    ainds = np.array(ainds)
    ainds = ainds[[1, 0]]

    # qdo main interpolation
    result = interp2d(aIn, ainds, order=1, mode='nearest', cval=np.NAN,
                      prefilter=False)
    # 1st-order spline is just bilinear interpolation

    # post-process replacing originals for any "exact" coordinate matches
    if fixmdis:
        bXexact = abs(xinds - np.round(xinds, 0)) < 1e-6
        iXoutExact = np.arange(nx)[bXexact]
        iXinExact = [int(round(ix)) for ix in xinds[iXoutExact]]

        bYexact = abs(yinds - np.round(yinds, 0)) < 1e-6
        iYoutExact = np.arange(ny)[bYexact]
        iYinExact = [int(round(iy)) for iy in yinds[iYoutExact]]

        for (i, ixOut) in enumerate(iXoutExact):
            for (j, iyOut) in enumerate(iYoutExact):
                result[iyOut, ixOut] = aIn[iYinExact[j], iXinExact[i]]

    return result


def read_and_regrid(sSrcFilename, sVarname, lons2, lats2):
    """
    Function for reading and regridding cmor compliant input data.

    Parameters
    ----------
    sSrcFilename : str
        filename
    sVarname : str
        xxxxx
    lons2 : xxxx
        xxxxx
    lats2 : xxxxx
        xxxxxx
    """

    npts = len(lons2)
    nrows = len(lats2)

    nt = len(Dataset(sSrcFilename, 'r').variables['time'][:])
    data_rg = np.zeros((nt, nrows, npts))
    logger.debug('Number of data times in file %i', nt)

    # read data
    srcDataset = Dataset(sSrcFilename, 'r', format='NETCDF3')
    srcData = srcDataset.variables[sVarname]

    # grid of input data
    lats = srcDataset.variables['lat'][:]
    lons = srcDataset.variables['lon'][:]

    # create mask (missing values)
    data = np.ma.masked_equal(srcData, srcData._FillValue)

    for iT in range(nt):    # range over fields in the file
        data_rg[iT, :, :] = regrid(data[iT, :, :], lons, lats, lons2, lats2,
                                   False, xCyclic=360.0)

    rgmasked = np.ma.masked_invalid(data_rg)
    np.ma.set_fill_value(rgmasked, 0.0)

    return(np.ma.filled(rgmasked))


def crem_calc(pointers):
    """
    Main program for calculating Cloud Regime Error Metric following equation
    4 in Williams and Webb (2009) (WW09) from CMOR-compliant netCDF data.

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
    CREMpd : xxxx
        present-day cloud regime error metric of WW09.
    rCREMpd : xxxx
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
    nregimes = np.array([7, 7, 6])

    # -----------------------------------------------------------

    # Section to re-grid onto 2.5 degr lat long grid.
    # Note this has been tested with regular lat-long grids - other grid
    # types may need changes to the regrid subroutine.

    sUsrnames = ['albisccp', 'pctisccp', 'cltisccp', 'rsut', 'rsutcs', 'rlut',
                 'rlutcs', 'snc', 'sic']

    sVarnames = sUsrnames[:]    # names used for nc vars.

    if not pointers['snc_nc']:
        sVarnames[7] = 'snw'

    # target grid spec
    npts = 144
    nrows = 72
    zx = -1.25
    dx = 2.5
    zy = -91.25
    dy = 2.5

    lons2 = np.array([zx + dx * (i + 1.0) for i in range(npts)])
    lats2 = np.array([zy + dy * (j + 1.0) for j in range(nrows)])

    # Read in and regrid input data

    logger.debug('Reading and regridding albisccp_nc')
    albisccp_data = read_and_regrid(pointers['albisccp_nc'], sVarnames[0],
                                    lons2, lats2)
#    E.add_to_filelist(pointers['albisccp_nc'])
    logger.debug('Reading and regridding pctisccp_nc')
    pctisccp_data = read_and_regrid(pointers['pctisccp_nc'], sVarnames[1],
                                    lons2, lats2)
#    E.add_to_filelist(pointers['pctisccp_nc'])
    logger.debug('Reading and regridding cltisccp_nc')
    cltisccp_data = read_and_regrid(pointers['cltisccp_nc'], sVarnames[2],
                                    lons2, lats2)
#    E.add_to_filelist(pointers['cltisccp_nc'])
    logger.debug('Reading and regridding rsut_nc')
    rsut_data = read_and_regrid(pointers['rsut_nc'], sVarnames[3],
                                lons2, lats2)
#    E.add_to_filelist(pointers['rsut_nc'])
    logger.debug('Reading and regridding rsutcs_nc')
    rsutcs_data = read_and_regrid(pointers['rsutcs_nc'], sVarnames[4],
                                  lons2, lats2)
#    E.add_to_filelist(pointers['rsutcs_nc'])
    logger.debug('Reading and regridding rlut_nc')
    rlut_data = read_and_regrid(pointers['rlut_nc'], sVarnames[5],
                                lons2, lats2)
#    E.add_to_filelist(pointers['rlut_nc'])
    logger.debug('Reading and regridding rlutcs_nc')
    rlutcs_data = read_and_regrid(pointers['rlutcs_nc'], sVarnames[6],
                                  lons2, lats2)
#    E.add_to_filelist(pointers['rlutcs_nc'])
    logger.debug('Reading and regridding sic_nc')
    sic_data = read_and_regrid(pointers['sic_nc'], sVarnames[8],
                               lons2, lats2)
#    E.add_to_filelist(pointers['sic_nc'])
    if not pointers['snc_nc']:
        logger.debug('Reading and regridding snw_nc')
        snc_data = read_and_regrid(pointers['snw_nc'], sVarnames[7],
                                   lons2, lats2)
#        E.add_to_filelist(pointers['snw_nc'])
    else:
        logger.debug('Reading and regridding snc_nc')
        snc_data = read_and_regrid(pointers['snc_nc'], sVarnames[7],
                                   lons2, lats2)
#        E.add_to_filelist(pointers['snc_nc'])

    # -----------------------------------------------------------

    # Set up storage arrays
    model_rfo = np.zeros((3, 7))
    model_ncf = np.zeros((3, 7))
    rCREMpd = np.zeros((3, 7))
    model_rfo[:] = 999.9
    model_ncf[:] = 999.9
    rCREMpd[:] = 999.9

    # Normalize data used for assignment to regimes to be in the range 0-1
    pctisccp_data = pctisccp_data / 100000.0
    cltisccp_data = cltisccp_data / 100.0

    # Calculate cloud forcing
    swcf_data = rsutcs_data - rsut_data
    lwcf_data = rlutcs_data - rlut_data

    logger.debug('Assigning data to observational cloud regimes')

    # loop over 3 regions (tropics, extra tropics, snow/ice)
    for region in range(3):

        # Set up validity mask for region

        mask = pctisccp_data.copy()
        if region == 0:
            mask[:, (lats2 < -20) | (lats2 > 20), :] = np.NAN
        elif region == 1:
            mask[:, (lats2 >= -20) & (lats2 <= 20), :] = np.NAN
            mask[(snc_data >= 0.1) | (sic_data >= 0.1)] = np.NAN
        elif region == 2:
            mask[:, (lats2 >= -20) & (lats2 <= 20), :] = np.NAN
            mask[(snc_data < 0.1) & (sic_data < 0.1)] = np.NAN

        mask[cltisccp_data == 0.0] = np.NAN

        points = np.isfinite(mask)
        npoints = len(mask[points])  # Number of valid data points in region

        group = np.zeros(npoints)
        ed = np.zeros((npoints, nregimes[region]))

        swcf_data_pts = swcf_data[points]
        lwcf_data_pts = lwcf_data[points]

        # Assign model data to observed regimes

        for i in range(nregimes[region]):
            ed[:, i] = ((albisccp_data[points] - obs_alb[region, i]) ** 2) + \
                       ((pctisccp_data[points] - obs_pct[region, i]) ** 2) + \
                       ((cltisccp_data[points] - obs_clt[region, i]) ** 2)

        group[:] = np.argmin(ed, axis=1)

        for i in range(nregimes[region]):
            mem = (group == i)

            count = len(group[mem])

            if count > 0:

                model_rfo[region, i] = float(count) / float(npoints)
                model_ncf[region, i] = np.average(swcf_data_pts[mem]) \
                    * solar_weights[region] +                         \
                    np.average(lwcf_data_pts[mem])
            else:
                logger.info("Model does not reproduce all observed cloud "
                            "regimes.")
                logger.info("Cannot calculate CREM. Abort.")
                sys.exit()
                model_rfo[region, i] = 0.0
                model_ncf[region, i] = 0.0

    # Calculation of eq 3 in WW09
    for region in range(3):
        rCREMpd[region, 0:nregimes[region]] = area_weights[region] * \
            (((model_ncf[region, 0:nregimes[region]] -
               obs_ncf[region, 0:nregimes[region]]) *
             obs_rfo[region, 0:nregimes[region]]) ** 2 +
             ((model_rfo[region, 0:nregimes[region]] -
              obs_rfo[region, 0:nregimes[region]]) *
             obs_ncf[region, 0:nregimes[region]]) ** 2) ** 0.5

    # Calculation of eq 4 in WW09
    CREMpd = ((np.sum(rCREMpd[0, :] ** 2) + np.sum(rCREMpd[1, :] ** 2) +
               np.sum(rCREMpd[2, 0:5] ** 2)) / 20.0) ** 0.5

    # A perfect CREMpd with respect to ISCCP would be 0.0
    # An estimate of observational uncertainty (obtained by calculating CREMpd
    # wrt MODIS/ERBE) is 0.96 (i.e. models with CREMpd less than 0.96 may be
    # regarded as within observational uncertainty overall, although not
    # necessarily for every regime)'.
    # Interrogation of the rCREMpd array from this program will indicate which
    # regimes contribute most to the total CREMpd (elements ordered as Table 3
    # of WW09)'

    return CREMpd, rCREMpd

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
