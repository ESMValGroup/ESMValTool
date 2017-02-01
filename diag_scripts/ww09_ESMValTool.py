"""
;;#############################################################################
;; Cloud Regime Error Metrics (CREM)
;; Author: Keith Williams (Metoffice, UK)
;; ESA-CMUG project
;;#############################################################################
;; Description
;;    Calculates the Cloud Regime Error Metric (CREM) following Williams and
;;    Webb (2009, Clim. Dyn.)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    none
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;
;; Modification history
;;    20151117-A_laue_ax: added parameters for call to "write_references"
;;    20151113-A_laue_ax: added creation of directory for plots if needed
;;                        (code was crashing if directory does not exist)
;;    20151029-A_laue_ax: added output of acknowledgements + processed files
;;                        to log-file
;;    20150903-A_laue_ax: ESMValTool implementation.
;;    20150521-A_will_ke: CREM routines written.
;;
;;#############################################################################
"""

# Basic Python packages

import sys   # debug
import os

import ConfigParser
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.interpolation import map_coordinates as interp2d
from netCDF4 import Dataset
#from scipy.io import netcdf as nc

# ESMValTool defined Python packages
from esmval_lib import ESMValProject
from auxiliary import info


def main(project_info):

    # print(">>>>>>>> entering ww09_ESMValTool.py <<<<<<<<<<<<")

    # create instance of a wrapper that allows easy access to data
    E = ESMValProject(project_info)

    config_file = E.get_configfile()
    plot_dir = E.get_plot_dir()
    verbosity = E.get_verbosity()
    plot_type = E.get_graphic_format()
    diag_script = E.get_diag_script_name()

    res = E.write_references(diag_script,              # diag script name
                             ["A_will_ke"],            # authors
                             [""],                     # contributors
                             ["D_Williams09climdyn"],  # diag_references
                             ["E_isccp_d1"],           # obs_references
                             ["P_cmug"],               # proj_references
                             project_info,
                             verbosity,
                             False)

    modelconfig = ConfigParser.ConfigParser()
    modelconfig.read(config_file)

    # create list of model names (plot labels)
    models = []
    for model in E.project_info['MODELS']:
        models.append(model.split_entries()[1])

    # number of models
    nummod = len(E.project_info['MODELS'])
    crems = np.empty(nummod)

    # get filenames of preprocessed climatological mean files (all models)
    fn_alb = get_climo_filenames(E, variable='albisccp')
    fn_pct = get_climo_filenames(E, variable='pctisccp')
    fn_clt = get_climo_filenames(E, variable='cltisccp')
    fn_su  = get_climo_filenames(E, variable='rsut')
    fn_suc = get_climo_filenames(E, variable='rsutcs')
    fn_lu  = get_climo_filenames(E, variable='rlut')
    fn_luc = get_climo_filenames(E, variable='rlutcs')
    fn_snc = get_climo_filenames(E, variable='snc')
    fn_sic = get_climo_filenames(E, variable='sic')

    if not fn_snc:
        print("no data for variable snc found, using variable snw instead")
        fn_snw = get_climo_filenames(E, variable='snw')

    # loop over models and calulate CREM

    for i in range(nummod):

        if fn_snc:
            snc = fn_snc[i]
        else:
            snc = ""

        pointers = {'albisccp_nc': fn_alb[i],
                    'pctisccp_nc': fn_pct[i],
                    'cltisccp_nc': fn_clt[i],
                    'rsut_nc'    : fn_su[i],
                    'rsutcs_nc'  : fn_suc[i],
                    'rlut_nc'    : fn_lu[i],
                    'rlutcs_nc'  : fn_luc[i],
                    'snc_nc'     : snc,
                    'sic_nc'     : fn_sic[i]}

        if not fn_snc:
            pointers['snw_nc'] = fn_snw[i]

        # calculate CREM

        (CREMpd, __) = crem_calc(E, pointers)

        crems[i] = CREMpd

    print("------------------------------------")
    print(crems)
    print("------------------------------------")

    # plot results

    fig = plt.figure()
    ypos = np.arange(nummod)
    plt.barh(ypos, crems)
    plt.yticks(ypos + 0.5, models)
    plt.xlabel('Cloud Regime Error Metric')

    # draw observational uncertainties (dashed red line)
    plt.plot([0.96, 0.96], [0, nummod], 'r--')

    # if needed, create directory for plots
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    plt.savefig(plot_dir + 'ww09_metric_multimodel.' + plot_type)
    print("Wrote " + plot_dir + "ww09_metric_multimodel." + plot_type)


def get_climo_filenames(E, variable):

    import projects
    import os

    res = []

    for currDiag in E.project_info['DIAGNOSTICS']:
        variables = currDiag.get_variables()
        field_types = currDiag.get_field_types()
        mip = currDiag.get_var_attr_mip()
        exp = currDiag.get_var_attr_exp()
        for idx in range(len(variables)):
            for model in E.project_info['MODELS']:
                currProject = getattr(vars()['projects'],
                                      model.split_entries()[0])()
                fullpath = currProject.get_cf_fullpath(E.project_info,
                                                       model,
                                                       field_types[idx],
                                                       variables[idx],
                                                       mip[idx],
                                                       exp[idx])
                if variable == variables[idx] and os.path.isfile(fullpath):
                    res.append(fullpath)

    return res

# Reading and Regridding functions (scroll down for main program)


def regrid(aIn, xIn, yIn, xOut, yOut, fixmdis=True, xCyclic=0.0):

    """
    Function for regridding onto 2.5 degree lat-long grid as the ISCCP
    obs data used for comparison was stored on.
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
        aiNew = range(nx) + [0]
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
    """

    npts = len(lons2)
    nrows = len(lats2)

#    nt = len(Dataset(sSrcFilename, 'r', format='NETCDF3').variables['time'][:])
    nt = len(Dataset(sSrcFilename, 'r').variables['time'][:])
#    f = nc.netcdf_file(sSrcFilename, 'r')
#    time = f.variables['time'][:]
#    nt = time.shape[0]  # number of time steps in netCDF
    data_rg = np.zeros((nt, nrows, npts))
    print 'Number of data times in file ', nt

    # read data
    srcDataset = Dataset(sSrcFilename, 'r', format='NETCDF3')
    srcData = srcDataset.variables[sVarname]
#    srcData = f.variables[sVarname][:, :, :]

    # grid of input data
    lats = srcDataset.variables['lat'][:]
    lons = srcDataset.variables['lon'][:]
#    lats = f.variables['lat'][:]
#    lons = f.variables['lon'][:]

    # create mask (missing values)
#    data = np.ma.masked_equal(srcData, f.variables[sVarname]._FillValue)
    data = np.ma.masked_equal(srcData, srcData._FillValue)

    for iT in range(nt):    # range over fields in the file
        data_rg[iT, :, :] = regrid(data[iT, :, :], lons, lats, lons2, lats2,
                                   xCyclic=360.0)

    rgmasked = np.ma.masked_invalid(data_rg)
    np.ma.set_fill_value(rgmasked, 0.0)

#    return data_rg
    return(np.ma.filled(rgmasked))


def crem_calc(E, pointers):

    """
    Main program for calculating Cloud Regime Error Metric following equation
    4 in Williams and Webb (2009) (WW09) from CMOR-compliant netCDF data.

    Inputs:
    pointers - Dictionary of paths to the required netCDF files with
               the following keys:
               albisccp_nc
               pctisccp_nc
               cltisccp_nc
               rsut_nc
               rsutcs_nc
               rlut_nc
               rlutcs_nc
               snc_nc
               sic_nc

    For CMIP5, snc is in the CMIP5 table 'day'. All other variables
    are in the CMIP5 table 'cfday'. A minimum of 2 years, and ideally 5
    years, of data are required. The observational regime characteristics
    were calculated for the period Mar 1985 - Feb 1990.

    If snc is not available then snw can be used instead. In this case
    pointers[snc_nc] should be set to None and snw_nc set.

    Outputs:
    CREMpd is the present-day cloud regime error metric of WW09.
    rCREMpd is the component from each regime.
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
    obs_ncf = np.array([[-10.14, -25.45, -5.80, -27.40, -16.83, -48.45, -55.84],
                        [-13.67, -58.28, -36.26, -25.34, -64.27, -56.91, -11.63],
                        [-3.35, -16.66, -13.76, -8.63, -12.17, 1.45, 999.9]])

    area_weights = np.array([0.342, 0.502, 0.156])   # aw in eq 3 of WW09
    solar_weights = np.array([1.000, 0.998, 0.846])  # weighting for swcf to
                                                     # account for lack of ISCCP
                                                     # diagnostics during polar
                                                     # night (p153 of WW09)

    nregimes = np.array([7, 7, 6])  # number of regimes in each region
                                    # (Table 3 of WW09)

    #-----------------------------------------------------------

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

    print 'Reading and regridding albisccp_nc'
    albisccp_data = read_and_regrid(pointers['albisccp_nc'], sVarnames[0], lons2, lats2)
    E.add_to_filelist(pointers['albisccp_nc'])
    print 'Reading and regridding pctisccp_nc'
    pctisccp_data = read_and_regrid(pointers['pctisccp_nc'], sVarnames[1], lons2, lats2)
    E.add_to_filelist(pointers['pctisccp_nc'])
    print 'Reading and regridding cltisccp_nc'
    cltisccp_data = read_and_regrid(pointers['cltisccp_nc'], sVarnames[2], lons2, lats2)
    E.add_to_filelist(pointers['cltisccp_nc'])
    print 'Reading and regridding rsut_nc'
    rsut_data = read_and_regrid(pointers['rsut_nc'], sVarnames[3], lons2, lats2)
    E.add_to_filelist(pointers['rsut_nc'])
    print 'Reading and regridding rsutcs_nc'
    rsutcs_data = read_and_regrid(pointers['rsutcs_nc'], sVarnames[4], lons2, lats2)
    E.add_to_filelist(pointers['rsutcs_nc'])
    print 'Reading and regridding rlut_nc'
    rlut_data = read_and_regrid(pointers['rlut_nc'], sVarnames[5], lons2, lats2)
    E.add_to_filelist(pointers['rlut_nc'])
    print 'Reading and regridding rlutcs_nc'
    rlutcs_data = read_and_regrid(pointers['rlutcs_nc'], sVarnames[6], lons2, lats2)
    E.add_to_filelist(pointers['rlutcs_nc'])
    print 'Reading and regridding sic_nc'
    sic_data = read_and_regrid(pointers['sic_nc'], sVarnames[8], lons2, lats2)
    E.add_to_filelist(pointers['sic_nc'])
    if not pointers['snc_nc']:
        print 'Reading and regridding snw_nc'
        snc_data = read_and_regrid(pointers['snw_nc'], sVarnames[7], lons2, lats2)
        E.add_to_filelist(pointers['snw_nc'])
    else:
        print 'Reading and regridding snc_nc'
        snc_data = read_and_regrid(pointers['snc_nc'], sVarnames[7], lons2, lats2)
        E.add_to_filelist(pointers['snc_nc'])

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

    print 'Assigning data to observational cloud regimes'
    for region in range(3):  # over 3 regions (tropics, extra tropics, snow/ice)

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
                                       * solar_weights[region] +      \
                                       np.average(lwcf_data_pts[mem])
            else:
                print("Model does not reproduce all observed cloud regimes.")
                print("Cannot calculate CREM. Abort.")
                sys.exit()
                model_rfo[region, i] = 0.0
                model_ncf[region, i] = 0.0

    # Calculation of eq 3 in WW09
    for region in range(3):
        rCREMpd[region, 0:nregimes[region]] = area_weights[region] * \
            (((model_ncf[region, 0:nregimes[region]] - obs_ncf[region, 0:nregimes[region]])
            * obs_rfo[region, 0:nregimes[region]]) ** 2 +
            ((model_rfo[region, 0:nregimes[region]] - obs_rfo[region, 0:nregimes[region]])
            * obs_ncf[region, 0:nregimes[region]]) ** 2) ** 0.5

    # Calculation of eq 4 in WW09
    CREMpd = ((np.sum(rCREMpd[0, :] ** 2) + np.sum(rCREMpd[1, :] ** 2)
               + np.sum(rCREMpd[2, 0:5] ** 2)) / 20.0) ** 0.5

    # A perfect CREMpd with respect to ISCCP would be 0.0
    # An estimate of observational uncertainty (obtained by calculating CREMpd
    # wrt MODIS/ERBE) is 0.96 (i.e. models with CREMpd less than 0.96 may be
    # regarded as within observational uncertainty overall, although not
    # necessarily for every regime)'.
    # Interrogation of the rCREMpd array from this program will indicate which
    # regimes contribute most to the total CREMpd (elements ordered as Table 3
    # of WW09)'

    return CREMpd, rCREMpd
