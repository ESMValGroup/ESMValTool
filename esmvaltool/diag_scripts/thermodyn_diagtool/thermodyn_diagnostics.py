# pylint: disable-msg=R0801
"""Tool for computation of some aspects of climate system thermodynamics.

Author
Valerio Lembo
(Meteorological Institute, Hamburg University - valerio.lembo@uni-hamburg.de)

Contributors
Frank   Lunkeit
(Meteorological Insitute, Hamburg University - f.lunkeit@uni-hamburg.de)
Nikolay Koldunov
(MARUM/AWI, nikolay.koldunov@awi.de, Germany)

Project
CRC - TRR 181 "Energy transfers in Atmosphere and Ocean"

#############################################################################

PREREQUISITES

The program shares the same prerequisites with the overall ESMValTool
architecture (see http://esmvaltool.readthedocs.io/en/latest/install.html)

ADDITIONAL REQUIREMENTS
- NCO: NCO command operators are required for attributes and coordinates
       manipulation on NetCDF files. The programme has been tested with
       the following  versions of NCO precompiled libraries:
       * Unix:   nco-4.6.7
       * Mac-OS: nco-4.6.6
- PyAstronomy: Python package PyAstronomy contains a function "zeocross1d"
               which is needed for the computation of transport maxima
               locations. We tested the method with v. 0.12.0 of the
               package, available in pip repositories.

USAGE

1: Obtain the datasets: the programme accepts the following variables as
   input for the computations:
     Monthly mean resolution or higher:
     - TOA shortwave radiation downwards;
     - TOA shortwave radiation upwards;
     - TOA longwave radiation upwards (OLR);
     - Surface shortwave radiation downwards;
     - Surface shortwave radiation upwards;
     - Surface longwave radiation downwards;
     - Surface longwave radiation upwards;
     - Surface turbulent latent heat fluxes;
     - Surface turbulent sensible heat fluxes;
     - Surface temperature;
     - Specific humidity;
     - Near-surface (or 10m) zonal velocity;
     - Near-surface (or 10m) meridional velocity;
     Daily mean resolution or higher:
     - Near-surface temperature;
     - Air temperature (on pressure levels);
     - Horizontal velocity (on pressure levels);
     - Meridional velocity (on pressure levels);
     - Vertical velocity (on pressure levels);
   Data on lonlat grid are accepted, with CMOR-compliant coordinate system.
   The pre-processing modules of ESMValTool scheme will take care of
   converting known grids and recognized datasets to CMOR standards. For a
   a list of known formats, see
      http://esmvaltool.readthedocs.io/en/latest/running.html#tab-obs-dat
   Whether you need to use non-recognized datasets (provided that they are
   globally gridded), ask for support!

2: A configuration template is available in the ESMValTool release. Set your
  own paths to local directories here. Input datasets are read in MODELPATH,
  MODELPATH2, OBSPATH or OBSPATH2 output datasets are stored in WORKPATH,
  plots in PLOTPATH (refer to the manual for ESMValTool).

3: Set the namelist with your own input datasets. Here you can also set the
   length of the dataset you want to subset, the time resolution, and
   whether you are working on a MacOS or Linux machine.
4: Run the tool by typing:
         python main.py nml/namelist_lembo17.xml
5: After the data preprocessing, choose the options as required by the
   program;

OUTPUT

In the output directory, NetCDF files containing annual mean fields and
global mean time series for budgets and material entropy production time
series.
Log files containing the values of the LEC components are stored in
subdirectories, identified by the year for which they are computed. LEC
quantities are stored in different files, one for each year.
- energy_gns_* contains global, NH and SH time mean fields as a function of
  wavenumbers;
- energy_mean_* contains time mean (lon,wave,plev) fields;
- energy_ts3d_* contains (time,lon,wave,plev) fields;
- <MODELNAME>_energy_*.ps is a diagram flux displaying annual mean
  components of the LEC;

In the plot directory figures are stored as PNG files in each model subdir.

For the budgets (toab: TOA energy budget, atmb: atmospheric energy budget,
surb: surface energy budget, latent: latent energy budget, wmass: water mass
budget) there can be found:
 - climatological mean fields (*_climap.png);
 - time series evolution of global, NH and SH mean fields (*_timeser.png)
 - meridional transports (*_transp.pdf);
 - scatter plots of annual peak magnitudes vs. locations (*_scatpeak.png)
 - LEC diagrams for each year (*_lec_diagram.png)

For the multi-model ensembles scatter plots of EB mean values vs. their
interannual variability are provided (scatters_variability.png), a summary
panel of the main thermodynamic quantities averaged for each model
(scatters_summary.png), and the zonal mean meridional heat transports for
each model (meridional_transp.png).

For the material entropy production (sver: vertical, shor: horizontal
(indirect method), sevap: evaporation component, slatpr: ice-rain
phase changes, slatps: vapor-snow phase changes, spotp: potential energy,
srain: rainfall, ssens: sensible heat, ssnow: snowfall) climatological mean
fields are shown.

In the working directory, a log file (lembo17_log.txt) is produced,
containing the global mean values of all the quantities.

N.B.: multi-model ensembles and means are performed on the results of the
      analysis. No multi-model mean analysis is performed on input fields.
      It is thus a deliberate choice not to use the multi-model mean tools
      provided by the postprocessor.

SOFTWARE TREE DESCRIPTION

The tool is divided in three modules; one for the
computation of energy and water mass budgets (and related meridional
transports), one for the Lorenz Energy Cycle (LEC), one for the material
entropy production.

Module Dependencies:
  - MODULE 1: Budgets and Transports --> none;
  - MODULE 2: Lorenz Energy Cycle    --> none;
  - MODULE 3: Material Entropy Production:
                                 * Indirect method --> requires Budgets
                                 * Direct method   --> requires Budgets and
                                                       LEC

- MODULE 1
Earth's energy budgets from radiative and heat fluxes at Top-of-Atmosphere,
at the surface and in the atmosphere (as a residual).
If required by the user, water mass and latent energy budgets are computed
from rainfall, snowfall and evaporation fluxes.
If required by the user, computations are also separately performed over
land and oceans. A land-sea mask (either binary or percentual) has to be
provided.
Meridional transports, magnitude and location of the peaks in each
hemisphere (only for heat transports) are also computed.
The baroclinic efficiency is computed from TOA energy budgets, emission
temperature (in turn retrieved from OLR) and near-surface temperature.

- MODULE 2

The Lorenz Energy Cycle (LEC) is computed in spectral components from near-
surface temperatures, temperatures and the three components of velocities
over pressure levels.
The storage and conversion terms are directly computed, the sources and
sinks are retrieved as residuals.
Components are grouped into a zonal mean, stationary and transient eddy
part.

- MODULE 3

The material entropy production is computed using the indirect method, the
direct method or both (following Lucarini et al., 2014).
For the indirect method a vertical and a horizontal component are provided.
For the direct method, all components are combined, related to the
hydrological cycle (attributable to evaporation, rainfall and snowfall
precipitation, phase changes and potential energy of the droplet), to the
sensible heat fluxes and to kinetic energy dissipation. For the latter the
LEC computation is required, given that the strength of the LEC can be
considered as equal to the kinetic energy dissipated to heating.




#############################################################################

20170803-A_lemb_va: Modified header with description and caveats
20170629-A_kold_ni: Atmospheric budgets diagnostics written
20180524-A_lemb_va: first complete working thermodynamics diagnostics


#############################################################################
"""

# pylint: disable-msg=C0412
import os
from shutil import move
import warnings
# New packages for version 2.0 of ESMValTool
import logging
import esmvaltool.diag_scripts.shared as e
from cdo import Cdo
from netCDF4 import Dataset
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
# Locally used modules
from esmvaltool.diag_scripts.thermodyn_diagtool import mkthe,\
                                                       fourier_coefficients,\
                                                       lorenz_cycle,\
                                                       plot_script
matplotlib.use('Agg')
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
logger = logging.getLogger(os.path.basename(__file__))

SIGMAINV = 17636684.3034 	# inverse of the Stefan-Boltzmann constant
L_C = 2501000 		        # latent heat of condensation
LC_SUB = 2835000 		    # latent heat of sublimation
L_S = 334000		            # latent heat of solidification
GRAV = 9.81		            # gravity acceleration


# pylint: disable=C0302
# pylint: disable-msg=C0302
# pylint: disable-msg=R0914
# Two hundreds and fiftyfour is reasonable in this case.
# pylint: disable-msg=R0915
# Nine hundreds and twentyone is reasonable in this case.
# pylint: disable=too-many-branches
# Sixtytwo is reasonable in this case.
# pylint: disable-msg=unused-argument
# flake8: noqa
def main(cfg):
    """Execute the program.

    Argument cfg, containing directory paths, preprocessed input dataset
    filenames and user-defined options, is passed by ESMValTool preprocessor.
    """
    logger.info('Entering the diagnostic tool')
    # Load paths
    wdir_up = cfg['work_dir']
    pdir_up = cfg['plot_dir']
    logger.info('Work directory: %s \n', wdir_up)
    logger.info('Plot directory: %s \n', pdir_up)
    plotsmod = plot_script
    data = e.Datasets(cfg)
    logger.debug(data)
    models = data.get_info_list('dataset')
    model_names = list(set(models))
    model_names.sort()
    logger.info(model_names)
    varnames = data.get_info_list('short_name')
    curr_vars = list(set(varnames))
    logger.debug(curr_vars)
    # load user-defined options
    lsm = str(cfg['lsm'])
    wat = str(cfg['wat'])
    lec = str(cfg['lec'])
    entr = str(cfg['entr'])
    met = str(cfg['met'])
    flags = [wat, entr, met]
    # Initialize multi-model arrays
    modnum = len(model_names)
    te_all = np.zeros(modnum)
    toab_all = np.zeros([modnum, 2])
    toab_oc_all = np.zeros(modnum)
    toab_la_all = np.zeros(modnum)
    atmb_all = np.zeros([modnum, 2])
    atmb_oc_all = np.zeros(modnum)
    atmb_la_all = np.zeros(modnum)
    surb_all = np.zeros([modnum, 2])
    surb_oc_all = np.zeros(modnum)
    surb_la_all = np.zeros(modnum)
    wmb_all = np.zeros([modnum, 2])
    wmb_oc_all = np.zeros(modnum)
    wmb_la_all = np.zeros(modnum)
    latent_all = np.zeros([modnum, 2])
    latent_oc_all = np.zeros(modnum)
    latent_la_all = np.zeros(modnum)
    baroc_eff_all = np.zeros(modnum)
    lec_all = np.zeros([modnum, 2])
    horzentr_all = np.zeros([modnum, 2])
    vertentr_all = np.zeros([modnum, 2])
    matentr_all = np.zeros([modnum, 2])
    irrevers_all = np.zeros(modnum)
    diffentr_all = np.zeros([modnum, 2])
    logger.info("Entering main loop\n")
    i_m = 0
    for model in model_names:
        # Load paths to individual models output and plotting directories
        wdir = os.path.join(wdir_up, model)
        pdir = os.path.join(pdir_up, model)
        if not os.path.exists(wdir):
            os.makedirs(wdir)
        # Reading file names for the specific model
        filenames = data.get_info_list('filename', dataset=model)
        logger.info('Processing model: %s \n', model)
        hfls_file = filenames[0]
        hfss_file = filenames[1]
        prsn_file = filenames[4]
        rlds_file = filenames[6]
        rlus_file = filenames[7]
        rsds_file = filenames[9]
        rsus_file = filenames[11]
        ta_file = filenames[13]
        ts_file = filenames[15]
        # Read path to land-sea mask
        head = Dataset(ta_file)
        info = getattr(head, 'metadata')
        attr = info.split()
        diction = {attr[i_m]: attr[i_m + 1] for i_m in range(len(attr) - 1)}
        sftlf_fx = str(diction['{sftlf:'])
        sftlf_fx = sftlf_fx.replace("}", "")
        aux_file = wdir + '/aux.nc'
        te_ymm_file, te_gmean_constant, _, _ = auxiliary(model, wdir,
                                                         filenames, flags)
        te_all[i_m] = te_gmean_constant
        # Compute energy budgets
        logger.info('Computing energy budgets\n')
        # Compute energy budgets
        eb_gmean, eb_file, toab_ymm_file = comp_budgets(model, wdir, aux_file,
                                                        filenames)
        toab_all[i_m, 0] = np.nanmean(eb_gmean[0])
        toab_all[i_m, 1] = np.nanstd(eb_gmean[0])
        atmb_all[i_m, 0] = np.nanmean(eb_gmean[2])
        atmb_all[i_m, 1] = np.nanstd(eb_gmean[2])
        surb_all[i_m, 0] = np.nanmean(eb_gmean[1])
        surb_all[i_m, 1] = np.nanstd(eb_gmean[1])
        logger.info('Atmospheric energy budget: %s\n', atmb_all[i_m, 0])
        logger.info('Done\n')
        # Baroclinic efficiency
        baroc_eff_all[i_m] = comp_baroceff(model, wdir, aux_file,
                                           toab_ymm_file, te_ymm_file)
        logger.info('Running the plotting module for the budgets\n')
        plotsmod.balances(wdir_up, pdir, [eb_file[0], eb_file[1], eb_file[2]],
                          ['toab', 'atmb', 'surb'], model)
        oname = '{}/{}_{}_timeser.png'.format(pdir, model, 'toab')
        logger.info('Done\n')
        # Water mass budget
        if wat in {'y', 'yes'}:
            logger.info('Computing water mass and latent energy budgets\n')
            wm_gmean, wm_file = comp_wmbudg(model, wdir, aux_file, filenames,
                                            flags)
            wmb_all[i_m, 0] = np.nanmean(wm_gmean[0])
            wmb_all[i_m, 1] = np.nanstd(wm_gmean[0])
            logger.info('Water mass budget: %s\n', wmb_all[i_m, 0])
            latent_all[i_m, 0] = np.nanmean(wm_gmean[1])
            latent_all[i_m, 1] = np.nanstd(wm_gmean[1])
            logger.info('Latent energy budget: %s\n', latent_all[i_m, 0])
            logger.info('Done\n')
            logger.info('Plotting the water mass and latent energy budgets\n')
            plotsmod.balances(wdir_up, pdir, [wm_file[0], wm_file[1]],
                              ['wmass', 'latent'], model)
            logger.info('Done\n')
        else:
            pass
        # Compute budgets over oceans and land separately.
        if lsm in {'y', 'yes'}:
            logger.info('Computing energy budgets over land and oceans\n')
            toab_oc_gmean, toab_la_gmean = comp_landoc_budg(model, wdir,
                                                            eb_file[0],
                                                            sftlf_fx, 'toab')
            toab_oc_all[i_m] = toab_oc_gmean
            toab_la_all[i_m] = toab_la_gmean
            logger.info('TOA energy budget over oceans: %s\n',
                        toab_oc_gmean)
            logger.info('TOA energy budget over land: %s\n', toab_la_gmean)
            atmb_oc_gmean, atmb_la_gmean = comp_landoc_budg(model, wdir,
                                                            eb_file[1],
                                                            sftlf_fx, 'atmb')
            atmb_oc_all[i_m] = atmb_oc_gmean
            atmb_la_all[i_m] = atmb_la_gmean
            logger.info('Atmospheric energy budget over oceans: %s\n',
                        atmb_oc_gmean)
            logger.info('Atmospheric energy budget over land: %s\n',
                        atmb_la_gmean)
            surb_oc_gmean, surb_la_gmean = comp_landoc_budg(model, wdir,
                                                            eb_file[2],
                                                            sftlf_fx, 'surb')
            surb_oc_all[i_m] = surb_oc_gmean
            surb_la_all[i_m] = surb_la_gmean
            logger.info('Surface energy budget over oceans: %s\n',
                        surb_oc_gmean)
            logger.info('Surface energy budget over land: %s\n',
                        surb_la_gmean)
            logger.info('Done\n')
            if wat in {'y', 'yes'}:
                logger.info('Computing water mass and latent energy'
                            ' budgets over land and oceans\n')
                wmb_oc_gmean, wmb_la_gmean = comp_landoc_budg(model, wdir,
                                                              wm_file[0],
                                                              sftlf_fx, 'wmb')
                wmb_oc_all[i_m] = wmb_oc_gmean
                wmb_la_all[i_m] = wmb_la_gmean
                logger.info('Water mass budget over oceans: %s\n',
                            wmb_oc_gmean)
                logger.info('Water mass budget over land: %s\n',
                            wmb_la_gmean)
                latent_oc_gmean, latent_la_gmean = comp_landoc_budg(
                    model, wdir, wm_file[0], sftlf_fx, 'latent')
                latent_oc_all[i_m] = latent_oc_gmean
                latent_la_all[i_m] = latent_la_gmean
                logger.info('Water mass budget over oceans: %s\n',
                            latent_oc_gmean)
                logger.info('Water mass budget over land: %s\n',
                            latent_la_gmean)
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        # Compute the Lorenz Energy Cycle.
        if lec in {'y', 'yes'}:
            logger.info('Computation of the Lorenz Energy'
                        'Cycle (year by year)\n')
            lect = lec_preproc(model, wdir, pdir, filenames)
            lec_all[i_m, 0] = np.nanmean(lect)
            lec_all[i_m, 1] = np.nanstd(lect)
            logger.info('Intensity of the annual mean Lorenz Energy '
                        'Cycle: %s\n', lec_all[i_m, 0])
            logger.info('Done\n')
        else:
            pass
        # Compute the material entropy production
        if entr in {'y', 'yes'}:
            if met in {'1', '3'}:
                _, _, te_file, _ = auxiliary(model, wdir, filenames, flags)
                logger.info('Computation of the material entropy production '
                            'with the indirect method\n')
                indentr_list = [rlds_file, rlus_file, rsds_file, rsus_file,
                                te_file, eb_file[0], ts_file]
                horzentr_mean, vertentr_mean, vertentr_file = comp_indentr(
                    model, wdir, indentr_list,
                    aux_file, eb_gmean[0])
                horzentr_all[i_m, 0] = np.nanmean(horzentr_mean)
                horzentr_all[i_m, 1] = np.nanstd(horzentr_mean)
                vertentr_all[i_m, 0] = np.nanmean(vertentr_mean)
                vertentr_all[i_m, 1] = np.nanstd(vertentr_mean)
                logger.info('Horizontal component of the material entropy '
                            'production: %s\n', horzentr_all[i_m, 0])
                logger.info('Vertical component of the material entropy '
                            'production: %s\n', vertentr_all[i_m, 0])
                logger.info('Done\n')
                logger.info('Running the plotting module for the material '
                            'entropy production (indirect method)\n')
                plotsmod.entropy(pdir, vertentr_file, 'sver',
                                 'Vertical entropy production', model)
                logger.info('Done\n')
            if met in {'2', '3'}:
                _, _, _, aux_files = auxiliary(model, wdir, filenames, flags)
                htop_file = aux_files[1]
                prr_file = aux_files[2]
                tabl_file = aux_files[3]
                tasvert_file = aux_files[4]
                tcloud_file = aux_files[5]
                tcolumn_file = aux_files[6]
                tlcl_file = aux_files[7]
                logger.info('Computation of the material entropy '
                            'production with the direct method\n')
                logger.info('1. Sensible heat fluxes\n')
                infile_list = [hfss_file, tabl_file, ts_file]
                sensentr_mean, sensentr_file = comp_sensentr(
                    model, wdir, infile_list, aux_file)
                logger.info('Material entropy production associated with '
                            'sens. heat fluxes: %s\n', sensentr_mean)
                logger.info('Done\n')
                logger.info('2. Hydrological cycle\n')
                logger.info('2.1 Evaporation fluxes\n')
                infile_list = [hfls_file, ts_file]
                evapentr_mean, evapentr_file = comp_evapentr(model, wdir,
                                                             infile_list,
                                                             aux_file)
                logger.info('Material entropy production associated with '
                            'evaporation fluxes: %s\n', evapentr_mean)
                logger.info('Done\n')
                infile_mask = [prr_file, prsn_file, tlcl_file]
                prrmask_file, prsnmask_file = mask_precip(model, wdir,
                                                          infile_mask)
                logger.info('2.2 Rainfall precipitation\n')
                infile_rain = [prrmask_file, tcloud_file]
                rainentr_mean, rainentr_file = comp_rainentr(model, wdir,
                                                             infile_rain,
                                                             aux_file)
                logger.info('Material entropy production associated with '
                            'rainfall: %s\n', rainentr_mean)
                logger.info('Done\n')
                logger.info('2.3 Snowfall precipitation\n')
                infile_rain = [prsnmask_file, tcloud_file]
                snowentr_mean, latsnow_file, snowentr_file = comp_snowentr(
                    model, wdir, infile_rain, aux_file)
                logger.info('Material entropy production associated with '
                            'snowfall: %s\n', snowentr_mean)
                logger.info('Done\n')
                logger.info('2.4 Melting of snow at the surface \n')
                meltentr_mean, meltentr_file = comp_meltentr(
                    model, wdir, latsnow_file, aux_file)
                logger.info('Material entropy production associated with snow '
                            'melting: %s\n', meltentr_mean)
                logger.info('Done\n')
                logger.info('2.5 Potential energy of the droplet\n')
                infile_pot = [htop_file, prrmask_file, prsnmask_file,
                              tcolumn_file]
                potentr_mean, potentr_file = comp_potentr(
                    model, wdir, infile_pot, aux_file)
                logger.info('Material entropy production associated with '
                            'potential energy of the droplet: %s\n',
                            potentr_mean)
                logger.info('Done\n')
                logger.info('3. Kinetic energy dissipation\n')
                minentr_mean = comp_kinentr(aux_file, tasvert_file, lect, lec)
                matentr = (float(sensentr_mean) - float(evapentr_mean)
                           + float(rainentr_mean) + float(snowentr_mean)
                           + float(potentr_mean) + float(minentr_mean)
                           - float(meltentr_mean))
                logger.info('Material entropy production with '
                            'the direct method: %s\n', matentr)
                matentr_all[i_m, 0] = matentr
                if met in {'3'}:
                    diffentr = (float(np.nanmean(vertentr_mean))
                                + float(np.nanmean(horzentr_mean))
                                - matentr)
                    logger.info('Difference between the two '
                                'methods: %s\n', diffentr)
                    diffentr_all[i_m, 0] = diffentr
                else:
                    pass
                irrevers = ((matentr - float(minentr_mean)) /
                            float(minentr_mean))
                logger.info('Degree of irreversibility of the '
                            'system: %s\n', irrevers)
                irrevers_all[i_m] = irrevers
                logger.info('Running the plotting module for the material '
                            'entropy production (direct method)\n')
                entr_list = [sensentr_file, evapentr_file, rainentr_file,
                             snowentr_file, meltentr_file, potentr_file]
                plotsmod.init_plotentr(model, pdir, entr_list)
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        logger.info('Done for model: %s \n', model)
        i_m = i_m + 1
    # Produce multi-model ensemble plots (if more than one model are taken
    # into account).
    logger.info('I will now start multi-model plots')
    logger.info('Meridional heat transports\n')
    plotsmod.plot_mm_transp(model_names, wdir_up, pdir_up)

    logger.info('Scatter plots')
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    axi = plt.subplot(321)
    title = '(a) TOA vs. atmospheric energy budget'
    xlabel = 'R_t [W m-2]'
    ylabel = 'F_a [W m-2]'
    varlist = [toab_all[:, 0], atmb_all[:, 0]]
    plotsmod.plot_mm_scatter(axi, varlist, title, xlabel, ylabel)
    
    axi = plt.subplot(322)
    title = '(b) Baroclinic efficiency vs. Intensity of LEC'
    xlabel = 'Eta'
    ylabel = 'W [W/m2]'
    varlist = [baroc_eff_all, lec_all[:, 0]]
    plotsmod.plot_mm_scatter(axi, varlist, title, xlabel, ylabel)

    axi = plt.subplot(323)
    title = '(c) Vertical vs. horizontal component'
    xlabel = 'S_hor [W m-2 K-1]'
    ylabel = 'S_ver [W m-2 K-1]'
    varlist = [horzentr_all[:, 0], vertentr_all[:, 0]]
    plotsmod.plot_mm_scatter_spec(axi, varlist, title, xlabel, ylabel)

    indentr_all = horzentr_all[:, 0] + vertentr_all[:, 0]
    axi = plt.subplot(324)
    axi.set_figsize = (50, 50)
    plt.scatter(indentr_all, matentr_all[:, 0], c=colors, alpha=1)
    plt.scatter(np.nanmean(indentr_all), np.nanmean(matentr_all[:, 0]),
                c='red')
    s_l, _, _, _, _ = stats.linregress(indentr_all, matentr_all[:, 0])
    plotsmod.plot_ellipse(semimin=np.nanstd(indentr_all),
                          semimaj=np.nanstd(matentr_all[:, 0]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(indentr_all),
                          y_cent=np.nanmean(matentr_all[:, 0]),
                          a_x=axi)
    plt.title('(d) Indirect vs. direct method', fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('S_ind [W m-2 K-1]', fontsize=14)
    plt.ylabel('S_dir [W m-2 K-1]', fontsize=14)
    d_x = 0.01 * (max(indentr_all) - min(indentr_all))
    d_y = 0.01 * (max(matentr_all[:, 0]) - min(matentr_all[:, 0]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (indentr_all[i_m], matentr_all[i_m, 0]),
                     xytext=(indentr_all[i_m] + d_x,
                             matentr_all[i_m, 0] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(325)
    axi.set_figsize = (50, 50)
    plt.scatter(te_all, indentr_all, c=colors, alpha=1)
    plt.scatter(np.nanmean(te_all), np.nanmean(indentr_all), c='red')
    s_l, _, _, _, _ = stats.linregress(te_all, indentr_all)
    plotsmod.plot_ellipse(semimaj=np.nanstd(te_all),
                          semimin=np.nanstd(indentr_all),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(te_all),
                          y_cent=np.nanmean(indentr_all),
                          a_x=axi)
    plt.title('(e) Indirect method vs. emission temperature', fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('T_E [K]', fontsize=14)
    plt.ylabel('S_mat [W m-2 K-1]', fontsize=14)
    d_x = 0.01 * (max(te_all) - min(te_all))
    d_y = 0.01 * (max(indentr_all) - min(indentr_all))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (te_all[i_m], indentr_all[i_m]),
                     xytext=(te_all[i_m] + d_x, (indentr_all[i_m]) + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(326)
    axi.set_figsize = (50, 50)
    plt.scatter(te_all, baroc_eff_all, c=colors, alpha=1)
    plt.scatter(np.nanmean(te_all), np.nanmean(baroc_eff_all), c='red')
    s_l, _, _, _, _ = stats.linregress(te_all, baroc_eff_all)
    plotsmod.plot_ellipse(semimaj=np.nanstd(te_all),
                          semimin=np.nanstd(baroc_eff_all),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(te_all),
                          y_cent=np.nanmean(baroc_eff_all),
                          a_x=axi)
    plt.title('(f) Baroclinic efficiency vs. emission temperature',
              fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('T_E [K]', fontsize=14)
    plt.ylabel('Eta', fontsize=14)
    d_x = 0.01 * (max(te_all) - min(te_all))
    d_y = 0.01 * (max(baroc_eff_all) - min(baroc_eff_all))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (te_all[i_m], baroc_eff_all[i_m]),
                     xytext=(te_all[i_m] + d_x, baroc_eff_all[i_m] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.grid()
    oname = pdir_up + '/scatters_summary.png'
    plt.savefig(oname)
    plt.subplots_adjust(hspace=.3)
    # Scatter plot of climatological mean values vs. interannual
    # variability for each model
    logger.info('Scatter plots for inter-annual variability of '
                'some quantities')
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    colors = (0, 0, 0)
    axi = plt.subplot(221)
    axi.set_figsize = (50, 50)
    plt.scatter(toab_all[:, 0], toab_all[:, 1], c=colors, alpha=1)
    plt.scatter(np.nanmean(toab_all[:, 0]), np.nanmean(toab_all[:, 1]),
                c='red')
    s_l, _, _, _, _ = stats.linregress(toab_all[:, 0], toab_all[:, 1])
    plotsmod.plot_ellipse(semimaj=np.nanstd(toab_all[:, 0]),
                          semimin=np.nanstd(toab_all[:, 1]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(toab_all[:, 0]),
                          y_cent=np.nanmean(toab_all[:, 1]), a_x=axi)
    plt.title('(a) TOA energy budget', fontsize=14)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('R_t [W m-2]', fontsize=14)
    plt.ylabel('Sigma (R_t) [W m-2]', fontsize=14)
    d_x = 0.01 * (max(toab_all[:, 0]) - min(toab_all[:, 0]))
    d_y = 0.01 * (max(toab_all[:, 1]) - min(toab_all[:, 1]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (toab_all[i_m, 0], toab_all[i_m, 1]),
                     xytext=(toab_all[i_m, 0] + d_x, toab_all[i_m, 1] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.ylim(bottom=0)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(222)
    axi.set_figsize = (50, 50)
    plt.scatter(atmb_all[:, 0], atmb_all[:, 1], c=colors, alpha=1)
    plt.scatter(np.nanmean(atmb_all[:, 0]), np.nanmean(atmb_all[:, 1]),
                c='red')
    s_l, _, _, _, _ = stats.linregress(atmb_all[:, 0], atmb_all[:, 1])
    plotsmod.plot_ellipse(semimaj=np.nanstd(atmb_all[:, 0]),
                          semimin=np.nanstd(atmb_all[:, 1]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(atmb_all[:, 0]),
                          y_cent=np.nanmean(atmb_all[:, 1]), a_x=axi)
    plt.title('(b) Atmospheric energy budget', fontsize=14)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('F_a [W m-2]', fontsize=14)
    plt.ylabel('Sigma (F_a) [W m-2]', fontsize=14)
    d_x = 0.01 * (max(atmb_all[:, 0]) - min(atmb_all[:, 0]))
    d_y = 0.01 * (max(atmb_all[:, 1]) - min(atmb_all[:, 1]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (atmb_all[i_m, 0], atmb_all[i_m, 1]),
                     xytext=(atmb_all[i_m, 0] + d_x, atmb_all[i_m, 1] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.ylim(bottom=0)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(223)
    axi.set_figsize = (50, 50)
    plt.scatter(surb_all[:, 0], surb_all[:, 1], c=colors, alpha=1)
    plt.scatter(np.nanmean(surb_all[:, 0]), np.nanmean(surb_all[:, 1]),
                c='red')
    s_l, _, _, _, _ = stats.linregress(surb_all[:, 0], surb_all[:, 1])
    plotsmod.plot_ellipse(semimaj=np.nanstd(surb_all[:, 0]),
                          semimin=np.nanstd(surb_all[:, 1]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(surb_all[:, 0]),
                          y_cent=np.nanmean(surb_all[:, 1]), a_x=axi)
    plt.title('(c) Surface energy budget', fontsize=14)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('F_s [W m-2]', fontsize=14)
    plt.ylabel('Sigma (F_s) [W m-2]', fontsize=14)
    d_x = 0.01 * (max(surb_all[:, 0]) - min(surb_all[:, 0]))
    d_y = 0.01 * (max(surb_all[:, 1]) - min(surb_all[:, 1]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (surb_all[i_m, 0], surb_all[i_m, 1]),
                     xytext=(surb_all[i_m, 0] + d_x, surb_all[i_m, 1] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.ylim(bottom=0)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(224)
    axi.set_figsize = (50, 50)
    plt.errorbar(x=atmb_all[:, 0], y=surb_all[:, 0],
                 xerr=atmb_all[:, 1], yerr=surb_all[:, 1],
                 fmt='none', ecolor=colors)
    plt.scatter(np.nanmean(atmb_all[:, 0]), np.nanmean(surb_all[:, 0]),
                c='red')
    s_l, _, _, _, _ = stats.linregress(surb_all[:, 0], surb_all[:, 1])
    plotsmod.plot_ellipse(semimaj=np.nanstd(atmb_all[:, 0]),
                          semimin=np.nanstd(surb_all[:, 0]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(atmb_all[:, 0]),
                          y_cent=np.nanmean(surb_all[:, 0]), a_x=axi)
    plt.title('(d) Atmospheric vs. Surface budget', fontsize=14)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('F_a [W m-2]', fontsize=14)
    plt.ylabel('F_s [W m-2]', fontsize=14)
    d_x = 0.01 * (max(atmb_all[:, 0]) - min(atmb_all[:, 0]))
    d_y = 0.01 * (max(surb_all[:, 0]) - min(surb_all[:, 0]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (atmb_all[i_m, 0], surb_all[i_m, 0]),
                     xytext=(atmb_all[i_m, 0] + d_x, surb_all[i_m, 0] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    plt.savefig(pdir_up + '/scatters_variability.png')
    plt.close(fig)
    logger.info("The diagnostic has finished. Now closing...\n")


def auxiliary(model, wdir, filelist, flags):
    """Compute auxiliary fields or perform time averaging of existing fields.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - filelist: a list of file names containing the input fields;
    - flags: (wat: a flag for the water mass budget module (y or n),
              entr: a flag for the material entropy production (y or n);
              met: a flag for the material entropy production method
              (1: indirect, 2, direct, 3: both));

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo
    mkth = mkthe
    wat = flags[0]
    entr = flags[1]
    met = flags[2]
    hfss_file = filelist[1]
    hus_file = filelist[2]
    ps_file = filelist[5]
    rlut_file = filelist[8]
    tas_file = filelist[14]
    ts_file = filelist[15]
    uas_file = filelist[17]
    vas_file = filelist[19]
    # Compute monthly mean fields from 2D surface daily fields
    aux_file = wdir + '/aux.nc'
    cdo.selvar('tas', input=tas_file, output=aux_file)
    move(aux_file, tas_file)
    tasmn_file = wdir + '/{}_tas_mm.nc'.format(model)
    cdo.selvar('tas', input='-monmean {}'.format(tas_file),
               option='-b F32', output=tasmn_file)
    cdo.selvar('uas', input=uas_file, output=aux_file)
    move(aux_file, uas_file)
    uasmn_file = wdir + '/{}_uas_mm.nc'.format(model)
    cdo.selvar('uas', input='-monmean {}'.format(uas_file),
               option='-b F32', output=uasmn_file)
    cdo.selvar('vas', input=vas_file, output=aux_file)
    move(aux_file, vas_file)
    vasmn_file = wdir + '/{}_vas_mm.nc'.format(model)
    cdo.selvar('vas', input='-monmean {}'.format(vas_file),
               option='-b F32', output=vasmn_file)
    logger.info('Computing auxiliary variables\n')
    # emission temperature
    te_file = wdir + '/{}_te.nc'.format(model)
    cdo.sqrt(input="-sqrt -mulc,{} {}".format(SIGMAINV, rlut_file),
             output=te_file)
    te_ymm_file = wdir + '/{}_te_ymm.nc'.format(model)
    cdo.yearmonmean(input=te_file, output=te_ymm_file)
    te_gmean_file = wdir + '/{}_te_gmean.nc'.format(model)
    cdo.timmean(input='-fldmean {}'.format(te_ymm_file), output=te_gmean_file)
    f_l = Dataset(te_gmean_file)
    te_gmean_constant = f_l.variables['rlut'][0, 0, 0]
    logger.info('Global mean emission temperature: %s\n',
                te_gmean_constant)
    # temperature of the atmosphere-surface interface
    tasvert_file = wdir + '/{}_tvertavg.nc'.format(model)
    removeif(tasvert_file)
    cdo.fldmean(input='-mulc,0.5 -add {} -selvar,tas {}'
                .format(ts_file, tasmn_file), options='-b F32',
                output=tasvert_file)
    # evaporation from latent heat fluxes at the surface
    if wat in {'y', 'yes'} and entr in {'n', 'no'}:
        evspsbl_file, prr_file = comp_wfluxes(model, wdir, filelist)
        aux_files = [evspsbl_file, prr_file]
    elif entr in {'y', 'yes'}:
        if met in {'2', '3'}:
            evspsbl_file, prr_file = comp_wfluxes(model, wdir, filelist)
            mk_list = [ts_file, hus_file, ps_file, uasmn_file, vasmn_file,
                       hfss_file, te_file]
            tabl_file, tlcl_file, htop_file = mkth.mkthe_main(wdir, mk_list,
                                                              model)
            # Working temperatures for the hydrological cycle
            tcloud_file = (wdir + '/{}_tcloud.nc'.format(model))
            removeif(tcloud_file)
            cdo.mulc('0.5', input='-add {} {}'.format(tlcl_file, te_file),
                     options='-b F32', output=tcloud_file)
            tcolumn_file = (wdir + '/{}_t_vertav_pot.nc'.format(model))
            removeif(tcolumn_file)
            cdo.mulc('0.5', input='-add {} {}'.format(ts_file, tcloud_file),
                     options='-b F32', output=tcolumn_file)
            # Working temperatures for the kin. en. diss. (updated)
            tasvert_file = (wdir + '/{}_tboundlay.nc'.format(model))
            removeif(tasvert_file)
            cdo.fldmean(input='-mulc,0.5 -add {} {}'
                        .format(ts_file, tabl_file), options='-b F32',
                        output=tasvert_file)
            aux_files = [evspsbl_file, htop_file, prr_file, tabl_file,
                         tasvert_file, tcloud_file, tcolumn_file,
                         tlcl_file]
        else:
            pass
    else:
        pass
    return te_ymm_file, te_gmean_constant, te_file, aux_files


def comp_baroceff(model, wdir, aux_file, toab_file, te_file):
    """Compute the baroclinic efficiency of the atmosphere.

    The function computes the baroclinic efficiency of the atmosphere, i.e.
    the efficiency of the meridional heat transports from the low latitudes,
    where there is a net energy gain, towards the high latitudes, where there
    is a net energy loss (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - aux_file: the name of a dummy aux. file to be used for computations;
    - toab_file: a file containing the annual mean TOA energy budgets
      (time,lon,lat);
    - te_file: a file containing the annual mean emission temperature
      (time,lon,lat);

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo
    removeif(aux_file)
    gain_file = wdir + '/{}_maskGain.nc'.format(model)
    cdo.gtc('0', input=toab_file, output=gain_file)
    loss_file = wdir + '/{}_maskLoss.nc'.format(model)
    cdo.ltc('0', input=toab_file, output=loss_file)
    toabgain_file = wdir + '/{}_toabGain.nc'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'.format(toab_file, gain_file),
                   output=toabgain_file)
    toabloss_file = wdir + '/{}_toabLoss.nc'.format(model)
    cdo.setrtomiss('0,1000', input='-mul {} {}'.format(toab_file, loss_file),
                   output=toabloss_file)
    tegain_file = wdir + '/{}_teGain.nc'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'.format(te_file, gain_file),
                   output=tegain_file)
    teloss_file = wdir + '/{}_teLoss.nc'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'.format(te_file, loss_file),
                   output=teloss_file)
    tegainm_file = wdir + '/{}_teGainm.nc'.format(model)
    cdo.div(input='-fldmean {} -fldmean -div {} {} '
            .format(toabgain_file, toabgain_file, tegain_file),
            output=tegainm_file)
    telossm_file = wdir + '/{}_teLossm.nc'.format(model)
    cdo.div(input='-fldmean {} -fldmean -div {} {} '
            .format(toabloss_file, toabloss_file, teloss_file),
            output=telossm_file)
    aux_baroceff_file = (wdir + '/{}_aux_barocEff.nc'.format(model))
    cdo.sub(input='-reci {} -reci {}'.format(telossm_file, tegainm_file),
            output=aux_baroceff_file)
    baroceff_file = wdir + '/{}_barocEff.nc'.format(model)
    cdo.div(input='{} -mulc,0.5 -add -reci {} -reci {}'
            .format(aux_baroceff_file, tegainm_file, telossm_file),
            output=baroceff_file)
    f_l = Dataset(baroceff_file)
    baroc = f_l.variables['toab'][0, 0, 0]
    logger.info('Baroclinic efficiency (Lucarini et al., 2011): %s\n', baroc)
    return baroc


def comp_budgets(model, wdir, aux_file, filelist):
    """Compute radiative budgets from radiative and heat fluxes.

    The function computes TOA and surface energy budgets from radiative and
    heat fluxes, then writes the annual mean to the log info file and write
    the (lat,lon) annual mean fields to a NetCDF file, as well as the time
    series of the annual mean globally averaged fields.

    toab = rsdt - rsut - rlut
    surb = rsds + rlds - rsus - rlus - hfls - hfss
    atmb = toab - atmb

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - aux_file: the name of a dummy aux. file to be used for computations;
    - filelist: a list of file names containing the input fields;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo
    hfls_file = filelist[0]
    hfss_file = filelist[1]
    rlds_file = filelist[6]
    rlus_file = filelist[7]
    rlut_file = filelist[8]
    rsds_file = filelist[9]
    rsdt_file = filelist[10]
    rsus_file = filelist[11]
    rsut_file = filelist[12]
    toab_file = wdir + '/{}_toab.nc'.format(model)
    toab_gmean_file = wdir + '/{}_toab_gmean.nc'.format(model)
    surb_file = wdir + '/{}_surb.nc'.format(model)
    aux_surb_file = wdir + '/{}_aux_surb.nc'.format(model)
    surb_gmean_file = wdir + '/{}_surb_gmean.nc'.format(model)
    atmb_file = wdir + '/{}_atmb.nc'.format(model)
    atmb_gmean_file = wdir + '/{}_atmb_gmean.nc'.format(model)
    removeif(aux_file)
    cdo.sub(input="-sub {} {} {}".format(rsdt_file, rsut_file, rlut_file),
            output=aux_file)
    toab_gmean = write_eb('rsdt', 'toab', aux_file, toab_file, toab_gmean_file)
    logger.info('TOA energy budget: %s\n', np.nanmean(toab_gmean))
    toab_ymm_file = wdir + '/{}_toab_ymm.nc'.format(model)
    cdo.yearmonmean(input=toab_file, output=toab_ymm_file)
    # Surface energy budget
    removeif(aux_file)
    cdo.add(input=" {} {}".format(rsds_file, rlds_file), output=aux_surb_file)
    cdo.sub(input="-sub -sub -sub {} {} {} {} {}"
            .format(aux_surb_file, rsus_file, rlus_file, hfls_file, hfss_file),
            output=aux_file)
    surb_gmean = write_eb('rsds', 'surb', aux_file, surb_file, surb_gmean_file)
    logger.info('Surface energy budget: %s\n', np.nanmean(surb_gmean))
    # Atmospheric energy budget
    removeif(aux_file)
    cdo.sub(input="{} {}".format(toab_file, surb_file), output=aux_file)
    atmb_gmean = write_eb('toab', 'atmb', aux_file, atmb_file, atmb_gmean_file)
    eb_gmean = [toab_gmean, atmb_gmean, surb_gmean]
    eb_file = [toab_file, atmb_file, surb_file]
    return eb_gmean, eb_file, toab_ymm_file


def comp_entr(filelist, nin, nout, entr_file, entr_mean_file):
    """Obtain the entropy dividing some energy by some working temperature.

    This function ingests an energy and a related temperature, then writes
    (time,lat,lon) entropy fluxes and entropy flux annual mean values to NC
    files.

    Arguments:
    - filelist: a list of file containing the name of the energy file, of the
      temperature file and of an auxiliary file needed for computation;
    - nin: the variable name of the input energy fields;
    - nout: the variable name to attribute to the entropy flux in the NC file;
    - entr_file: the name of the file containing the 3D entropy fluxes;
    - entr_mean_file: the name of the file containing the global annual mean
      entropy value;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    ch_name = '{},{}'.format(nin, nout)
    en_file = filelist[0]
    tem_file = filelist[1]
    aux_file = filelist[2]
    aux2_file = 'aux.nc'
    removeif(aux_file)
    removeif(aux2_file)
    cdo.timmean(input='-yearmonmean -monmean -div {} {}'
                .format(en_file, tem_file), options='-b F32', output=aux_file)
    entr_gmean = write_eb(nin, nout, aux2_file, entr_file, entr_mean_file)
    cdo.chname(ch_name, input=aux_file, options='-b F32', output=entr_file)
    os.remove(aux2_file)
    return entr_gmean


def comp_evapentr(model, wdir, infile, aux_file):
    """Compute entropy production related to evaporation fluxes.

    The function computes the material entropy production related to
    evaporation fluxes, as part of the material entropy production
    obtained with the direct method (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of file containing hfls and ts, respectively
      (with dimensions (time,lat,lon);
    - aux_file: the name of a dummy aux. file to be used for computations;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    evapentr_file = wdir + '/{}_evap_entr.nc'.format(model)
    evapentr_mean_file = wdir + '/{}_evapEntropy_gmean.nc'.format(model)
    flist = [infile[0], infile[1], aux_file]
    evapentr_gmean = comp_entr(flist, 'tabl', 'sevap', evapentr_file,
                               evapentr_mean_file)
    evapentr_gmean = masktonull(evapentr_gmean)
    return evapentr_gmean, evapentr_file


def comp_indentr(model, wdir, infile, aux_file, toab_gmean):
    """Compute the material entropy production with the indirect method.

    The function computes the material entropy production with the indirect
    method, isolating a vertical and a horizontal component
    (after Lucarini et al., 2011). The outputs are stored in terms of global
    mean time series, and in terms of (lat,lon) fields for each year to a NC
    file.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of files, containing each the fields rlds, rlus, rsds,
      rsus, emission temperature (te), TOA energy budget (toab) and ts;
    - toab_file: a file containing the annual mean TOA energy budgets
      (time,lon,lat);
    - aux_file: the name of a dummy aux. file to be used for computations;
    - toab_gmean: the climatological annaul mean TOA energy budget;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    horzentropy_file = wdir + '/{}_horizEntropy.nc'.format(model)
    vertentropy_file = wdir + '/{}_verticalEntropy.nc'.format(model)
    vertentropy_mean_file = wdir + '/{}_vertEntropy_gmean.nc'.format(model)
    horzentropy_mean_file = wdir + '/{}_horizEntropy_gmean.nc'.format(model)
    removeif(aux_file)
    cdo.yearmonmean(input='-mulc,-1 -div -subc,{}  {}  {}'
                    .format(np.nanmean(toab_gmean), infile[5], infile[4]),
                    output=aux_file)
    horzentr_mean = write_eb('toab', 'shor', aux_file, horzentropy_file,
                             horzentropy_mean_file)
    cdo.yearmonmean(input=' -add {} -sub {} -add {} {}'
                    .format(infile[0], infile[2], infile[1], infile[3]),
                    output=aux_file)
    cdo.mulc('-1', input='-mul  -sub -yearmonmean -reci {} \
             -yearmonmean -reci {} {}'
             .format(infile[6], infile[4], aux_file),
             output=vertentropy_file)
    vertentr_mean = write_eb('ts', 'sver', aux_file, vertentropy_file,
                             vertentropy_mean_file)
    return horzentr_mean, vertentr_mean, vertentropy_file


def comp_kinentr(aux_file, tasvert_file, lect, lec):
    """Compute the material entropy production from kin. energy dissipation

    The function computes the material entropy production associated with the
    kinetic energy dissipation, through the intensity of the LEC.

    Arguments:
    - aux_file: the name of a dummy aux. file to be used for computations;
    - tasvert_file: a file containing the vertically integrated boundary layer
      temperature;
    - lect: an array containing the annual mean LEC intensity;
    - lec: a flag marking whether the LEC has been previously computed or not

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    removeif(aux_file)
    if lec in {'y', 'yes'}:
        cdo.yearmonmean(input=tasvert_file, output=aux_file)
        f_l = Dataset(aux_file)
        tabl_mean = f_l.variables['ts'][:, 0, 0]
        minentr_mean = np.nanmean(lect / tabl_mean)
        logger.info('Material entropy production associated with '
                    'kinetic energy dissipation: %s\n',
                    minentr_mean)
        minentr_mean = masktonull(minentr_mean)
    else:
        minentr_mean = 0.010
        logger.info('I cannot compute the material entropy '
                    'production without the LEC...\n')
        logger.info('I will assign a given value for the material '
                    'entropy production attributed to LEC '
                    '(0.01 W/m2*K)\n')
    return minentr_mean


def comp_landoc_budg(model, wdir, infile, mask, name):
    """Compute budgets separately on land and oceans.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: the file containing the original budget field as (time,lat,lon);
    - mask: the file containing the land-sea mask;
    - name: the variable name as in the input file;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    ocean_file = wdir + '/{}_{}_ocean.nc'.format(model, name)
    oc_gmean_file = wdir + '/{}_{}_oc_gmean.nc'.format(model, name)
    land_file = wdir + '/{}_{}_land.nc'.format(model, name)
    la_gmean_file = wdir + '/{}_{}_la_gmean.nc'.format(model, name)
    aux_file = wdir + 'aux.nc'
    removeif(aux_file)
    cdo.mul(input='{} -eqc,0 {}'.format(infile, mask), output=ocean_file)
    cdo.timmean(input='-fldmean {}'.format(ocean_file),
                output=oc_gmean_file)
    f_l = Dataset(oc_gmean_file)
    oc_gmean = f_l.variables[name][0, 0, 0]
    cdo.sub(input='{} {}'.format(infile, ocean_file), output=land_file)
    cdo.setctomiss('0', input=ocean_file, output=aux_file)
    move(aux_file, ocean_file)
    cdo.setctomiss('0', input=land_file, output=aux_file)
    move(aux_file, land_file)
    cdo.timmean(input='-fldmean {}'.format(land_file),
                output=la_gmean_file)
    f_l = Dataset(la_gmean_file)
    la_gmean = f_l.variables[name][0, 0, 0]
    return oc_gmean, la_gmean


def comp_meltentr(model, wdir, latsnow_file, aux_file):
    """Compute entropy production related to snow melting at the ground.

    The function computes the material entropy production related to snow
    melting at the ground, as part of the material entropy production
    obtained with the direct method (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: the latent energy associated with snowfall precipitation;
    - aux_file: the name of a dummy aux. file to be used for computations;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    removeif(aux_file)
    latmelt_file = (wdir + '/{}_latentEnergy_snowmelt.nc'.format(model))
    meltentr_file = (wdir + '/{}_snowmelt_entr.nc'.format(model))
    meltentr_mean_file = wdir + '/{}_snowmeltEntropy_gmean.nc'.format(model)
    cdo.mulc(str(L_S), input='-divc,{} {}'.format(str(LC_SUB), latsnow_file),
             options='-b F32', output=latmelt_file)
    cdo.timmean(input='-yearmonmean -monmean -setmisstoc,0 -divc,273.15 {}'
                .format(latmelt_file), options='-b F32', output=aux_file)
    cdo.chname('prsn,smelt', input=aux_file, options='-b F32',
               output=meltentr_file)
    cdo.fldmean(input=meltentr_file, options='-b F32',
                output=meltentr_mean_file)
    f_l = Dataset(meltentr_mean_file)
    meltentr_gmean = f_l.variables['smelt'][0, 0, 0]
    meltentr_gmean = masktonull(meltentr_gmean)
    return meltentr_gmean, meltentr_file


def comp_potentr(model, wdir, infile, aux_file):
    """Compute entropy production related to potential energy of the droplet.

    The function computes the material entropy production related to the
    potential energy of the snowfall or rainfall droplet. This term must be
    part of a material entropy production budget, even though it does take part
    to the energy exchanges of a model "normally".

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of files containing the height of the bondary layer top
      (htop), the masked rainfall precipitation (prrmask), the masked snowfall
      precipitation (prsnmask), the temperature of the vertical column between
      the cloud top and the ground (tcolumn);
    - aux_file: the name of a dummy aux. file to be used for computations;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    removeif(aux_file)
    htop_file = infile[0]
    prrmask_file = infile[1]
    prsnmask_file = infile[2]
    tcolumn_file = infile[3]
    poten_file = wdir + '/{}_potentialEnergy_droplet.nc'.format(model)
    potentr_file = wdir + '/{}_pot_drop_entr.nc'.format(model)
    potentr_mean_file = wdir + '/{}_potentialdropEnergy_gmean.nc'.format(model)
    cdo.mulc(GRAV, input='-mul {} -add {} {}'
             .format(htop_file, prrmask_file, prsnmask_file),
             options='-b F32', output=poten_file)
    flist = [poten_file, tcolumn_file, aux_file]
    potentr_gmean = comp_entr(flist, 'htop', 'spotp', potentr_file,
                              potentr_mean_file)
    potentr_gmean = masktonull(potentr_gmean)
    return potentr_gmean, potentr_file


def comp_rainentr(model, wdir, infile, aux_file):
    """Compute entropy production related to rainfall precipitation.

    The function computes the material entropy production related to rainfall
    precipitation, as part of the material entropy production obtained with the
    direct method (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of file containing the masked rainfall precipitation
      (prrmask) and the temperature of the cloud (tcloud);
    - aux_file: the name of a dummy aux. file to be used for computations;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    prrmask_file = infile[0]
    removeif(aux_file)
    latrain_file = wdir + '/{}_latentEnergy_rain.nc'.format(model)
    rainentr_file = wdir + '/{}_rain_entr.nc'.format(model)
    rainentr_mean_file = wdir + '/{}_rainEntropy_gmean.nc'.format(model)
    cdo.mulc(str(L_C), input='-setmisstoc,0 {}'.format(prrmask_file),
             options='-b F32', output=latrain_file)
    flist = [latrain_file, infile[1], aux_file]
    rainentr_gmean = comp_entr(flist, 'prr', 'srain', rainentr_file,
                               rainentr_mean_file)
    rainentr_gmean = masktonull(rainentr_gmean)
    return rainentr_gmean, rainentr_file


def comp_sensentr(model, wdir, infile, aux_file):
    """Compute entropy production related to sensible heat fluxes.

    The function computes the material entropy production related to sensible
    heat fluxes, as part of the material entropy production obtained with the
    direct method (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of file containing hfss, the temperature at the boundary
    layer top (tabl), ts, respectively (with dimensions (time,lat,lon);
    - aux_file: the name of a dummy aux. file to be used for computations;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    difftemp_file = wdir + '/{}_difftemp_bl.nc'.format(model)
    sensentr_file = (wdir + '/{}_sens_entr.nc'.format(model))
    sensentr_mean_file = wdir + '/{}_sensEntropy_gmean.nc'.format(model)
    cdo.reci(input='-sub -reci {}  -reci {}'.format(infile[1], infile[2]),
             options='-b F32', output=difftemp_file)
    flist = [infile[0], difftemp_file, aux_file]
    sensentr_gmean = comp_entr(flist, 'tabl', 'ssens', sensentr_file,
                               sensentr_mean_file)
    sensentr_gmean = masktonull(sensentr_gmean)
    return sensentr_gmean, sensentr_file


def comp_snowentr(model, wdir, infile, aux_file):
    """Compute entropy production related to snowfall precipitation.

    The function computes the material entropy production related to snowfall
    precipitation, as part of the material entropy production obtained with the
    direct method (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of file containing the masked snowfall precipitation
      (prsnmask) and the temperature of the cloud (tcloud);
    - aux_file: the name of a dummy aux. file to be used for computations;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    prsnmask_file = infile[0]
    removeif(aux_file)
    latsnow_file = wdir + '/{}_latentEnergy_snow.nc'.format(model)
    snowentr_file = wdir + '/{}_snow_entr.nc'.format(model)
    snowentr_mean_file = wdir + '/{}_snowEntropy_gmean.nc'.format(model)
    cdo.mulc(str(LC_SUB), input='-setmisstoc,0 {}'.format(prsnmask_file),
             options='-b F32', output=latsnow_file)
    flist = [latsnow_file, infile[1], aux_file]
    snowentr_gmean = comp_entr(flist, 'prsn', 'srain', snowentr_file,
                               snowentr_mean_file)
    snowentr_gmean = masktonull(snowentr_gmean)
    return snowentr_gmean, latsnow_file, snowentr_file


def comp_wfluxes(model, wdir, filelist):
    """Compute auxiliary fields and perform time averaging of existing fields.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - filelist: a list of file names containing the input fields;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo
    hfls_file = filelist[0]
    pr_file = filelist[3]
    prsn_file = filelist[4]
    aux_file = wdir + '/aux.nc'
    evspsbl_file = (wdir + '/{}_evspsbl.nc'.format(model))
    cdo.divc(str(L_C), input="{}".format(hfls_file),
             output=evspsbl_file)
    # Rainfall precipitation
    prr_file = wdir + '/{}_prr.nc'.format(model)
    cdo.sub(input="{} {}".format(pr_file, prsn_file),
            output=aux_file)
    cdo.chname('pr,prr', input=aux_file, output=prr_file)
    logger.info('Done\n')
    return evspsbl_file, prr_file


def comp_wmbudg(model, wdir, aux_file, filelist, flags):
    """Compute the water mass and latent energy budgets.

    This function computes the annual mean water mass and latent energy budgets
    from the evaporation and rainfall/snowfall precipitation fluxes and prints
    them to a NetCDF file.
    The globally averaged annual mean budgets are also provided and saved to
    a NetCDF file.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - aux_file: the name of a dummy aux. file to be used for computations;
    - filelist: a list of file names containing the input fields;
    - wat: a flag for the water mass budget module (y or n);
    - entr: a flag for the material entropy production (y or n);
    - met: a flag for the material entropy production method (1: indirect,
        2, direct, 3: both);

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo
    hfls_file = filelist[0]
    pr_file = filelist[3]
    prsn_file = filelist[4]
    _, _, _, aux_list = auxiliary(model, wdir, filelist, flags)
    evspsbl_file = aux_list[0]
    prr_file = aux_list[1]
    wmbudg_file = wdir + '/{}_wmb.nc'.format(model)
    wm_gmean_file = wdir + '/{}_wmb_gmean.nc'.format(model)
    latene_file = wdir + '/{}_latent.nc'.format(model)
    latene_gmean_file = wdir + '/{}_latent_gmean.nc'.format(model)
    removeif(aux_file)
    cdo.sub(input="{} {}".format(evspsbl_file, pr_file), output=aux_file)
    wmass_gmean = write_eb('hfls', 'wmb', aux_file, wmbudg_file, wm_gmean_file)
    # Latent energy budget
    removeif(aux_file)
    cdo.sub(input="{} -add -mulc,{} {} -mulc,{} {}"
            .format(hfls_file, str(LC_SUB), prsn_file,
                    str(L_C), prr_file),
            output=aux_file)
    latent_gmean = write_eb('hfls', 'latent', aux_file, latene_file,
                            latene_gmean_file)
    varlist = [wmass_gmean, latent_gmean]
    filelist = [wmbudg_file, latene_file]
    return varlist, filelist


def lec_preproc(model, wdir, pdir, filelist):
    """Preprocess fields for LEC computations and send it to lorenz program.

    This function computes the interpolation of ta, ua, va, wap daily fields to
    fill gaps using near-surface data, then computes the Fourier coefficients
    and performs the LEC computations. For every year, (lev,lat,wave) fields,
    global and hemispheric time series of each conversion and reservoir term
    of the LEC is provided.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - pdir: a new directory is created as a sub-directory of the plot directory
      to store tables of conversion/reservoir terms and the flux diagram for
      year;
    - filelist: a list of file names containing the input fields;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    fourc = fourier_coefficients
    lorenz = lorenz_cycle
    ta_file = filelist[13]
    tas_file = filelist[14]
    ua_file = filelist[16]
    uas_file = filelist[17]
    va_file = filelist[18]
    vas_file = filelist[19]
    wap_file = filelist[20]
    ldir = os.path.join(pdir, 'LEC_results')
    os.makedirs(ldir)
    maskorog = wdir + '/orog.nc'
    ua_file_mask = wdir + '/ua_fill.nc'
    va_file_mask = wdir + '/va_fill.nc'
    energy3_file = wdir + '/energy_short.nc'
    cdo.setmisstoc('0', input='-setmisstoc,1 -sub {} {}'
                   .format(ua_file, ua_file), options='-b F32',
                   output=maskorog)
    cdo.add(input=('-setmisstoc,0 -selvar,ua {} '
                   '-setmisstoc,0 -mul {} -selvar,ua {}')
            .format(ua_file, uas_file, maskorog), options='-b F32',
            output=ua_file_mask)
    cdo.add(input=('-setmisstoc,0 -selvar,va {} '
                   '-setmisstoc,0 -mul {} -selvar,ua {}')
            .format(va_file, vas_file, maskorog), options='-b F32',
            output=va_file_mask)
    cdo.setmisstoc('0', input=('-invertlat -sellevel,10000/90000 '
                               '-merge {} {} {} {}')
                   .format(ta_file, ua_file_mask,
                           va_file_mask, wap_file),
                   options='-b F32', output=energy3_file)
    yrs = cdo.showyear(input=energy3_file)
    yrs = str(yrs)
    yrs2 = yrs.split()
    y_i = 0
    lect = np.zeros(len(yrs2))
    for y_r in yrs2:
        y_r = int(filter(str.isdigit, y_r))
        enfile_yr = wdir + '/inputen.nc'
        tasfile_yr = wdir + '/tas_yr.nc'
        tadiag_file = wdir + '/ta_filled.nc'
        ncfile = wdir + '/fourier_coeff.nc'
        cdo.selyear(y_r, input=energy3_file, options='-b F32',
                    output=enfile_yr)
        cdo.selyear(y_r, input=tas_file, options='-b F32',
                    output=tasfile_yr)
        fourc.fourier_coeff(tadiag_file, ncfile, enfile_yr, tasfile_yr)
        diagfile = (ldir + '/{}_{}_lec_diagram.png'.format(model, y_r))
        logfile = (ldir + '/{}_{}_lec_table.txt'.format(model, y_r))
        lect[y_i] = lorenz.lorenz(wdir, model, y_r, ncfile, diagfile, logfile)
        y_i = y_i + 1
        os.remove(enfile_yr)
        os.remove(tasfile_yr)
        os.remove(tadiag_file)
        os.remove(ncfile)
    os.remove(ua_file_mask)
    os.remove(va_file_mask)
    os.remove(energy3_file)
    return lect


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def mask_precip(model, wdir, infile):
    """Mask precipitation according to the phase of the droplet.

    This function mask the rainfall and snowfall precipitation fields, as well
    as in dependency of the temperature of the cloud at the droplet formation.
    This allows to isolate some intermediate phase changes of the droplet life
    cycle in the atmosphere.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of input file, containing rainfall precipitation (prr) and
      prsn, respectively (dimensions (time,lat,lon));

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo()
    prr_file = infile[0]
    prsn_file = infile[1]
    tlcl_file = infile[2]
    # Prepare masks for snowfall and rainfall
    maskrain_file = wdir + '/{}_maskprecr.nc'.format(model)
    cdo.gtc('1.0E-7', input=prr_file, options=' -b F32', output=maskrain_file)
    masksnow_file = wdir + '/{}_maskprecs.nc'.format(model)
    cdo.gtc('1.0E-7', input=prsn_file, options=' -b F32', output=masksnow_file)
    prrmask_file = wdir + '/{}_prr_masked.nc'.format(model)
    cdo.mul(input='{} {}'.format(maskrain_file, prr_file), options='-b F32',
            output=prrmask_file)
    prsnmask_file = wdir + '/{}_prsn_masked.nc'.format(model)
    cdo.mul(input='{} {}'.format(masksnow_file, prsn_file), options='-b F32',
            output=prsnmask_file)
    # Temperatures of the rainfall and snowfall clouds
    tliq_file = wdir + '/{}_tliq.nc'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'
                   .format(tlcl_file, maskrain_file),
                   options='-b F32', output=tliq_file)
    tsol_file = wdir + '/{}_tsol.nc'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'
                   .format(tlcl_file, masksnow_file),
                   options='-b F32', output=tsol_file)
    tdegl_file = wdir + '/{}_tliqdeg.nc'.format(model)
    cdo.subc('273.15', input=tliq_file, options='-b F32', output=tdegl_file)
    tdegs_file = wdir + '/{}_tsoldeg.nc'.format(model)
    cdo.subc('273.15', input=tsol_file, options='-b F32', output=tdegs_file)
    # Mask for ice cloud and temperature for phase changes from ice to rain
    maskice_file = wdir + '/{}_maskice.nc'.format(model)
    cdo.ltc('0.0', input=tdegl_file, options='-b F32', output=maskice_file)
    ticer_file = wdir + '/{}_t_icerain_file'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'
                   .format(tliq_file, maskice_file),
                   options='-b F32', output=ticer_file)
    prrice_file = wdir + '/{}_prr_ice_file.nc'.format(model)
    cdo.mul(input='{} {}'.format(maskice_file, prr_file), options='-b F32',
            output=prrice_file)
    # Mask for vapor cloud and temperature for phase changes from vapor to snow
    maskvap_file = wdir + '/{}_maskvap.nc'.format(model)
    cdo.gtc('0.0', input=tdegs_file, options='-b F32', output=maskvap_file)
    tvaps_file = wdir + '/{}_t_vapsnow.nc'.format(model)
    cdo.setrtomiss('-1000,0', input='-mul {} {}'
                   .format(tsol_file, maskvap_file),
                   options='-b F32', output=tvaps_file)
    prsnvap_file = wdir + '/{}_prsn_vap.nc'.format(model)
    cdo.mul(input='{} {}'.format(maskvap_file, prsn_file), options='-b F32',
            output=prsnvap_file)
    return prrmask_file, prsnmask_file


def masktonull(value):
    """Replace missing values with zeros."""
    try:
        value = float(value)
    except Warning:
        value = 0
    return value


def write_eb(namein, nameout, aux_file, d3_file, gmean_file):
    """Change variable name in the NetCDF file and compute averages.

    Arguments:
    - namein: initial name of the variable;
    - nameout: final name of the variable;
    - aux_file: the name of an auxiliary file;
    - d3_file: the file containing (time,lat,lon) fields;
    - gmean_file: the name of a file where to put the annual and globally
      averaged fields;

    Author:
    Valerio Lembo, University of Hamburg (2019).
    """
    cdo = Cdo
    ch_name = '{},{}'.format(namein, nameout)
    cdo.chname(ch_name, input=aux_file, options='-b F32', output=d3_file)
    cdo.fldmean(input='-yearmonmean {}'.format(d3_file), output=gmean_file)
    f_l = Dataset(gmean_file)
    constant = f_l.variables[nameout][:, :, :]
    return constant


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
