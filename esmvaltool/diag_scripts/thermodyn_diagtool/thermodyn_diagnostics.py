"""MAIN PROGRAM.

The diagnostic tool for climate system thermodynamics.

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
surb: surface energy budget, latent: latent energy budget, wmb: water mass
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

# New packages for version 2.0 of ESMValTool
import logging
import os
import warnings

import matplotlib
import numpy as np
from netCDF4 import Dataset

import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger
# Locally used modules
from esmvaltool.diag_scripts.thermodyn_diagtool import (
    computations, lorenz_cycle, mkthe, plot_script)

matplotlib.use('Agg')
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Execute the program.

    Argument cfg, containing directory paths, preprocessed input dataset
    filenames and user-defined options, is passed by ESMValTool preprocessor.
    """
    lorenz = lorenz_cycle
    comp = computations
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
        os.makedirs(wdir)
        os.makedirs(pdir)
        # Reading file names for the specific model
        filenames = data.get_info_list('filename', dataset=model)
        logger.info('Processing model: %s \n', model)
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
        te_ymm_file, te_gmean_constant, _, _ = mkthe.init_mkthe(
            model, wdir, filenames, flags)
        te_all[i_m] = te_gmean_constant
        logger.info('Computing energy budgets\n')
        eb_gmean, eb_file, toab_ymm_file = comp.budgets(
            model, wdir, aux_file, filenames)
        toab_all[i_m, 0] = np.nanmean(eb_gmean[0])
        toab_all[i_m, 1] = np.nanstd(eb_gmean[0])
        atmb_all[i_m, 0] = np.nanmean(eb_gmean[1])
        atmb_all[i_m, 1] = np.nanstd(eb_gmean[1])
        surb_all[i_m, 0] = np.nanmean(eb_gmean[2])
        surb_all[i_m, 1] = np.nanstd(eb_gmean[2])
        logger.info('Global mean emission temperature: %s\n',
                    te_gmean_constant)
        logger.info('TOA energy budget: %s\n', toab_all[i_m, 0])
        logger.info('Atmospheric energy budget: %s\n', atmb_all[i_m, 0])
        logger.info('Surface energy budget: %s\n', surb_all[i_m, 0])
        logger.info('Done\n')
        baroc_eff_all[i_m] = comp.baroceff(model, wdir, aux_file,
                                           toab_ymm_file, te_ymm_file)
        logger.info('Baroclinic efficiency (Lucarini et al., 2011): %s\n',
                    baroc_eff_all[i_m])
        logger.info('Running the plotting module for the budgets\n')
        plotsmod.balances(wdir_up, pdir, [eb_file[0], eb_file[1], eb_file[2]],
                          ['toab', 'atmb', 'surb'], model)
        logger.info('Done\n')
        # Water mass budget
        if wat is True:
            logger.info('Computing water mass and latent energy budgets\n')
            _, _, _, aux_list = mkthe.init_mkthe(model, wdir, filenames, flags)
            wm_gmean, wm_file = comp.wmbudg(model, wdir, aux_file, filenames,
                                            aux_list)
            wmb_all[i_m, 0] = np.nanmean(wm_gmean[0])
            wmb_all[i_m, 1] = np.nanstd(wm_gmean[0])
            logger.info('Water mass budget: %s\n', wmb_all[i_m, 0])
            latent_all[i_m, 0] = np.nanmean(wm_gmean[1])
            latent_all[i_m, 1] = np.nanstd(wm_gmean[1])
            logger.info('Latent energy budget: %s\n', latent_all[i_m, 0])
            logger.info('Done\n')
            logger.info('Plotting the water mass and latent energy budgets\n')
            plotsmod.balances(wdir_up, pdir, [wm_file[0], wm_file[1]],
                              ['wmb', 'latent'], model)
            logger.info('Done\n')
        else:
            pass
        if lsm is True:
            logger.info('Computing energy budgets over land and oceans\n')
            toab_oc_gmean, toab_la_gmean = comp.landoc_budg(
                model, wdir, eb_file[0], sftlf_fx, 'toab')
            toab_oc_all[i_m] = toab_oc_gmean
            toab_la_all[i_m] = toab_la_gmean
            logger.info('TOA energy budget over oceans: %s\n', toab_oc_gmean)
            logger.info('TOA energy budget over land: %s\n', toab_la_gmean)
            atmb_oc_gmean, atmb_la_gmean = comp.landoc_budg(
                model, wdir, eb_file[1], sftlf_fx, 'atmb')
            atmb_oc_all[i_m] = atmb_oc_gmean
            atmb_la_all[i_m] = atmb_la_gmean
            logger.info('Atmospheric energy budget over oceans: %s\n',
                        atmb_oc_gmean)
            logger.info('Atmospheric energy budget over land: %s\n',
                        atmb_la_gmean)
            surb_oc_gmean, surb_la_gmean = comp.landoc_budg(
                model, wdir, eb_file[2], sftlf_fx, 'surb')
            surb_oc_all[i_m] = surb_oc_gmean
            surb_la_all[i_m] = surb_la_gmean
            logger.info('Surface energy budget over oceans: %s\n',
                        surb_oc_gmean)
            logger.info('Surface energy budget over land: %s\n', surb_la_gmean)
            logger.info('Done\n')
            if wat is True:
                logger.info('Computing water mass and latent energy'
                            ' budgets over land and oceans\n')
                wmb_oc_gmean, wmb_la_gmean = comp.landoc_budg(
                    model, wdir, wm_file[0], sftlf_fx, 'wmb')
                wmb_oc_all[i_m] = wmb_oc_gmean
                wmb_la_all[i_m] = wmb_la_gmean
                logger.info('Water mass budget over oceans: %s\n',
                            wmb_oc_gmean)
                logger.info('Water mass budget over land: %s\n', wmb_la_gmean)
                latent_oc_gmean, latent_la_gmean = comp.landoc_budg(
                    model, wdir, wm_file[1], sftlf_fx, 'latent')
                latent_oc_all[i_m] = latent_oc_gmean
                latent_la_all[i_m] = latent_la_gmean
                logger.info('Latent energy budget over oceans: %s\n',
                            latent_oc_gmean)
                logger.info('Latent energy budget over land: %s\n',
                            latent_la_gmean)
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        if lec is True:
            logger.info('Computation of the Lorenz Energy'
                        'Cycle (year by year)\n')
            lect = lorenz.preproc_lec(model, wdir, pdir, filenames)
            lec_all[i_m, 0] = np.nanmean(lect)
            lec_all[i_m, 1] = np.nanstd(lect)
            logger.info(
                'Intensity of the annual mean Lorenz Energy '
                'Cycle: %s\n', lec_all[i_m, 0])
            logger.info('Done\n')
        else:
            pass
        if entr is True:
            if met in {'1', '3'}:
                _, _, te_file, _ = mkthe.init_mkthe(model, wdir, filenames,
                                                    flags)
                logger.info('Computation of the material entropy production '
                            'with the indirect method\n')
                indentr_list = [
                    rlds_file, rlus_file, rsds_file, rsus_file, te_file,
                    eb_file[0], ts_file
                ]
                horzentr_mean, vertentr_mean, vertentr_file = comp.indentr(
                    model, wdir, indentr_list, aux_file, eb_gmean[0])
                horzentr_all[i_m, 0] = np.nanmean(horzentr_mean)
                horzentr_all[i_m, 1] = np.nanstd(horzentr_mean)
                vertentr_all[i_m, 0] = np.nanmean(vertentr_mean)
                vertentr_all[i_m, 1] = np.nanstd(vertentr_mean)
                logger.info(
                    'Horizontal component of the material entropy '
                    'production: %s\n', horzentr_all[i_m, 0])
                logger.info(
                    'Vertical component of the material entropy '
                    'production: %s\n', vertentr_all[i_m, 0])
                logger.info('Done\n')
                logger.info('Running the plotting module for the material '
                            'entropy production (indirect method)\n')
                plotsmod.entropy(pdir, vertentr_file, 'sver',
                                 'Vertical entropy production', model)
                logger.info('Done\n')
            if met in {'2', '3'}:
                matentr, irrevers, entr_list = comp.direntr(
                    logger, model, wdir, filenames, aux_file, lect, lec, flags)
                matentr_all[i_m, 0] = matentr
                if met in {'3'}:
                    diffentr = (float(np.nanmean(vertentr_mean)) + float(
                        np.nanmean(horzentr_mean)) - matentr)
                    logger.info('Difference between the two '
                                'methods: %s\n', diffentr)
                    diffentr_all[i_m, 0] = diffentr
                else:
                    pass
                logger.info('Degree of irreversibility of the '
                            'system: %s\n', irrevers)
                irrevers_all[i_m] = irrevers
                logger.info('Running the plotting module for the material '
                            'entropy production (direct method)\n')
                plotsmod.init_plotentr(model, pdir, entr_list)
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        logger.info('Done for model: %s \n', model)
        i_m = i_m + 1
    logger.info('I will now start multi-model plots')
    logger.info('Meridional heat transports\n')
    plotsmod.plot_mm_transp(model_names, wdir_up, pdir_up)
    logger.info('Scatter plots')
    summary_varlist = [
        atmb_all, baroc_eff_all, horzentr_all, lec_all, matentr_all, te_all,
        toab_all, vertentr_all
    ]
    plotsmod.plot_mm_summaryscat(pdir_up, summary_varlist)
    logger.info('Scatter plots for inter-annual variability of'
                ' some quantities')
    eb_list = [toab_all, atmb_all, surb_all]
    plotsmod.plot_mm_ebscatter(pdir_up, eb_list)
    logger.info("The diagnostic has finished. Now closing...\n")


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
