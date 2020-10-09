r"""MAIN PROGRAM.

TheDiaTo - The diagnostic tool for climate system thermodynamics.

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

SOFTWARE DESCRIPTION

The tool consists of three modules; one for the
computation of energy budgets and transports, one for the water mass budgets
(and related meridional transports), one for the Lorenz Energy Cycle (LEC),
one for the material entropy production.

The first module is run by default, the others are optional. If the lsm option
is set to true, the module 1 and the module 2 will be run with additional
separate results over land and oceans. The land-sea mask is provided by the
ESMValTool preprocessor.

- MODULE 1 (default)
Earth's energy budgets from radiative and heat fluxes at Top-of-Atmosphere,
at the surface and in the atmosphere (as a residual).
Meridional transports, magnitude and location of the peaks in each
hemisphere (only for heat transports) are also computed.
The baroclinic efficiency is computed from TOA energy budgets, emission
temperature (in turn retrieved from OLR) and near-surface temperature.

- MODULE 2 (optional)
Water mass and latent energy budgets and meridional transports are computed
from latent heat fluxes, snowfall and rainfall precipitation fluxes. Magnitude
and location of the peaks in each hemisphere (only for heat transports) are
also computed, as for module 1.

- MODULE 3 (optional)
The Lorenz Energy Cycle (LEC) is computed in spectral components from near-
surface temperatures, temperatures and the three components of velocities
over pressure levels.
The storage and conversion terms are directly computed, the sources and
sinks are retrieved as residuals.
Components are grouped into a zonal mean, stationary and transient eddy
part.

- MODULE 4 (optional)
The material entropy production is computed using the indirect method, the
direct method or both (following Lucarini et al., 2014).
For the indirect method a vertical and a horizontal component are provided.
For the direct method, all components are combined, related to the
hydrological cycle (attributable to evaporation, rainfall and snowfall
precipitation, phase changes and potential energy of the droplet), to the
sensible heat fluxes and to kinetic energy dissipation. For the latter the
LEC computation is required, given that the strength of the LEC can be
considered as equal to the kinetic energy dissipated to heating. If the option
for module 3 is set to false, a reference value for the material entropy
production related to the kinetic energy dissipation is provided.

PREREQUISITES

The program shares the same prerequisites with the overall ESMValTool
architecture
(see https://docs.esmvaltool.org/en/latest/quickstart/installation.html)

USAGE

1: Obtain the datasets: the program accepts the following variables as
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
     Daily mean resolution or higher:
     - Near-surface temperature;
     - Near-surface (or 10m) zonal velocity;
     - Near-surface (or 10m) meridional velocity;
     - Air temperature (on pressure levels);
     - Horizontal velocity (on pressure levels);
     - Meridional velocity (on pressure levels);
     - Vertical velocity (on pressure levels);
   Data on lonlat grid are accepted, with CMOR-compliant coordinate system.
   The pre-processing modules of ESMValTool scheme will take care of
   converting known grids and recognized datasets to CMOR standards. For a
   a list of known formats, see
      https://docs.esmvaltool.org/en/latest/input.html#observations

2: A configuration template is available in the ESMValTool release. Set your
  own paths to local directories here. Input datasets are read in MODELPATH,
  MODELPATH2, OBSPATH or OBSPATH2 output datasets are stored in WORKPATH,
  plots in PLOTPATH (refer to the manual for ESMValTool).

3: Go to the recipe file in ~/recipes/recipe_thermodyn_diagtool.yml.
   Set the namelist with the datasets that you neeed, following the ESMValTool
   naming convention. Here you can also set the length of the dataset you want
   to subset.
   In the 'scripts' section, set the options with the modules that you want the
   program to use:
       - wat: if set to true, the program will compute the water mass and
              latent energy budget,
       - lec: if set to true, the program will compute the Lorenz Energy Cycle
              (LEC) averaged on each year;
       - entr: if set to true, the program will compute the material entropy
               production (MEP);
       - met: if set to 1, the program will compute the MEP with the indirect
              method, if set to 2 with the direct method, if set to 3, both
              methods will be computed and compared with each other;
4: Run the tool by typing:
         esmvaltool run esmvaltool/recipes/recipe_thermodyn_diagtool.yml

OUTPUT

The output directory contains the following NetCDF files:
    - (output directory):
        atmos_transp_mean_<model_name>.nc
        latent_transp_mean_<model_name>.nc
        ocean_transp_mean_<model_name>.nc
        total_transp_mean_<model_name>.nc
        wmb_transp_mean_<model_name>.nc

        contain annual mean meridional sections of heat transports in the
        atmosphere, oceans, and as a total; latent energy transports and water
        mass transports;

    - (output directory)/<model_name>:
        <model-name>_atmb.nc
        (<model-name>_latent.nc; if wat is set to true)
        <model-name>_surb.nc
        <model-name>_toab.nc
        (<model-name>_wmb.nc; is wat is set to true)

        contain annual mean 2D fields of energy budget, latent heat and water
        mass budgets;

        <model-name>_barocEff.nc

        contains the evolution of annual mean baroclinic efficiency
        (Lucarini et al., 2011).

        (if entr is set to true):
        <model-name>_evap_entr.nc (if met is set to 2 or 3)
        <model-name>_horizEntropy.nc (if met is set to 1 or 3)
        <model-name>_pot_drop_entr.nc (if met is set to 2 or 3)
        <model-name>_rain_entr.nc (if met is set to 2 or 3)
        <model-name>_sens_entr.nc (if met is set to 2 or 3)
        <model-name>_snow_entr.nc (if met is set to 2 or 3)
        <model-name>_snowmelt_entr.nc (if met is set to 2 or 3)
        <model-name>_verticalEntropy.nc (if met is set to 1 or 3)
        contain the evolution of annual mean components of the material entropy
        production.

    - (plots directory):
        meridional_transp.png: contains the model inter-comparison of total,
        atmospheric and oceanic of meridional sections in zonally averaged
        meridional heat transports;
        scatters_summary.png: contains the scatter plots of
        model intercomparisons of various metrics retrieved in the program;
        scattes_variability: contains scatter plots of model intercomparisons
        between TOA, atmospheric and surface global mean energy budgets and
        their inter-annual variability;

    - (plots directory)/<model-name>:
        <model-name>_atmb_timeser.png: the atmospheric budget annual mean
        global and hemispheric time series;
        <model-name>_energy_climap.png: the TOA, atmospheric and surface
        climatological mean fields;
        <model-name>_heat_transp.png: the meridional sections of total,
        atmospheric and oceanic meridional heat transports (implied from energy
        budgets);
        <model-name>_latent_climap.png: the climatological mean latent heat
        field;
        <model-name>_latent_timeser.png: the latent heat annual mean global and
        hemispheric evolutions;
        <model-name>_latent_transp.png: the meridional section of annual mean
        meridional latent heat transport;
        <model-name>_scatpeak.png: the scatter plots of atmospheric vs. oceanic
        peak magnitude in both hemispheres;
        <model-name>_sevap_climap.png: the annual mean field of material
        entropy production due to evaporation;
        <model-name>_smelt_climap.png: the annual mean field of material
        entropy production due to snow melting;
        <model-name>_spotp_climap.png: the annual mean field of material
        entropy production due to potential energy of the droplet;
        <model-name>_srain_climap.png: the annual mean field of material
        entropy production due to rainfall precipitation;
        <model-name>_ssens_climap.png: the annual mean field of material
        entropy production due to sensible heat fluxes;
        <model-name>_ssnow_climap.png: the annual mean field of material
        entropy production due to snowfall precipitation;
        <model-name>_surb_timeser.png: the surface budget annual mean
        global and hemispheric time series;
        <model-name>_sver_climap.png: the annual mean field of vertical
        material entropy production through the indirect method;
        <model-name>_toab_timeser.png: the TOA budget annual mean
        global and hemispheric time series;
        <model-name>_wmb_climap.png: the climatological mean water mass budget
        field;
        <model-name>_wmb_timeser.png: the water mass annual mean global and
        hemispheric evolutions;
        <model-name>_wmb_transp.png: the meridional section of annual mean
        meridional water mass transport;

    - (plots directory)/<model-name>/LEC_results:
        <model-name>_<year>_lec_diagram.png: the flux diagram for the annual
        mean LEC cycle in a specific year;
        <model-name>_<year>_lec_table.txt: the table containing the storage and
        conversion terms for the annual mean LEC cycle in a specific year;

The file log.txt in the '$WORK_PATH/recipe_thermodyn_diagtool_date_hour/run'
sub-directory contains the values for the metrics and all useful information
for immediate model intercomparison.


#############################################################################

20170803-lembo_valerio: modified header with description and caveats
20170629-koldunov_nikolay: atmospheric budgets diagnostics written
20180524-lembo_valerio: first complete working thermodynamics diagnostics
20190325-lembo_valerio: complete updated version for ESMValTool v2.0b

#############################################################################
"""

# New packages for version 2.0 of ESMValTool
import logging
import os
import warnings

import numpy as np

import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared import ProvenanceLogger
from esmvaltool.diag_scripts.thermodyn_diagtool import (computations,
                                                        lorenz_cycle, mkthe,
                                                        plot_script,
                                                        provenance_meta)

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
logger = logging.getLogger(os.path.basename(__file__))


def compute_water_mass_budget(cfg, wdir_up, pdir, model, wdir, input_data,
                              flags, aux_file):
    logger.info('Computing water mass and latent energy budgets\n')
    aux_list = mkthe.init_mkthe_wat(model, wdir, input_data, flags)
    wm_gmean, wm_file = computations.wmbudg(model, wdir, aux_file, input_data,
                                            aux_list)
    wm_time_mean = np.nanmean(wm_gmean[0])
    wm_time_std = np.nanstd(wm_gmean[0])
    logger.info('Water mass budget: %s\n', wm_time_mean)
    latent_time_mean = np.nanmean(wm_gmean[1])
    latent_time_std = np.nanstd(wm_gmean[1])
    logger.info('Latent energy budget: %s\n', latent_time_mean)
    logger.info('Done\n')
    logger.info('Plotting the water mass and latent energy budgets\n')
    plot_script.balances(cfg, wdir_up, pdir, [wm_file[0], wm_file[1]],
                         ['wmb', 'latent'], model)
    logger.info('Done\n')
    for filen in aux_list:
        os.remove(filen)
    return (wm_file, wm_time_mean, wm_time_std, latent_time_mean,
            latent_time_std)


def compute_land_ocean(model, wdir, file, sftlf_fx, name):
    ocean_mean, land_mean = computations.landoc_budg(model, wdir, file,
                                                     sftlf_fx, name)
    logger.info('%s budget over oceans: %s\n', name, ocean_mean)
    logger.info('%s budget over land: %s\n', name, land_mean)
    return (ocean_mean, land_mean)


def main(cfg):
    """Execute the program.

    Argument cfg, containing directory paths, preprocessed input dataset
    filenames and user-defined options, is passed by ESMValTool preprocessor.
    """
    provlog = ProvenanceLogger(cfg)
    lorenz = lorenz_cycle
    comp = computations
    logger.info('Entering the diagnostic tool')
    # Load paths
    wdir_up = cfg['work_dir']
    pdir_up = cfg['plot_dir']
    input_data = cfg['input_data'].values()
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
    flags = [wat, lec, entr, met]
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
        aux_file = wdir + '/aux.nc'
        te_ymm_file, te_gmean_constant, te_file = mkthe.init_mkthe_te(
            model, wdir, input_data)
        te_all[i_m] = te_gmean_constant
        logger.info('Computing energy budgets\n')
        in_list, eb_gmean, eb_file, toab_ymm_file = comp.budgets(
            model, wdir, aux_file, input_data)
        prov_rec = provenance_meta.get_prov_map(
            ['TOA energy budgets', model],
            [in_list[4], in_list[6], in_list[7]])
        provlog.log(eb_file[0], prov_rec)
        prov_rec = provenance_meta.get_prov_map(
            ['atmospheric energy budgets', model], [
                in_list[0], in_list[1], in_list[2], in_list[3], in_list[4],
                in_list[5], in_list[6], in_list[7], in_list[8]
            ])
        provlog.log(eb_file[1], prov_rec)
        prov_rec = provenance_meta.get_prov_map(
            ['surface energy budgets', model], [
                in_list[0], in_list[1], in_list[2], in_list[3], in_list[5],
                in_list[7]
            ])
        provlog.log(eb_file[2], prov_rec)
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
        plotsmod.balances(cfg, wdir_up, pdir,
                          [eb_file[0], eb_file[1], eb_file[2]],
                          ['toab', 'atmb', 'surb'], model)
        logger.info('Done\n')
        # Water mass budget
        if wat == 'True':
            (wm_file,
             wmb_all[i_m, 0],
             wmb_all[i_m, 1],
             latent_all[i_m, 0],
             latent_all[i_m, 1]) = compute_water_mass_budget(
                 cfg, wdir_up, pdir, model, wdir, input_data, flags, aux_file)
        if lsm == 'True':
            sftlf_fx = e.select_metadata(input_data,
                                         short_name='sftlf',
                                         dataset=model)[0]['filename']
            logger.info('Computing energy budgets over land and oceans\n')
            toab_oc_all[i_m], toab_la_all[i_m] = compute_land_ocean(
                model, wdir, eb_file[0], sftlf_fx, 'toab')
            atmb_oc_all[i_m], atmb_la_all[i_m] = compute_land_ocean(
                model, wdir, eb_file[1], sftlf_fx, 'atmb')
            surb_oc_all[i_m], surb_la_all[i_m] = compute_land_ocean(
                model, wdir, eb_file[2], sftlf_fx, 'surb')
            if wat == 'True':
                logger.info('Computing water mass and latent energy'
                            ' budgets over land and oceans\n')
                wmb_oc_all[i_m], wmb_la_all[i_m] = compute_land_ocean(
                    model, wdir, wm_file[0], sftlf_fx, 'wmb')
                latent_oc_all[i_m], latent_la_all[i_m] = compute_land_ocean(
                    model, wdir, wm_file[1], sftlf_fx, 'latent')
            logger.info('Done\n')
        if lec == 'True':
            logger.info('Computation of the Lorenz Energy '
                        'Cycle (year by year)\n')
            _, _ = mkthe.init_mkthe_lec(model, wdir, input_data)
            lect = lorenz.preproc_lec(model, wdir, pdir, input_data)
            lec_all[i_m, 0] = np.nanmean(lect)
            lec_all[i_m, 1] = np.nanstd(lect)
            logger.info(
                'Intensity of the annual mean Lorenz Energy '
                'Cycle: %s\n', lec_all[i_m, 0])
            logger.info('Done\n')
        else:
            lect = np.repeat(2.0, len(eb_gmean[0]))
            lec_all[i_m, 0] = 2.0
            lec_all[i_m, 1] = 0.2
        if entr == 'True':
            if met in {'1', '3'}:
                logger.info('Computation of the material entropy production '
                            'with the indirect method\n')
                indentr_list = [te_file, eb_file[0]]
                horz_mn, vert_mn, horzentr_file, vertentr_file = comp.indentr(
                    model, wdir, indentr_list, input_data, aux_file,
                    eb_gmean[0])
                listind = [horzentr_file, vertentr_file]
                provenance_meta.meta_indentr(cfg, model, input_data, listind)
                horzentr_all[i_m, 0] = np.nanmean(horz_mn)
                horzentr_all[i_m, 1] = np.nanstd(horz_mn)
                vertentr_all[i_m, 0] = np.nanmean(vert_mn)
                vertentr_all[i_m, 1] = np.nanstd(vert_mn)
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
                    logger, model, wdir, input_data, aux_file, te_file, lect,
                    flags)
                provenance_meta.meta_direntr(cfg, model, input_data, entr_list)
                matentr_all[i_m, 0] = matentr
                if met in {'3'}:
                    diffentr = (float(np.nanmean(vert_mn)) +
                                float(np.nanmean(horz_mn)) - matentr)
                    logger.info('Difference between the two '
                                'methods: %s\n', diffentr)
                    diffentr_all[i_m, 0] = diffentr
                logger.info('Degree of irreversibility of the '
                            'system: %s\n', irrevers)
                irrevers_all[i_m] = irrevers
                logger.info('Running the plotting module for the material '
                            'entropy production (direct method)\n')
                plotsmod.init_plotentr(model, pdir, entr_list)
                logger.info('Done\n')
            os.remove(te_file)
        os.remove(te_ymm_file)
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
