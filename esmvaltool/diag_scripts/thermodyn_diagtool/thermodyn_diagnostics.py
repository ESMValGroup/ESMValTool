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
# from esmvaltool.diag_scripts.shared import run_diagnostic
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
# from esmvaltool.diag_scripts.thermodyn_diagtool import fourier_coefficients
# from esmvaltool.diag_scripts.thermodyn_diagtool import lorenz_cycle
# from esmvaltool.diag_scripts.thermodyn_diagtool import plot_script

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
    workdir = cfg['work_dir']
    plotdir = cfg['plot_dir']
    logger.info('Work directory: %s \n', workdir)
    logger.info('Plot directory: %s \n', plotdir)
    diagworkdir_up = workdir
    cdo = Cdo()
    fourc = fourier_coefficients
    lorenz = lorenz_cycle
    mkth = mkthe
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
    e_b = str(cfg['eb'])
    lsm = str(cfg['lsm'])
    wat = str(cfg['wat'])
    lec = str(cfg['lec'])
    entr = str(cfg['entr'])
    met = str(cfg['met'])
    logger.info(e_b)
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
    wmass_all = np.zeros([modnum, 2])
    wmass_oc_all = np.zeros(modnum)
    wmass_la_all = np.zeros(modnum)
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
    for model_name in model_names:
        # Load paths to individual models output and plotting directories
        diagworkdir = os.path.join(diagworkdir_up, model_name)
        plotpath2 = os.path.join(plotdir, model_name)
        if not os.path.exists(diagworkdir):
            os.makedirs(diagworkdir)
        # Clean output directories, if already existing
        if not os.path.exists(plotpath2):
            os.makedirs(plotpath2)
        else:
            for dirs, files in os.walk(diagworkdir):
                for name in files:
                    file_path = os.path.join(diagworkdir, name)
                for name in dirs:
                    file_path = os.path.join(diagworkdir, name)
            for dirs, files in os.walk(plotpath2):
                for name in files:
                    file_path = os.path.join(plotpath2, name)
                    os.remove(file_path)
                for name in dirs:
                    file_path = os.path.join(plotpath2, name)
                    os.rmdir(file_path)
        # Reading file names for the specific model
        filenames = data.get_info_list('filename', dataset=model_name)
        logger.info('Processing model: %s \n', model_name)
        hfls_file = filenames[0]
        hfss_file = filenames[1]
        hus_file = filenames[2]
        pr_file = filenames[3]
        prsn_file = filenames[4]
        ps_file = filenames[5]
        rlds_file = filenames[6]
        rlus_file = filenames[7]
        rlut_file = filenames[8]
        rsds_file = filenames[9]
        rsdt_file = filenames[10]
        rsus_file = filenames[11]
        rsut_file = filenames[12]
        ta_file = filenames[13]
        tas_file = filenames[14]
        ts_file = filenames[15]
        ua_file = filenames[16]
        uas_file = filenames[17]
        va_file = filenames[18]
        vas_file = filenames[19]
        wap_file = filenames[20]
        # Read path to land-sea mask
        head = Dataset(ta_file)
        info = getattr(head, 'metadata')
        attr = info.split()
        diction = {attr[i_m]: attr[i_m + 1] for i_m in range(len(attr) - 1)}
        sftlf_fx = str(diction['{sftlf:'])
        sftlf_fx = sftlf_fx.replace("}", "")
        # Compute monthly mean fields from 2D surface daily fields
        aux_file = diagworkdir + '/{}_aux.nc'.format(model_name)
        cdo.selvar('tas', input=tas_file, output=aux_file)
        move(aux_file, tas_file)
        tasmn_file = diagworkdir + '/{}_tas_mm.nc'.format(model_name)
        cdo.selvar('tas', input='-monmean {}'.format(tas_file),
                   option='-b F32', output=tasmn_file)
        cdo.selvar('uas', input=uas_file, output=aux_file)
        move(aux_file, uas_file)
        uasmn_file = diagworkdir + '/{}_uas_mm.nc'.format(model_name)
        cdo.selvar('uas', input='-monmean {}'.format(uas_file),
                   option='-b F32', output=uasmn_file)
        cdo.selvar('vas', input=vas_file, output=aux_file)
        move(aux_file, vas_file)
        vasmn_file = diagworkdir + '/{}_vas_mm.nc'.format(model_name)
        cdo.selvar('vas', input='-monmean {}'.format(vas_file),
                   option='-b F32', output=vasmn_file)
        logger.info('Computing auxiliary variables\n')
        # emission temperature
        te_file = diagworkdir + '/{}_te.nc'.format(model_name)
        cdo.sqrt(input="-sqrt -mulc,{} {}".format(SIGMAINV, rlut_file),
                 output=te_file)
        te_ymm_file = diagworkdir + '/{}_te_ymm.nc'.format(model_name)
        cdo.yearmonmean(input=te_file, output=te_ymm_file)
        te_gmean_file = diagworkdir + '/{}_te_gmean.nc'.format(model_name)
        cdo.timmean(input='-fldmean {}'.format(te_ymm_file),
                    output=te_gmean_file)
        f_l = Dataset(te_gmean_file)
        te_gmean_constant = f_l.variables['rlut'][0, 0, 0]
        logger.info('Global mean emission temperature: %s\n',
                    te_gmean_constant)
        te_all[i_m] = te_gmean_constant
        # temperature of the atmosphere-surface interface
        tasvert_file = diagworkdir + '/{}_tvertavg.nc'.format(model_name)
        removeif(tasvert_file)
        cdo.fldmean(input='-mulc,0.5 -add {} -selvar,tas {}'
                    .format(ts_file, tasmn_file),
                    options='-b F32', output=tasvert_file)
        # evaporation from latent heat fluxes at the surface
        if (wat in {'y', 'yes'} or met in {'2', '3'}):
            evspsbl_file = (diagworkdir
                            + '/{}_evspsbl.nc'.format(model_name))
            cdo.divc(str(L_C), input="{}".format(hfls_file),
                     output=evspsbl_file)
            # Rainfall precipitation
            prr_file = diagworkdir + '/{}_prr.nc'.format(model_name)
            cdo.sub(input="{} {}".format(pr_file, prsn_file),
                    output=aux_file)
            cdo.chname('pr,prr', input=aux_file, output=prr_file)
            logger.info('Done\n')
        else:
            pass
        if entr in {'y', 'yes'}:
            if met in {'2', '3'}:
                file_list = [ts_file, hus_file, ps_file, uasmn_file,
                             vasmn_file, hfss_file, te_file]
                mkth.mkthe_main(diagworkdir, file_list, model_name)
                tlcl_temp = diagworkdir + '/tlcl.nc'
                tlcl_file = diagworkdir + '/{}_tlcl.nc'.format(model_name)
                cdo.setrtomiss('400,1e36', input=tlcl_temp,
                               output=tlcl_file)
                tabl_temp = diagworkdir + '/tabl.nc'
                tabl_file = diagworkdir + '/{}_tabl.nc'.format(model_name)
                cdo.setrtomiss('400,1e36', input=tabl_temp,
                               output=tabl_file)
                htop_temp = diagworkdir + '/htop.nc'
                htop_file = diagworkdir + '/{}_htop.nc'.format(model_name)
                cdo.setrtomiss('12000,1e36', input=htop_temp,
                               output=htop_file)
                # Working temperatures for the hydrological cycle
                tcloud_file = (diagworkdir
                               + '/{}_tcloud.nc'.format(model_name))
                removeif(tcloud_file)
                cdo.mulc('0.5',
                         input='-add {} {}'.format(tlcl_file, te_file),
                         options='-b F32', output=tcloud_file)
                tcolumn_file = (diagworkdir
                                + '/{}_t_vertav_pot.nc'.format(model_name))
                removeif(tcolumn_file)
                cdo.mulc('0.5', input='-add {} {}'.
                         format(ts_file, tcloud_file), options='-b F32',
                         output=tcolumn_file)
                # Working temperatures for the kin. en. diss. (updated)
                tasvert_file = (diagworkdir
                                + '/{}_tboundlay.nc'.format(model_name))
                removeif(tasvert_file)
                cdo.fldmean(input='-mulc,0.5 -add {} {}'
                            .format(ts_file, tabl_file), options='-b F32',
                            output=tasvert_file)
            else:
                pass
        # Compute energy budgets
        if e_b in {'y', 'yes'}:
            logger.info('Computing energy budgets\n')
            # TOA energy budget
            toab_file = diagworkdir + '/{}_toab.nc'.format(model_name)
            removeif(toab_file)
            removeif(aux_file)
            cdo.sub(input="-sub {} {} {}".format(rsdt_file, rsut_file,
                                                 rlut_file), output=aux_file)
            cdo.chname('rsdt,toab', input=aux_file, options='-b F32',
                       output=toab_file)
            toab_gmean_file = (diagworkdir
                               + '/{}_toab_gmean.nc'.format(model_name))
            cdo.fldmean(input='-yearmonmean {}'.format(toab_file),
                        output=toab_gmean_file)
            f_l = Dataset(toab_gmean_file)
            toab_gmean_constant = f_l.variables['toab'][:, :, :]
            toab_all[i_m, 0] = np.nanmean(toab_gmean_constant)
            toab_all[i_m, 1] = np.nanstd(toab_gmean_constant)
            logger.info('TOA energy budget: %s\n', toab_all[i_m, 0])
            toab_ymm_file = diagworkdir + '/{}_toab_ymm.nc'.format(model_name)
            cdo.yearmonmean(input=toab_file, output=toab_ymm_file)
            # Surface energy budget
            surb_file = diagworkdir + '/{}_surb.nc'.format(model_name)
            aux_surb_file = diagworkdir + '/{}_aux_surb.nc'.format(model_name)
            removeif(surb_file)
            removeif(aux_file)
            cdo.add(input=" {} {}".format(rsds_file, rlds_file),
                    output=aux_surb_file)
            cdo.sub(input="-sub -sub -sub {} {} {} {} {}"
                    .format(aux_surb_file, rsus_file, rlus_file,
                            hfls_file, hfss_file),
                    output=aux_file)
            cdo.chname('rsds,surb', input=aux_file, options='-b F32',
                       output=surb_file)
            surb_gmean_file = (diagworkdir
                               + '/{}_surb_gmean.nc'.format(model_name))
            cdo.fldmean(input='-yearmonmean {}'.format(surb_file),
                        output=surb_gmean_file)
            f_l = Dataset(surb_gmean_file)
            surb_gmean_constant = f_l.variables['surb'][:, :, :]
            surb_all[i_m, 0] = np.nanmean(surb_gmean_constant)
            surb_all[i_m, 1] = np.nanstd(surb_gmean_constant)
            logger.info('Surface energy budget: %s\n', surb_all[i_m, 0])
            # Atmospheric energy budget
            atmb_file = diagworkdir + '/{}_atmb.nc'.format(model_name)
            removeif(atmb_file)
            removeif(aux_file)
            cdo.sub(input="{} {}".format(toab_file, surb_file),
                    output=aux_file)
            cdo.chname('toab,atmb', input=aux_file, options='-b F32',
                       output=atmb_file)
            atmb_gmean_file = (diagworkdir
                               + '/{}_atmb_gmean.nc'.format(model_name))
            cdo.fldmean(input='-yearmonmean {}'.format(atmb_file),
                        output=atmb_gmean_file)
            f_l = Dataset(atmb_gmean_file)
            atmb_gmean_constant = f_l.variables['atmb'][:, :, :]
            atmb_all[i_m, 0] = np.nanmean(atmb_gmean_constant)
            atmb_all[i_m, 1] = np.nanstd(atmb_gmean_constant)
            logger.info('Atmospheric energy budget: %s\n', atmb_all[i_m, 0])
            logger.info('Done\n')
            # Baroclinic efficiency
            maskgain_file = diagworkdir + '/{}_maskGain.nc'.format(model_name)
            cdo.gtc('0', input=toab_ymm_file, output=maskgain_file)
            maskloss_file = diagworkdir + '/{}_maskLoss.nc'.format(model_name)
            cdo.ltc('0', input=toab_ymm_file, output=maskloss_file)
            toabgain_file = diagworkdir + '/{}_toabGain.nc'.format(model_name)
            cdo.setrtomiss('-1000,0', input='-mul {} {}'
                           .format(toab_ymm_file, maskgain_file),
                           output=toabgain_file)
            toabloss_file = diagworkdir + '/{}_toabLoss.nc'.format(model_name)
            cdo.setrtomiss('0,1000', input='-mul {} {}'
                           .format(toab_ymm_file, maskloss_file),
                           output=toabloss_file)
            tegain_file = diagworkdir + '/{}_teGain.nc'.format(model_name)
            cdo.setrtomiss('-1000,0', input='-mul {} {}'
                           .format(te_ymm_file, maskgain_file),
                           output=tegain_file)
            teloss_file = diagworkdir + '/{}_teLoss.nc'.format(model_name)
            cdo.setrtomiss('-1000,0', input='-mul {} {}'
                           .format(te_ymm_file, maskloss_file),
                           output=teloss_file)
            tegainm_file = diagworkdir + '/{}_teGainm.nc'.format(model_name)
            cdo.div(input='-fldmean {} -fldmean -div {} {} '
                    .format(toabgain_file, toabgain_file, tegain_file),
                    output=tegainm_file)
            telossm_file = diagworkdir + '/{}_teLossm.nc'.format(model_name)
            cdo.div(input='-fldmean {} -fldmean -div {} {} '
                    .format(toabloss_file, toabloss_file, teloss_file),
                    output=telossm_file)
            aux_baroceff_file = (diagworkdir +
                                 '/{}_aux_barocEff.nc'.format(model_name))
            cdo.sub(input='-reci {} -reci {}'
                    .format(telossm_file, tegainm_file),
                    output=aux_baroceff_file)
            baroceff_file = diagworkdir + '/{}_barocEff.nc'.format(model_name)
            cdo.div(input='{} -mulc,0.5 -add -reci {} -reci {}'
                    .format(aux_baroceff_file, tegainm_file, telossm_file),
                    output=baroceff_file)
            f_l = Dataset(baroceff_file)
            baroceff = f_l.variables['toab'][0, 0, 0]
            logger.info('Baroclinic efficiency (Lucarini et al., 2011): %s\n',
                        baroceff)
            baroc_eff_all[i_m] = baroceff
        else:
            pass
        # Water mass budget
        if wat in {'y', 'yes'}:
            logger.info('Computing water mass and latent energy budgets\n')
            wmassbudget_file = (diagworkdir +
                                '/{}_wmassBudget.nc'.format(model_name))
            removeif(wmassbudget_file)
            removeif(aux_file)
            cdo.sub(input="{} {}".format(evspsbl_file, pr_file),
                    output=aux_file)
            cdo.chname('hfls,wmass', input=aux_file, options='-b F32',
                       output=wmassbudget_file)
            wmass_gmean_file = (diagworkdir +
                                '/{}_wmass_gmean.nc'.format(model_name))
            cdo.fldmean(input='-yearmonmean {}'.format(wmassbudget_file),
                        output=wmass_gmean_file)
            f_l = Dataset(wmass_gmean_file)
            wmass_gmean_constant = f_l.variables['wmass'][:, :, :]
            wmass_all[i_m, 0] = np.nanmean(wmass_gmean_constant)
            wmass_all[i_m, 1] = np.nanstd(wmass_gmean_constant)
            logger.info('Water mass budget: %s\n', wmass_all[i_m, 0])
            # Latent energy budget
            latentenergy_file = (diagworkdir +
                                 '/{}_latentEnergy.nc'.format(model_name))
            removeif(latentenergy_file)
            removeif(aux_file)
            cdo.sub(input="{} -add -mulc,{} {} -mulc,{} {}"
                    .format(hfls_file, str(LC_SUB), prsn_file,
                            str(L_C), prr_file),
                    output=aux_file)
            cdo.chname('hfls,latent',
                       input=aux_file, options='-b F32',
                       output=latentenergy_file)
            laten_gmean_file = (diagworkdir +
                                '/{}_latenergy_gmean.nc'.format(model_name))
            cdo.fldmean(input='-yearmonmean {}'.format(latentenergy_file),
                        output=laten_gmean_file)
            f_l = Dataset(laten_gmean_file)
            laten_gmean_constant = f_l.variables['latent'][:, :, :]
            latent_all[i_m, 0] = np.nanmean(laten_gmean_constant)
            latent_all[i_m, 1] = np.nanstd(laten_gmean_constant)
            logger.info('Latent energy budget: %s\n', latent_all[i_m, 0])
            logger.info('Done\n')
        else:
            pass
        # Compute budgets over oceans and land separately.
        if lsm in {'y', 'yes'}:
            if e_b in {'y', 'yes'}:
                logger.info('Computing energy budgets over land and oceans\n')
                toab_ocean_file = (diagworkdir +
                                   '/{}_toab_ocean.nc'.format(model_name))
                cdo.mul(input='{} -eqc,0 {}'.format(toab_file, sftlf_fx),
                        output=toab_ocean_file)
                toab_oc_gmean_file = (diagworkdir +
                                      '/{}_toab_oc_gmean.nc'
                                      .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(toab_ocean_file),
                            output=toab_oc_gmean_file)
                f_l = Dataset(toab_oc_gmean_file)
                toab_oc_gmean_constant = f_l.variables['toab'][0, 0, 0]
                logger.info('TOA energy budget over oceans: %s\n',
                            toab_oc_gmean_constant)
                toab_oc_all[i_m] = toab_oc_gmean_constant
                toab_land_file = (diagworkdir +
                                  '/{}_toab_land.nc'.format(model_name))
                cdo.sub(input='{} {}'.format(toab_file, toab_ocean_file),
                        output=toab_land_file)
                cdo.setctomiss('0', input=toab_ocean_file, output=aux_file)
                move(aux_file, toab_ocean_file)
                cdo.setctomiss('0', input=toab_land_file, output=aux_file)
                move(aux_file, toab_land_file)
                toab_la_gmean_file = (diagworkdir +
                                      '/{}_toab_la_gmean.nc'
                                      .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(toab_land_file),
                            output=toab_la_gmean_file)
                f_l = Dataset(toab_la_gmean_file)
                toab_la_gmean_constant = f_l.variables['toab'][0, 0, 0]
                logger.info('TOA energy budget over land: %s\n',
                            toab_la_gmean_constant)
                toab_la_all[i_m] = toab_la_gmean_constant
                atmb_ocean_file = (diagworkdir +
                                   '/{}_atmb_ocean.nc'.format(model_name))
                cdo.mul(input='{} -eqc,0 {}'.format(atmb_file, sftlf_fx),
                        output=atmb_ocean_file)
                atmb_oc_gmean_file = (diagworkdir +
                                      '/{}_atmb_oc_gmean.nc'
                                      .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(atmb_ocean_file),
                            output=atmb_oc_gmean_file)
                f_l = Dataset(atmb_oc_gmean_file)
                atmb_oc_gmean_constant = f_l.variables['atmb'][0, 0, 0]
                logger.info('Atmospheric energy budget over oceans: %s\n',
                            atmb_oc_gmean_constant)
                atmb_oc_all[i_m] = atmb_oc_gmean_constant
                atmb_land_file = (diagworkdir +
                                  '/{}_atmb_land.nc'.format(model_name))
                cdo.sub(input='{} {}'.format(atmb_file, atmb_ocean_file),
                        output=atmb_land_file)
                aux_file = diagworkdir + '/{}_aux.nc'.format(model_name)
                cdo.setctomiss('0', input=atmb_ocean_file, output=aux_file)
                move(aux_file, atmb_ocean_file)
                cdo.setctomiss('0', input=atmb_land_file, output=aux_file)
                move(aux_file, atmb_land_file)
                atmb_la_gmean_file = (diagworkdir +
                                      '/{}_atmb_la_gmean.nc'
                                      .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(atmb_land_file),
                            output=atmb_la_gmean_file)
                f_l = Dataset(atmb_la_gmean_file)
                atmb_la_gmean_constant = f_l.variables['atmb'][0, 0, 0]
                logger.info('Atmospheric energy budget over land: %s\n',
                            atmb_la_gmean_constant)
                atmb_la_all[i_m] = atmb_la_gmean_constant
                surb_ocean_file = (diagworkdir +
                                   '/{}_surb_ocean.nc'.format(model_name))
                cdo.mul(input='{} -eqc,0 {}'.format(surb_file, sftlf_fx),
                        output=surb_ocean_file)
                surb_oc_gmean_file = (diagworkdir +
                                      '/{}_surb_oc_gmean.nc'
                                      .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(surb_ocean_file),
                            output=surb_oc_gmean_file)
                f_l = Dataset(surb_oc_gmean_file)
                surb_oc_gmean_constant = f_l.variables['surb'][0, 0, 0]
                logger.info('Surface energy budget over oceans: %s\n',
                            surb_oc_gmean_constant)
                surb_oc_all[i_m] = surb_oc_gmean_constant
                surb_land_file = (diagworkdir +
                                  '/{}_surb_land.nc'.format(model_name))
                cdo.sub(input='{} {}'.format(surb_file, surb_ocean_file),
                        output=surb_land_file)
                aux_file = diagworkdir + '/{}_aux.nc'.format(model_name)
                cdo.setctomiss('0', input=surb_ocean_file, output=aux_file)
                move(aux_file, surb_ocean_file)
                cdo.setctomiss('0', input=surb_land_file, output=aux_file)
                move(aux_file, surb_land_file)
                surb_la_gmean_file = (diagworkdir +
                                      '/{}_surb_la_gmean.nc'
                                      .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(surb_land_file),
                            output=surb_la_gmean_file)
                f_l = Dataset(surb_la_gmean_file)
                surb_la_gmean_constant = f_l.variables['surb'][0, 0, 0]
                logger.info('Surface energy budget over land: %s\n',
                            surb_la_gmean_constant)
                surb_la_all[i_m] = surb_la_gmean_constant
                logger.info('Done\n')
            else:
                pass
            if wat in {'y', 'yes'}:
                logger.info('Computing water mass and latent energy'
                            ' budgets over land and oceans\n')
                toab_ocean_file = (diagworkdir +
                                   '/{}_toab_ocean.nc'.format(model_name))
                wmassbudget_ocean_file = (diagworkdir +
                                          '/{}_wmassBudget_ocean.nc'
                                          .format(model_name))
                cdo.mul(input='{} -eqc,0 {}'
                        .format(wmassbudget_file, sftlf_fx),
                        output=wmassbudget_ocean_file)
                wmassbudget_land_file = (diagworkdir +
                                         '/{}_wmassBudget_land.nc'
                                         .format(model_name))
                cdo.sub(input='{} {}'
                        .format(wmassbudget_file, wmassbudget_ocean_file),
                        output=wmassbudget_land_file)
                cdo.setctomiss('0', input=wmassbudget_ocean_file,
                               output=aux_file)
                move(aux_file, wmassbudget_ocean_file)
                cdo.setctomiss('0', input=wmassbudget_land_file,
                               output=aux_file)
                move(aux_file, wmassbudget_land_file)
                wmass_oc_gmean_file = (diagworkdir +
                                       '/{}_wmass_oc_gmean.nc'
                                       .format(model_name))
                cdo.timmean(input='-fldmean {}'
                            .format(wmassbudget_ocean_file),
                            output=wmass_oc_gmean_file)
                f_l = Dataset(wmass_oc_gmean_file)
                wmass_oc_gmean_constant = f_l.variables['wmass'][0, 0, 0]
                logger.info('Water mass budget over oceans: %s\n',
                            wmass_oc_gmean_constant)
                wmass_oc_all[i_m] = wmass_oc_gmean_constant
                wmass_la_gmean_file = (diagworkdir +
                                       '/{}_wmass_la_gmean.nc'
                                       .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(wmassbudget_land_file),
                            output=wmass_la_gmean_file)
                f_l = Dataset(wmass_la_gmean_file)
                wmass_la_gmean_constant = f_l.variables['wmass'][0, 0, 0]
                logger.info('Water mass budget over land: %s\n',
                            wmass_la_gmean_constant)
                wmass_la_all[i_m] = wmass_la_gmean_constant
                latentenergy_ocean_file = (diagworkdir +
                                           '/{}_latentEnergy_ocean.nc'
                                           .format(model_name))
                cdo.mul(input='{} -eqc,0 {}'
                        .format(latentenergy_file, sftlf_fx),
                        output=latentenergy_ocean_file)
                latentenergy_land_file = (diagworkdir +
                                          '/{}_latentEnergy_land.nc'
                                          .format(model_name))
                cdo.sub(input='{} {}'
                        .format(latentenergy_file, latentenergy_ocean_file),
                        output=latentenergy_land_file)
                aux_file = diagworkdir + '/{}_aux.nc'.format(model_name)
                cdo.setctomiss('0', input=latentenergy_ocean_file,
                               output=aux_file)
                move(aux_file, latentenergy_ocean_file)
                cdo.setctomiss('0', input=latentenergy_land_file,
                               output=aux_file)
                move(aux_file, latentenergy_land_file)
                latent_oc_gmean_file = (diagworkdir +
                                        '/{}_latent_oc_gmean.nc'
                                        .format(model_name))
                cdo.timmean(input='-fldmean {}'
                            .format(latentenergy_ocean_file),
                            output=latent_oc_gmean_file)
                f_l = Dataset(latent_oc_gmean_file)
                latent_oc_gmean_constant = f_l.variables['latent'][0, 0, 0]
                logger.info('Latent energy budget over oceans: %s\n',
                            latent_oc_gmean_constant)
                latent_oc_all[i_m] = latent_oc_gmean_constant
                latent_la_gmean_file = (diagworkdir +
                                        '/{}_latent_la_gmean.nc'
                                        .format(model_name))
                cdo.timmean(input='-fldmean {}'.format(latentenergy_land_file),
                            output=latent_la_gmean_file)
                f_l = Dataset(latent_la_gmean_file)
                latent_la_gmean_constant = f_l.variables['latent'][0, 0, 0]
                logger.info('Latent energy budget over land: %s\n',
                            latent_la_gmean_constant)
                latent_la_all[i_m] = latent_la_gmean_constant
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        # Compute the Lorenz Energy Cycle.
        if lec in {'y', 'yes'}:
            logger.info('Computation of the Lorenz Energy'
                        'Cycle (year by year)\n')
            lecpath = os.path.join(plotpath2, 'LEC_results')
            if not os.path.exists(lecpath):
                os.makedirs(lecpath)
            maskorog = diagworkdir + '/orog.nc'
            cdo.setmisstoc('0', input='-setmisstoc,1 -sub {} {}'
                           .format(ua_file, ua_file), options='-b F32',
                           output=maskorog)
            ua_file_mask = diagworkdir + '/ua_fill.nc'
            removeif(ua_file_mask)
            cdo.add(input=('-setmisstoc,0 -selvar,ua {} '
                           '-setmisstoc,0 -mul {} -selvar,ua {}')
                    .format(ua_file, uas_file, maskorog), options='-b F32',
                    output=ua_file_mask)
            va_file_mask = diagworkdir + '/va_fill.nc'
            removeif(va_file_mask)
            cdo.add(input=('-setmisstoc,0 -selvar,va {} '
                           '-setmisstoc,0 -mul {} -selvar,ua {}')
                    .format(va_file, vas_file, maskorog), options='-b F32',
                    output=va_file_mask)
            energy3_file = diagworkdir + '/energy_short.nc'
            removeif(energy3_file)
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
                enfile_yr = diagworkdir + '/inputen.nc'
                tasfile_yr = diagworkdir + '/tas_yr.nc'
                tadiag_file = diagworkdir + '/ta_filled.nc'
                ncfile = diagworkdir + '/fourier_coeff.nc'
                cdo.selyear(y_r, input=energy3_file, options='-b F32',
                            output=enfile_yr)
                cdo.selyear(y_r, input=tas_file, options='-b F32',
                            output=tasfile_yr)
                fourc.fourier_coeff(tadiag_file, ncfile, enfile_yr, tasfile_yr)
                diagfile = (lecpath +
                            '/{}_{}_lec_diagram.png'.format(model_name, y_r))
                logfile = (lecpath +
                           '/{}_{}_lec_table.txt'.format(model_name, y_r))
                lect[y_i] = lorenz.lorenz(diagworkdir, model_name,
                                          y_r, ncfile, diagfile, logfile)
                y_i = y_i + 1
                removeif(ncfile)
                removeif(enfile_yr)
                removeif(tasfile_yr)
            lec_all[i_m, 0] = np.nanmean(lect)
            lec_all[i_m, 1] = np.nanstd(lect)
            logger.info('Intensity of the annual mean Lorenz Energy '
                        'Cycle: %s\n', lec_all[i_m, 0])
            logger.info('Done\n')
            os.remove(ua_file_mask)
            os.remove(va_file_mask)
            os.remove(energy3_file)
        else:
            pass
        # Compute the material entropy production
        if entr in {'y', 'yes'}:
            if met in {'1', '3'}:
                logger.info('Computation of the material entropy production '
                            'with the indirect method\n')
                # Horizonzal material entropy production
                horizentropy_file = (diagworkdir +
                                     '/{}_horizEntropy.nc'.format(model_name))
                removeif(horizentropy_file)
                removeif(aux_file)
                cdo.yearmonmean(input='-mulc,-1 -div -subc,{}  {}  {}'
                                .format(np.nanmean(toab_gmean_constant),
                                        toab_file, te_file), output=aux_file)
                cdo.chname('toab,shor', input=aux_file, options='-b F32',
                           output=horizentropy_file)
                horizentropy_mean_file = (diagworkdir +
                                          '/{}_horizEntropy_gmean.nc'
                                          .format(model_name))
                removeif(horizentropy_mean_file)
                cdo.fldmean(input=horizentropy_file,
                            options='-b F32', output=horizentropy_mean_file)
                f_l = Dataset(horizentropy_mean_file)
                horzentr_mean = f_l.variables['shor'][:, :, :]
                horzentr_all[i_m, 0] = np.nanmean(horzentr_mean)
                horzentr_all[i_m, 1] = np.nanstd(horzentr_mean)
                logger.info('Horizontal component of the material entropy '
                            'production: %s\n', horzentr_all[i_m, 0])
                # Vertical material entropy production
                verticalentropy_file = (diagworkdir +
                                        '/{}_verticalEntropy.nc'
                                        .format(model_name))
                cdo.yearmonmean(input=' -add {} -sub {} -add {} {}'
                                .format(rlds_file, rsds_file,
                                        rlus_file, rsus_file),
                                output=aux_file)
                cdo.mulc('-1', input='-mul  -sub -yearmonmean -reci {} \
                         -yearmonmean -reci {} {}'
                         .format(ts_file, te_file, aux_file),
                         output=verticalentropy_file)
                cdo.chname('ts,sver', input=verticalentropy_file,
                           options='-b F32', output=aux_file)
                move(aux_file, verticalentropy_file)
                verticalentropy_mean_file = (diagworkdir +
                                             '/{}_vertEntropy_gmean.nc'
                                             .format(model_name))
                removeif(verticalentropy_mean_file)
                cdo.fldmean(input=verticalentropy_file,
                            options='-b F32', output=verticalentropy_mean_file)
                cdo.chname('ts,sver', input=verticalentropy_mean_file,
                           options='-b F32', output=aux_file)
                move(aux_file, verticalentropy_mean_file)
                f_l = Dataset(verticalentropy_mean_file)
                vertentr_mean = f_l.variables['sver'][:, :, :]
                vertentr_all[i_m, 0] = np.nanmean(vertentr_mean)
                vertentr_all[i_m, 1] = np.nanstd(vertentr_mean)
                logger.info('Vertical component of the material entropy '
                            'production: %s\n', vertentr_all[i_m, 0])
                logger.info('Done\n')
            if met in {'2', '3'}:
                logger.info('Computation of the material entropy '
                            'production with the direct method\n')
                logger.info('1. Sensible heat fluxes\n')
                difftemp_file = (diagworkdir +
                                 '/{}_difftemp_bl.nc'.format(model_name))
                removeif(difftemp_file)
                cdo.sub(input='-reci {}  -reci {}'.format(tabl_file, ts_file),
                        options='-b F32', output=difftemp_file)
                sensentr_file = (diagworkdir +
                                 '/{}_sens_entr.nc'.format(model_name))
                removeif(sensentr_file)
                removeif(aux_file)
                cdo.timmean(input='-yearmonmean -monmean -mul {} {}'
                            .format(difftemp_file, hfss_file),
                            options='-b F32', output=aux_file)
                cdo.chname('tabl,ssens', input=aux_file, options='-b F32',
                           output=sensentr_file)
                sensentr_mean_file = (diagworkdir +
                                      '/{}_sensEntropy_gmean.nc'
                                      .format(model_name))
                removeif(sensentr_mean_file)
                cdo.fldmean(input=sensentr_file,
                            options='-b F32', output=sensentr_mean_file)
                f_l = Dataset(sensentr_mean_file)
                sensentr_mean = f_l.variables['ssens'][0, 0, 0]
                logger.info('Material entropy production associated with '
                            'sens. heat fluxes: %s\n', sensentr_mean)
                logger.info('Done\n')
                logger.info('2. Hydrological cycle\n')
                logger.info('2.1 Evaporation fluxes\n')
                evapentr_file = (diagworkdir +
                                 '/{}_evap_entr.nc'.format(model_name))
                removeif(evapentr_file)
                removeif(aux_file)
                cdo.timmean(input='-yearmonmean -monmean -div {} {}'
                            .format(hfls_file, ts_file), options='-b F32',
                            output=aux_file)
                cdo.chname('hfls,sevap', input=aux_file, options='-b F32',
                           output=evapentr_file)
                evapentr_mean_file = (diagworkdir +
                                      '/{}_evapEntropy_gmean.nc'
                                      .format(model_name))
                removeif(evapentr_mean_file)
                cdo.fldmean(input=evapentr_file, options='-b F32',
                            output=evapentr_mean_file)
                f_l = Dataset(evapentr_mean_file)
                evapentr_mean = f_l.variables['sevap'][0, 0, 0]
                logger.info('Material entropy production associated with '
                            'evaporation fluxes: %s\n', evapentr_mean)
                logger.info('Done\n')
                # Prepare masks for snowfall and rainfall
                maskrain_file = (diagworkdir +
                                 '/{}_maskprecr.nc'.format(model_name))
                removeif(maskrain_file)
                cdo.gtc('1.0E-7', input=prr_file, options=' -b F32',
                        output=maskrain_file)
                masksnow_file = (diagworkdir +
                                 '/{}_maskprecs.nc'.format(model_name))
                removeif(masksnow_file)
                cdo.gtc('1.0E-7', input=prsn_file, options=' -b F32',
                        output=masksnow_file)
                prrmask_file = (diagworkdir +
                                '/{}_prr_masked.nc'.format(model_name))
                removeif(prrmask_file)
                cdo.mul(input='{} {}'.format(maskrain_file, prr_file),
                        options='-b F32', output=prrmask_file)
                prsnmask_file = (diagworkdir +
                                 '/{}_prsn_masked.nc'.format(model_name))
                removeif(prsnmask_file)
                cdo.mul(input='{} {}'.format(masksnow_file, prsn_file),
                        options='-b F32', output=prsnmask_file)
                # Temperatures of the rainfall and snowfall clouds
                tliq_file = diagworkdir + '/{}_tliq.nc'.format(model_name)
                removeif(tliq_file)
                cdo.setrtomiss('-1000,0', input='-mul {} {}'
                               .format(tlcl_file, maskrain_file),
                               options='-b F32', output=tliq_file)
                tsol_file = diagworkdir + '/{}_tsol.nc'.format(model_name)
                removeif(tsol_file)
                cdo.setrtomiss('-1000,0', input='-mul {} {}'
                               .format(tlcl_file, masksnow_file),
                               options='-b F32', output=tsol_file)
                tdegl_file = diagworkdir + '/{}_tliqdeg.nc'.format(model_name)
                removeif(tdegl_file)
                cdo.subc('273.15', input=tliq_file, options='-b F32',
                         output=tdegl_file)
                tdegs_file = diagworkdir + '/{}_tsoldeg.nc'.format(model_name)
                removeif(tdegs_file)
                cdo.subc('273.15', input=tsol_file, options='-b F32',
                         output=tdegs_file)
                # Mask for ice cloud and original temperature
                # for phase changes from ice to rain
                maskice_file = (diagworkdir +
                                '/{}_maskice.nc'.format(model_name))
                removeif(maskice_file)
                cdo.ltc('0.0', input=tdegl_file, options='-b F32',
                        output=maskice_file)
                ticer_file = (diagworkdir +
                              '/{}_t_icerain_file'.format(model_name))
                removeif(ticer_file)
                cdo.setrtomiss('-1000,0', input='-mul {} {}'
                               .format(tliq_file, maskice_file),
                               options='-b F32', output=ticer_file)
                prrice_file = (diagworkdir +
                               '/{}_prr_ice_file.nc'.format(model_name))
                removeif(prrice_file)
                cdo.mul(input='{} {}'.format(maskice_file, prr_file),
                        options='-b F32', output=prrice_file)
                # Mask for water vapor cloud and original temperature
                # for phase changes from vapor to snow
                maskvap_file = (diagworkdir +
                                '/{}_maskvap.nc'.format(model_name))
                removeif(maskvap_file)
                cdo.gtc('0.0', input=tdegs_file, options='-b F32',
                        output=maskvap_file)
                tvaps_file = (diagworkdir +
                              '/{}_t_vapsnow.nc'.format(model_name))
                removeif(tvaps_file)
                cdo.setrtomiss('-1000,0', input='-mul {} {}'
                               .format(tsol_file, maskvap_file),
                               options='-b F32', output=tvaps_file)
                prsnvap_file = (diagworkdir +
                                '/{}_prsn_vap.nc'.format(model_name))
                removeif(prsnvap_file)
                cdo.mul(input='{} {}'.format(maskvap_file, prsn_file),
                        options='-b F32', output=prsnvap_file)
                logger.info('2.2 Rainfall precipitation\n')
                latrain_file = (diagworkdir +
                                '/{}_latentEnergy_rain.nc'.format(model_name))
                removeif(latrain_file)
                cdo.mulc(str(L_C), input=prrmask_file, options='-b F32',
                         output=latrain_file)
                rainentr_file = (diagworkdir +
                                 '/{}_rain_entr.nc'.format(model_name))
                removeif(rainentr_file)
                removeif(aux_file)
                cdo.timmean(input=('-yearmonmean -monmean -setmisstoc,0 '
                                   '-div {} {}').format(latrain_file,
                                                        tcloud_file),
                            options='-b F32', output=aux_file)
                cdo.chname('prr,srain', input=aux_file, options='-b F32',
                           output=rainentr_file)
                rainentr_mean_file = (diagworkdir +
                                      '/{}_rainEntropy_gmean.nc'
                                      .format(model_name))
                removeif(rainentr_mean_file)
                cdo.fldmean(input=rainentr_file,
                            options='-b F32', output=rainentr_mean_file)
                f_l = Dataset(rainentr_mean_file)
                rainentr_mean = f_l.variables['srain'][0, 0, 0]
                logger.info('Material entropy production associated with '
                            'rainfall: %s\n', rainentr_mean)
                logger.info('Done\n')
                logger.info('2.3 Snowfall precipitation\n')
                latsnow_file = (diagworkdir +
                                '/{}_latentEnergy_snow.nc'.format(model_name))
                removeif(latsnow_file)
                cdo.mulc(str(LC_SUB), input=prsnmask_file, options='-b F32',
                         output=latsnow_file)
                snowentr_file = (diagworkdir +
                                 '/{}_snow_entr.nc'.format(model_name))
                removeif(snowentr_file)
                removeif(aux_file)
                cdo.timmean(input=('-yearmonmean -monmean -setmisstoc,0 '
                                   '-div {} {}').format(latsnow_file,
                                                        tcloud_file),
                            options='-b F32', output=aux_file)
                cdo.chname('prsn,ssnow', input=aux_file, options='-b F32',
                           output=snowentr_file)
                snowentr_mean_file = (diagworkdir +
                                      '/{}_snowEntropy_gmean.nc'
                                      .format(model_name))
                removeif(snowentr_mean_file)
                cdo.fldmean(input=snowentr_file,
                            options='-b F32', output=snowentr_mean_file)
                f_l = Dataset(snowentr_mean_file)
                snowentr_mean = f_l.variables['ssnow'][0, 0, 0]
                logger.info('Material entropy production associated with '
                            'snowfall: %s\n', snowentr_mean)
                logger.info('Done\n')
                logger.info('2.4 Melting of snow at the surface \n')
                latmelt_file = (diagworkdir +
                                '/{}_latentEnergy_snowmelt.nc'
                                .format(model_name))
                removeif(latmelt_file)
                cdo.mulc(str(L_S),
                         input='-divc,{} {}'.format(str(LC_SUB), latsnow_file),
                         options='-b F32', output=latmelt_file)
                meltentr_file = (diagworkdir +
                                 '/{}_snowmelt_entr.nc'.format(model_name))
                removeif(meltentr_file)
                removeif(aux_file)
                cdo.timmean(input=('-yearmonmean -monmean -setmisstoc,0 '
                                   '-divc,273.15 {}').format(latmelt_file),
                            options='-b F32', output=aux_file)
                cdo.chname('prsn,smelt', input=aux_file, options='-b F32',
                           output=meltentr_file)
                meltentr_mean_file = (diagworkdir +
                                      '/{}_snowmeltEntropy_gmean.nc'
                                      .format(model_name))
                removeif(meltentr_mean_file)
                cdo.fldmean(input=meltentr_file,
                            options='-b F32', output=meltentr_mean_file)
                f_l = Dataset(meltentr_mean_file)
                meltentr_mean = f_l.variables['smelt'][0, 0, 0]
                logger.info('Material entropy production associated with snow '
                            'melting: %s\n', meltentr_mean)
                logger.info('Done\n')
                logger.info('2.5 Potential energy of the droplet\n')
                poten_file = (diagworkdir +
                              '/{}_potentialEnergy_droplet.nc'
                              .format(model_name))
                removeif(poten_file)
                cdo.mulc(GRAV, input='-mul {} -add {} {}'
                         .format(htop_file, prrmask_file, prsnmask_file),
                         options='-b F32', output=poten_file)
                potentr_file = (diagworkdir +
                                '/{}_pot_drop_entr.nc'.format(model_name))
                removeif(potentr_file)
                removeif(aux_file)
                cdo.timmean(input='-yearmonmean -monmean -div {} {}'
                            .format(poten_file, tcolumn_file),
                            options='-b F32', output=aux_file)
                cdo.chname('htop,spotp', input=aux_file, options='-b F32',
                           output=potentr_file)
                potentr_mean_file = (diagworkdir +
                                     '/{}_potentialdropEnergy_gmean.nc'
                                     .format(model_name))
                removeif(potentr_mean_file)
                cdo.fldmean(input=potentr_file,
                            options='-b F32', output=potentr_mean_file)
                f_i = Dataset(potentr_mean_file)
                potentr_mean = f_i.variables['spotp'][0, 0, 0]
                logger.info('Material entropy production associated with '
                            'potential energy of the droplet: %s\n',
                            potentr_mean)
                logger.info('Done\n')
                logger.info('3. Kinetic energy dissipation\n')
                if lec in {'y', 'yes'}:
                    cdo.yearmonmean(input=tasvert_file, output=aux_file)
                    f_l = Dataset(aux_file)
                    tabl_mean = f_l.variables['ts'][:, 0, 0]
                    minentr_mean = np.nanmean(lect / tabl_mean)
                    logger.info('Material entropy production associated with '
                                'kinetic energy dissipation: %s\n',
                                minentr_mean)
                else:
                    minentr_mean = 0.010
                    logger.info('I cannot compute the material entropy '
                                'production without the LEC...\n')
                    logger.info('I will assign a given value for the material '
                                'entropy production attributed to LEC '
                                '(0.01 W/m2*K)\n')
                sensentr_mean = masktonull(sensentr_mean)
                evapentr_mean = masktonull(evapentr_mean)
                rainentr_mean = masktonull(rainentr_mean)
                snowentr_mean = masktonull(snowentr_mean)
                meltentr_mean = masktonull(meltentr_mean)
                potentr_mean = masktonull(potentr_mean)
                minentr_mean = masktonull(minentr_mean)
                matentr = (float(sensentr_mean) - float(evapentr_mean) +
                           float(rainentr_mean) + float(snowentr_mean) +
                           float(potentr_mean) + float(minentr_mean) -
                           float(meltentr_mean))
                logger.info('Material entropy production with '
                            'the direct method: %s\n', matentr)
                matentr_all[i_m, 0] = matentr
                if met in {'3'}:
                    diffentr = (float(np.nanmean(vertentr_mean)) +
                                float(np.nanmean(horzentr_mean)) -
                                matentr)
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
            else:
                pass
        else:
            pass
        # Produce plots for the specific model
        logger.info('Running the plotting module for the budgets\n')
        if e_b in {'y', 'yes'}:
            plotsmod.balances(diagworkdir_up, plotpath2,
                              [toab_file, atmb_file, surb_file],
                              ['toab', 'atmb', 'surb'],
                              model_name)
            oname = '{}/{}_{}_timeser.png'.format(plotpath2,
                                                  model_name, 'toab')
            logger.info('Done\n')
        if wat in {'y', 'yes'}:
            logger.info('Running the plotting module for the water budgets\n')
            plotsmod.balances(diagworkdir_up, plotpath2,
                              [wmassbudget_file, latentenergy_file],
                              ['wmass', 'latent'],
                              model_name)
            logger.info('Done\n')
        else:
            pass
        if entr in {'y', 'yes'}:
            if met in {'1', '3'}:
                logger.info('Running the plotting module for the material '
                            'entropy production (indirect method)\n')
                plotsmod.entropy(plotpath2, verticalentropy_file, 'sver',
                                 'Vertical entropy production', model_name)
                logger.info('Done\n')
            elif met in {'2', '3'}:
                logger.info('Running the plotting module for the material '
                            'entropy production (direct method)\n')
                plotsmod.entropy(plotpath2, sensentr_file, 'ssens',
                                 'Sensible Heat entropy production',
                                 model_name)
                plotsmod.entropy(plotpath2, evapentr_file, 'sevap',
                                 'Evaporation entropy production',
                                 model_name)
                plotsmod.entropy(plotpath2, rainentr_file, 'srain',
                                 'Rainfall precipitation entropy production',
                                 model_name)
                plotsmod.entropy(plotpath2, snowentr_file, 'ssnow',
                                 'Snowfall precipitation entropy production',
                                 model_name)
                plotsmod.entropy(plotpath2, meltentr_file, 'smelt',
                                 'Snow melting entropy production', model_name)
                plotsmod.entropy(plotpath2, potentr_file, 'spotp',
                                 'Potential energy entropy production',
                                 model_name)
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        logger.info('Done for model: %s \n', model_name)
        i_m = i_m + 1
    # Produce multi-model ensemble plots (if more than one model are taken
    # into account).
    logger.info('I will now start multi-model plots')
    if e_b in {'y', 'yes'}:
        logger.info('Meridional heat transports\n')
        fig = plt.figure()
        fig.set_size_inches(12, 22)
        axi = plt.subplot(311)
        axi.set_figsize = (50, 50)
        for model_name in model_names:
            tot_transp_file = (diagworkdir_up +
                               '/total_transp_mean_{}.nc'.format(model_name))
            dataset = Dataset(tot_transp_file)
            name = 'total_{}'.format(model_name)
            toat = dataset.variables[name][:]
            lats = dataset.variables['lat'][:]
            plt.plot(np.array(lats), np.array(toat), color='black',
                     linewidth=1.)
        plt.title('(a) Total heat transports', fontsize=18)
        plt.xlabel('Latitude [deg]', fontsize=14)
        plt.ylabel('[W]', fontsize=14)
        plt.tight_layout()
        plt.ylim([-6.25E15, 6.25E15])
        plt.xlim(-90, 90)
        axi.tick_params(axis='both', which='major', labelsize=12)
        plt.grid()
        axi = plt.subplot(312)
        axi.set_figsize = (50, 50)
        for model_name in model_names:
            atm_transp_file = (diagworkdir_up +
                               '/atmos_transp_mean_{}.nc'.format(model_name))
            dataset = Dataset(atm_transp_file)
            name = 'atmos_{}'.format(model_name)
            atmt = dataset.variables[name][:]
            lats = dataset.variables['lat'][:]
            plt.plot(np.array(lats), np.array(atmt), color='black',
                     linewidth=1.)
        plt.title('(b) Atmospheric heat transports', fontsize=18)
        plt.xlabel('Latitude [deg]', fontsize=14)
        plt.ylabel('[W]', fontsize=14)
        plt.tight_layout()
        plt.ylim([-6.25E15, 6.25E15])
        plt.xlim(-90, 90)
        axi.tick_params(axis='both', which='major', labelsize=12)
        plt.grid()
        axi = plt.subplot(313)
        axi.set_figsize = (50, 50)
        for model_name in model_names:
            oce_transp_file = (diagworkdir_up +
                               '/ocean_transp_mean_{}.nc'.format(model_name))
            dataset = Dataset(oce_transp_file)
            name = 'ocean_{}'.format(model_name)
            surt = dataset.variables[name][:]
            lats = dataset.variables['lat'][:]
            plt.plot(np.array(lats), np.array(surt), color='black',
                     linewidth=1.)
        plt.title('(c) Oceanic heat transports', fontsize=18)
        plt.xlabel('Latitude [deg]', fontsize=14)
        plt.ylabel('[W]', fontsize=14)
        plt.tight_layout()
        plt.ylim([-3E15, 3E15])
        plt.xlim(-90, 90)
        axi.tick_params(axis='both', which='major', labelsize=12)
        plt.grid()
        oname = plotdir + '/meridional_transp.png'
        plt.savefig(oname)
        plt.close(fig)
    else:
        pass
    logger.info('Scatter plots')
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    colors = (0, 0, 0)
    axi = plt.subplot(321)
    axi.set_figsize = (50, 50)
    plt.scatter(toab_all[:, 0], atmb_all[:, 0], c=colors, alpha=1)
    plt.scatter(np.nanmean(toab_all[:, 0]), np.nanmean(atmb_all[:, 0]),
                c='red')
    s_l, interc, r_2, pval, stderr = stats.linregress(toab_all[:, 0],
                                                      atmb_all[:, 0])
    plotsmod.plot_ellipse(semimaj=np.nanstd(toab_all[:, 0]),
                          semimin=np.nanstd(atmb_all[:, 0]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(toab_all[:, 0]),
                          y_cent=np.nanmean(atmb_all[:, 0]), a_x=axi)
    plt.title('(a) TOA vs. atmospheric energy budget', fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('R_t [W m-2]', fontsize=14)
    plt.ylabel('F_a [W m-2]', fontsize=14)
    d_x = 0.01 * (max(toab_all[:, 0]) - min(toab_all[:, 0]))
    d_y = 0.01 * (max(atmb_all[:, 0]) - min(atmb_all[:, 0]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (toab_all[i_m, 0], atmb_all[i_m, 0]),
                     xytext=(toab_all[i_m, 0] + d_x, atmb_all[i_m, 0] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(322)
    axi.set_figsize = (50, 50)
    plt.scatter(baroc_eff_all, lec_all[:, 0], c=colors, alpha=1)
    plt.scatter(np.nanmean(baroc_eff_all), np.nanmean(lec_all[:, 0]), c='red')
    s_l, interc, r_2, pval, stderr = stats.linregress(baroc_eff_all,
                                                      lec_all[:, 0])
    plotsmod.plot_ellipse(semimin=np.nanstd(baroc_eff_all),
                          semimaj=np.nanstd(lec_all[:, 0]),
                          phi=np.arctan(s_l), x_cent=np.nanmean(baroc_eff_all),
                          y_cent=np.nanmean(lec_all[:, 0]), a_x=axi)
    plt.title('(b) Baroclinic efficiency vs. Intensity of LEC', fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('Eta', fontsize=14)
    plt.ylabel('W [W/m2]', fontsize=14)
    d_x = 0.01 * (max(baroc_eff_all) - min(baroc_eff_all))
    d_y = 0.01 * (max(lec_all[:, 0]) - min(lec_all[:, 0]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (baroc_eff_all[i_m], lec_all[i_m, 0]),
                     xytext=(baroc_eff_all[i_m] + d_x, lec_all[i_m, 0] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    axi = plt.subplot(323)
    axi.set_figsize = (50, 50)
    plt.scatter(horzentr_all[:, 0], vertentr_all[:, 0], c=colors, alpha=1)
    plt.scatter(np.nanmean(horzentr_all[:, 0]), np.nanmean(vertentr_all[:, 0]),
                c='red')
    s_l, interc, r_2, pval, stderr = stats.linregress(horzentr_all[:, 0],
                                                      vertentr_all[:, 0])
    plotsmod.plot_ellipse(semimin=np.nanstd(horzentr_all[:, 0]),
                          semimaj=np.nanstd(vertentr_all[:, 0]),
                          phi=np.arctan(s_l),
                          x_cent=np.nanmean(horzentr_all[:, 0]),
                          y_cent=np.nanmean(vertentr_all[:, 0]),
                          a_x=axi)
    xrang = abs(max(horzentr_all[:, 0]) - min(horzentr_all[:, 0]))
    yrang = abs(max(vertentr_all[:, 0]) - min(vertentr_all[:, 0]))
    plt.xlim(min(horzentr_all[:, 0]) - 0.1 * xrang,
             max(horzentr_all[:, 0]) + 0.1 * xrang)
    plt.ylim(min(vertentr_all[:, 0]) - 0.1 * yrang,
             max(vertentr_all[:, 0]) + 0.1 * yrang)
    x_x = np.linspace(min(horzentr_all[:, 0]) - 0.1 * xrang,
                      max(horzentr_all[:, 0]) + 0.1 * xrang,
                      10)
    y_y = np.linspace(min(vertentr_all[:, 0]) - 0.1 * yrang,
                      max(vertentr_all[:, 0]) + 0.1 * yrang,
                      10)
    x_m, y_m = np.meshgrid(x_x, y_y)
    z_m = x_m + y_m
    c_p = plt.contour(x_m, y_m, z_m, colors='black', linestyles='dashed',
                      linewidths=1.)
    plt.clabel(c_p, inline=True, inline_spacing=-4, fontsize=8)
    plt.title('(c) Vertical vs. horizontal component', fontsize=12)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('S_hor [W m-2 K-1]', fontsize=14)
    plt.ylabel('S_ver [W m-2 K-1]', fontsize=14)
    d_x = 0.01 * (max(horzentr_all[:, 0]) - min(horzentr_all[:, 0]))
    d_y = 0.01 * (max(vertentr_all[:, 0]) - min(vertentr_all[:, 0]))
    for i_m in np.arange(modnum):
        axi.annotate(str(i_m + 1), (horzentr_all[i_m, 0],
                                    vertentr_all[i_m, 0]),
                     xytext=(horzentr_all[i_m, 0] + d_x,
                             vertentr_all[i_m, 0] + d_y),
                     fontsize=12)
    axi.tick_params(axis='both', which='major', labelsize=12)
    plt.subplots_adjust(hspace=.3)
    plt.grid()
    indentr_all = horzentr_all[:, 0] + vertentr_all[:, 0]
    axi = plt.subplot(324)
    axi.set_figsize = (50, 50)
    plt.scatter(indentr_all, matentr_all[:, 0], c=colors, alpha=1)
    plt.scatter(np.nanmean(indentr_all), np.nanmean(matentr_all[:, 0]),
                c='red')
    s_l, interc, r_2, pval, stderr = stats.linregress(indentr_all,
                                                      matentr_all[:, 0])
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
    s_l, interc, r_2, pval, stderr = stats.linregress(te_all, indentr_all)
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
    s_l, interc, r_2, pval, stderr = stats.linregress(te_all, baroc_eff_all)
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
    oname = plotdir + '/scatters_summary.png'
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
    s_l, interc, r_2, pval, stderr = stats.linregress(toab_all[:, 0],
                                                      toab_all[:, 1])
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
    s_l, interc, r_2, pval, stderr = stats.linregress(atmb_all[:, 0],
                                                      atmb_all[:, 1])
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
    s_l, interc, r_2, pval, stderr = stats.linregress(surb_all[:, 0],
                                                      surb_all[:, 1])
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
    s_l, interc, r_2, pval, stderr = stats.linregress(surb_all[:, 0],
                                                      surb_all[:, 1])
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
    plt.savefig(plotdir + '/scatters_variability.png')
    plt.close(fig)
    logger.info("The diagnostic has finished. Now closing...\n")


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def masktonull(value):
    """Replace missing values with zeros."""
    try:
        value = float(value)
    except Warning:
        value = 0
    return value


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
