#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
;;#############################################################################
;; Diagnostic tool for several aspects of climate system thermodynamics
;;
;; Author
;; Valerio Lembo    
;;(Meteorological Institute, Hamburg University - valerio.lembo@uni-hamburg.de)
;;
;; Contributors
;; Frank   Lunkeit  
;; (Meteorological Insitute, Hamburg University - f.lunkeit@uni-hamburg.de)
;; Nikolay Koldunov 
;; (MARUM/AWI, nikolay.koldunov@awi.de, Germany)
;; 
;; Project
;; CRC - TRR 181 ”Energy transfers in Atmosphere and Ocean“
;;
;;#############################################################################
;; 
;; PREREQUISITES
;; 
;; The programme shares the same prerequisites with the overall ESMValTool 
;; architecture (see http://esmvaltool.readthedocs.io/en/latest/install.html)
;;
;; ADDITIONAL REQUIREMENTS
;; - NCO: NCO command operators are required for attributes and coordinates 
;;        manipulation on NetCDF files. The programme has been tested with
;;        the following  versions of NCO precompiled libraries:
;;        * Unix:   nco-4.6.7
;;        * Mac-OS: nco-4.6.6
;; - PyAstronomy: Python package PyAstronomy contains a function "zeocross1d"
;;                which is needed for the computation of transport maxima 
;;                locations. We tested the method with v. 0.12.0 of the 
;;                package, available in pip repositories.
;;
;; USAGE
;; 
;; 1: Obtain the datasets: the programme accepts the following variables as 
;;    input for the computations:
;;      Monthly mean resolution or higher:
;;      - TOA shortwave radiation downwards;
;;      - TOA shortwave radiation upwards;
;;      - TOA longwave radiation upwards (OLR);
;;      - Surface shortwave radiation downwards;
;;      - Surface shortwave radiation upwards;
;;      - Surface longwave radiation downwards;
;;      - Surface longwave radiation upwards;
;;      - Surface turbulent latent heat fluxes;
;;      - Surface turbulent sensible heat fluxes;
;;      - Surface temperature;
;;      - Specific humidity;
;;      - Near-surface (or 10m) zonal velocity;
;;      - Near-surface (or 10m) meridional velocity;
;;      Daily mean resolution or higher:
;;      - Near-surface temperature;
;;      - Air temperature (on pressure levels);
;;      - Horizontal velocity (on pressure levels);
;;      - Meridional velocity (on pressure levels);
;;      - Vertical velocity (on pressure levels);
;;    Data on lonlat grid are accepted, with CMOR-compliant coordinate system.
;;    The pre-processing modules of ESMValTool scheme will take care of 
;;    converting known grids and recognized datasets to CMOR standards. For a 
;;    a list of known formats, see 
;;       http://esmvaltool.readthedocs.io/en/latest/running.html#tab-obs-data
;;    Whether you need to use non-recognized datasets (provided that they are 
;;    globally gridded), ask for support! 
;;
;; 2: A configuration template is available in the ESMValTool release. Set your
;;   own paths to local directories here. Input datasets are read in MODELPATH,
;;   MODELPATH2, OBSPATH or OBSPATH2 output datasets are stored in WORKPATH, 
;;   plots in PLOTPATH (refer to the manual for ESMValTool).
;;
;; 3: Set the namelist with your own input datasets. Here you can also set the
;;    length of the dataset you want to subset, the time resolution, and 
;;    whether you are working on a MacOS or Linux machine.
;; 4: Run the tool by typing:
;;          python main.py nml/namelist_lembo17.xml
;; 5: After the data preprocessing, choose the options as required by the 
;;    program;
;;
;; OUTPUT
;;
;; In the output directory, NetCDF files containing annual mean fields and 
;; global mean time series for budgets and material entropy production time
;; series. 
;; Log files containing the values of the LEC components are stored in 
;; subdirectories, identified by the year for which they are computed. LEC 
;; quantities are stored in different files, one for each year.
;; - energy_gns_* contains global, NH and SH time mean fields as a function of 
;;   wavenumbers;
;; - energy_mean_* contains time mean (lon,wave,plev) fields;
;; - energy_ts3d_* contains (time,lon,wave,plev) fields; 
;; - <MODELNAME>_energy_*.ps is a diagram flux displaying annual mean 
;;   components of the LEC;
;;
;; In the plot directory figures are stored as PNG files in each model subdir. 
;; 
;; For the budgets (toab: TOA energy budget, atmb: atmospheric energy budget,
;; surb: surface energy budget, latent: latent energy budget, wmass: water mass
;; budget) there can be found:
;;  - climatological mean fields (*_climap.png);
;;  - time series evolution of global, NH and SH mean fields (*_timeser.png)
;;  - meridional transports (*_transp.pdf);
;;  - scatter plots of annual peak magnitudes vs. locations (*_scatpeak.png)
;;  - LEC diagrams for each year (*_lec_diagram.png)
;;
;; For the multi-model ensembles scatter plots of EB mean values vs. their
;; interannual variability are provided (scatters_variability.png), a summary
;; panel of the main thermodynamic quantities averaged for each model 
;; (scatters_summary.png), and the zonal mean meridional heat transports for
;; each model (meridional_transp.png).
;;
;; For the material entropy production (sver: vertical, shor: horizontal 
;; (indirect method), sevap: evaporation component, slatpr: ice-rain 
;; phase changes, slatps: vapor-snow phase changes, spotp: potential energy, 
;; srain: rainfall, ssens: sensible heat, ssnow: snowfall) climatological mean
;; fields are shown.
;;
;; In the working directory, a log file (lembo17_log.txt) is produced, 
;; containing the global mean values of all the quantities.
;; 
;; N.B.: multi-model ensembles and means are performed on the results of the
;;       analysis. No multi-model mean analysis is performed on input fields.
;;       It is thus a deliberate choice not to use the multi-model mean tools
;;       provided by the postprocessor.
;;
;; SOFTWARE TREE DESCRIPTION
;; 
;; The tool is divided in three modules; one for the 
;; computation of energy and water mass budgets (and related meridional 
;; transports), one for the Lorenz Energy Cycle (LEC), one for the material
;; entropy production.
;; 
;; Module Dependencies:
;;   - MODULE 1: Budgets and Transports --> none;
;;   - MODULE 2: Lorenz Energy Cycle    --> none;
;;   - MODULE 3: Material Entropy Production: 
;;                                  * Indirect method --> requires Budgets
;;                                  * Direct method   --> requires Budgets and
;;                                                        LEC
;; 
;; - MODULE 1
;; Earth's energy budgets from radiative and heat fluxes at Top-of-Atmosphere, 
;; at the surface and in the atmosphere (as a residual).
;; If required by the user, water mass and latent energy budgets are computed
;; from rainfall, snowfall and evaporation fluxes.
;; If required by the user, computations are also separately performed over 
;; land and oceans. A land-sea mask (either binary or percentual) has to be 
;; provided.
;; Meridional transports, magnitude and location of the peaks in each 
;; hemisphere (only for heat transports) are also computed.
;; The baroclinic efficiency is computed from TOA energy budgets, emission
;; temperature (in turn retrieved from OLR) and near-surface temperature.
;;
;; - MODULE 2
;;  
;; The Lorenz Energy Cycle (LEC) is computed in spectral components from near-
;; surface temperatures, temperatures and the three components of velocities 
;; over pressure levels.
;; The storage and conversion terms are directly computed, the sources and 
;; sinks are retrieved as residuals.
;; Components are grouped into a zonal mean, stationary and transient eddy
;; part.
;; 
;; - MODULE 3
;;
;; The material entropy production is computed using the indirect method, the 
;; direct method or both (following Lucarini et al., 2014). 
;; For the indirect method a vertical and a horizontal component are provided.
;; For the direct method, all components are combined, related to the 
;; hydrological cycle (attributable to evaporation, rainfall and snowfall 
;; precipitation, phase changes and potential energy of the droplet), to the 
;; sensible heat fluxes and to kinetic energy dissipation. For the latter the
;; LEC computation is required, given that the strength of the LEC can be
;; considered as equal to the kinetic energy dissipated to heating.
;;
;;
;;
;;
;;#############################################################################
;;
;; 20170803-A_lemb_va: Modified header with description and caveats
;; 20170629-A_kold_ni: Atmospheric budgets diagnostics written
;; 20180524-A_lemb_va: first complete working thermodynamics diagnostics
;;
;; 
;;#############################################################################
"""

#import ConfigParser
import os
from subprocess import *
from shutil import move
import commands
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
#import imp

#New packages for version 2.0 of ESMValTool
import logging
import esmvaltool.diag_scripts.shared as e 
from esmvaltool.diag_scripts.shared import(group_metadata, run_diagnostic,
					   select_metadata, sorted_metadata)
#import esmvaltool.diag_scripts.shared.names as n
#import cf_units
#import iris

from cdo import *
from netCDF4 import Dataset

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

#sys.path.append('./diag_scripts/aux/Thermodynamics/')
#Locally used classe
import srvfile_read as srv
import fluxogram as fluxogram
from mkthe import *
import lorenz_cycle as lecscript
import fourier_coefficients as fourc
#import plot_script as plotsmod
from plot_script import *

#import diag_scripts.aux.Thermodynamics.srvfile_read as srv
#import diag_scripts.aux.Thermodynamics.fluxogram as fluxogram
#import diag_scripts.aux.Thermodynamics.lorenz_cycle as lec

#import projects
#Some classes I need to debug the code"
#import time
#import sys
#Will be removed"

logger = logging.getLogger(os.path.basename(__file__))

sigmainv = 17636684.3034 	# inverse of the Stefan-Boltzmann constant
lc       =  2501000 		# latent heat of condensation
lcsub    =  2835000 		# latent heat of sublimation
ls       =   334000		# latent heat of solidification
grav     =        9.81		# gravity acceleration

#logger.info('Hello world')

#if __name__ == '__main__':
#	
#   with run_diagnostic() as config:
#       main(config)


def main(cfg):

    logger.info('Entering the diagnostic tool')
    
    #Load paths
    workdir = cfg['work_dir']
    plotdir = cfg['plot_dir']
    logger.info('Work directory: ',workdir)
    logger.info('Plots directory: ', plotdir)
    diagworkdir_up = workdir
    inputpath      = commands.getoutput("pwd")

    cdo = Cdo()
    mkthe = Mkthe()
    plotsmod = Plot_script() 

    data = e.Datasets(cfg)
    logger.debug(data)
    models = data.get_info_list('dataset')
    model_names = list(set(models))
    model_names.sort()
    logger.info(model_names)
    
    #Retrieve land-sea masks (TO BE CHECKED)
    #sftlf_fx = {}
#    for model in models:
#        sftlf_fx[model_name] = currProject.get_cf_lmaskfile(project_info, 
#                                                            model)
        
    varnames = data.get_info_list('short_name')
    currVars = list(set(varnames))
    logger.debug(currVars)

    eb=str(cfg['eb'])
    lsm=str(cfg['lsm'])
    wat=str(cfg['wat'])
    lec=str(cfg['lec'])
    entr=str(cfg['entr'])
    met=str(cfg['met'])
    logger.info(eb)
    
    #Initialize multi-model arrays
    modnum  = len(model_names)
    te_all  = np.zeros(modnum)
    toab_all= np.zeros([modnum, 2])
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
    barocEff_all = np.zeros(modnum)
    lec_all = np.zeros(modnum)
    horzentr_all = np.zeros([modnum, 2])
    vertentr_all = np.zeros([modnum, 2])
    matentr_all  = np.zeros([modnum, 2])
    irrevers_all = np.zeros(modnum)
    diffentr_all = np.zeros([modnum, 2])

    logger.info("Entering main loop\n")
    i=0
    for model_name in model_names:
        
        diagworkdir = os.path.join(diagworkdir_up, model_name)
        plotpath2   = os.path.join(plotdir, model_name)
        if not os.path.exists(diagworkdir):
            os.makedirs(diagworkdir)
        if not os.path.exists(plotpath2):
            os.makedirs(plotpath2)
        else:
            for root, dirs, files in os.walk(diagworkdir):
                for name in files:
                    file_path = os.path.join(diagworkdir, name)
                    #os.remove(file_path)
                for name in dirs:
                    file_path = os.path.join(diagworkdir, name)
                    #for root, dirs, files in os.walk(file_path):
                    #    for named in files:
                    #        filed_path = os.path.join(diagworkdir, name, named)
                            #os.remove(filed_path)
                    #os.rmdir(file_path)
            for root, dirs, files in os.walk(plotpath2):
                for name in files:
                    file_path = os.path.join(plotpath2, name)
                    #os.remove(file_path)
                for name in dirs:
                    file_path = os.path.join(plotpath2, name)
                    #os.rmdir(file_path)
        
        #Reading file names for the specific model
        filenames = data.get_info_list('filename',dataset = model_name)
        #logger.info(filenames)  
        
        print(model_name)
        logger.info('Processing model: {}\n'.format(model_name))

        hfls_file  = filenames[0]
        hfss_file  = filenames[1]
        hus_file   = filenames[2]
        pr_file    = filenames[3]
        prsn_file  = filenames[4]
        ps_file    = filenames[5]
        rlds_file  = filenames[6]
        rlus_file  = filenames[7]
        rlut_file  = filenames[8]
        rsds_file  = filenames[9]
        rsdt_file  = filenames[10]
        rsus_file  = filenames[11]
        rsut_file  = filenames[12]
#        ta_file    = filenames[13]
        tas_file   = filenames[13]
        ts_file    = filenames[14]
#        ua_file    = filenames[16]
        uas_file   = filenames[15]
#        va_file    = filenames[18]
        vas_file   = filenames[16]
#        wap_file   = filenames[20]

        aux_file = diagworkdir + '/{}_aux.nc'.format(model_name)
        cdo.selvar('tas', input = tas_file,output = aux_file)
        move(aux_file,tas_file)

        logger.info('Computing auxiliary variables\n')
        # emission temperature
        te_file = diagworkdir + '/{}_te.nc'.format(model_name)
        cdo.sqrt(input = "-sqrt -mulc,{} {}".format(sigmainv, rlut_file), 
                 output = te_file)
        te_ymm_file = diagworkdir + '/{}_te_ymm.nc'.format(model_name)
        cdo.yearmonmean(input=te_file, output=te_ymm_file)
        te_gmean_file = diagworkdir + '/{}_te_gmean.nc'.format(model_name)
        cdo.timmean(input='-fldmean {}'.format(te_ymm_file), 
                    output = te_gmean_file)
        fl = Dataset(te_gmean_file)
        te_gmean_constant = fl.variables['rlut'][0, 0, 0]
        logger.info('Global mean emission temperature: {}\n'
                  .format(te_gmean_constant))
        te_all[i] = te_gmean_constant
        # temperature of the atmosphere-surface interface
        tasvert_file = diagworkdir + '/{}_tvertavg.nc'.format(model_name)
        removeif(tasvert_file)
        cdo.fldmean(input = '-mulc,0.5 -add {} -selvar,tas {}'.format(ts_file,
                     tas_file), options = '-b F32', output = tasvert_file) 
        # evaporation
        if (wat in {'y','yes'} or met in {'2','3'}):
                    evspsbl_file = diagworkdir + '/{}_evspsbl.nc'.format(model_name)
                    cdo.divc(str(lc), input = "{}".format(hfls_file), 
                             output = evspsbl_file)
                    #Rainfall precipitation 
                    prr_file = diagworkdir + '/{}_prr.nc'.format(model_name)
                    cdo.sub(input = "{} {}".format(pr_file, prsn_file), output = aux_file)
                    cdo.chname('pr,prr',input = aux_file, output = prr_file)
                    logger.info('Done\n')
        else:
        	    pass
        if entr in {'y','yes'}:
                if met in {'2','3'}:
                    os.chdir(diagworkdir)
                    mkthe.mkthe_main(diagworkdir,ts_file,hus_file,tas_file,
                                      ps_file,uas_file,vas_file,hfss_file,
                                      te_file,model_name)
                    tlcl_file = diagworkdir + '/{}_tlcl.nc'.format(model_name)
                    cdo.setrtomiss('400,1e36',input = 'tlcl.nc ', output = tlcl_file) 
                    tabl_file = diagworkdir+'/{}_tabl.nc'.format(model_name)
                    cdo.setrtomiss('400,1e36',input = 'tabl.nc ', output = tabl_file) 
                    htop_file = diagworkdir + '/{}_htop.nc'.format(model_name)
                    cdo.setrtomiss('12000,1e36',input = 'htop.nc ', output = htop_file) 
                    #Working temperatures for the hydrological cycle
                    tcloud_file = diagworkdir + '/{}_tcloud.nc'.format(model_name)
                    removeif(tcloud_file)
                    cdo.mulc('0.5',input = '-add {} {}'.format(tlcl_file, te_file), 
				options = '-b F32', output = tcloud_file)
                    tcolumn_file = diagworkdir + '/{}_t_vertavg_pot.nc'.format(model_name)
                    removeif(tcolumn_file)
                    cdo.mulc('0.5', input = '-add {} {}'.
                             format(ts_file, tcloud_file),options = '-b F32', 
                             output = tcolumn_file)
                    #Working temperatures for the kin. en. diss. (updated)
                    tasvert_file = diagworkdir + '/{}_tboundlay.nc'.format(model_name)
                    removeif(tasvert_file)
                    cdo.fldmean(input='-mulc,0.5 -add {} {}'
                                .format(ts_file, tabl_file),options='-b F32',
                                output = tasvert_file)
                    os.chdir(inputpath)
                else:
                    pass
    
        #Energy budgets
        if eb in {'y','yes'}:
            logger.info('Computing energy budgets\n')
            #TOA energy budget
            toab_file = diagworkdir + '/{}_toab.nc'.format(model_name)
            removeif(toab_file)
            removeif(aux_file)
            cdo.sub(input = "-sub {} {} {}".format(rsdt_file, rsut_file, rlut_file),
                    output = aux_file)  
            cdo.chname('rsdt,toab', input = aux_file,options='-b F32',
                       output = toab_file)
            toab_gmean_file = diagworkdir + '/{}_toab_gmean.nc'.format(model_name)
            cdo.fldmean(input = '-yearmonmean {}'.format(toab_file), 
                        output = toab_gmean_file)        
            fl = Dataset(toab_gmean_file)
            toab_gmean_constant = fl.variables['toab'][:, :, :]
            toab_all[i,0] = np.nanmean(toab_gmean_constant)
            toab_all[i,1] = np.nanstd(toab_gmean_constant)
            logger.info('TOA energy budget: {}\n'.format(toab_all[i,0]))
            toab_ymm_file = diagworkdir + '/{}_toab_ymm.nc'.format(model_name)
            cdo.yearmonmean(input = toab_file, output=toab_ymm_file)
        	   #Surface energy budget 
            surb_file = diagworkdir + '/{}_surb.nc'.format(model_name)
            aux_surb_file = diagworkdir + '/{}_aux_surb.nc'.format(model_name)
            removeif(surb_file)
            removeif(aux_file)
            cdo.add(input = " {} {}".format(rsds_file, rlds_file), output = aux_surb_file) 
            cdo.sub(input = "-sub -sub -sub {} {} {} {} {}"
                    .format(aux_surb_file, rsus_file, rlus_file, hfls_file, hfss_file),
                    output = aux_file) 
            cdo.chname('rsds,surb',input = aux_file,options = '-b F32',
                       output = surb_file)
            surb_gmean_file = diagworkdir + '/{}_surb_gmean.nc'.format(model_name)
            cdo.fldmean(input = '-yearmonmean {}'.format(surb_file), 
                        output=surb_gmean_file)        
            fl = Dataset(surb_gmean_file)
            surb_gmean_constant = fl.variables['surb'][:, :, :]
            surb_all[i,0] = np.nanmean(surb_gmean_constant)
            surb_all[i,1] = np.nanstd(surb_gmean_constant)
            logger.info('Surface energy budget: {}\n'.format(surb_all[i,0]))
        	   #Atmospheric energy budget
            atmb_file = diagworkdir+'/{}_atmb.nc'.format(model_name)
            removeif(atmb_file)
            removeif(aux_file)
            cdo.sub(input = "{} {}".format(toab_file, surb_file), output = aux_file)  
            cdo.chname('toab,atmb',input = aux_file,options='-b F32',
                       output=atmb_file)
            atmb_gmean_file = diagworkdir + '/{}_atmb_gmean.nc'.format(model_name)
            cdo.fldmean(input = '-yearmonmean {}'.format(atmb_file), 
                        output=atmb_gmean_file)        
            fl = Dataset(atmb_gmean_file)
            atmb_gmean_constant = fl.variables['atmb'][:,:,:]
            atmb_all[i,0] = np.nanmean(atmb_gmean_constant)
            atmb_all[i,1] = np.nanstd(atmb_gmean_constant)
            logger.info('Atmospheric energy budget: {}\n'.format(atmb_all[i,0]))
            logger.info('Done\n')
            #Baroclinic efficiency
            maskGain_file = diagworkdir + '/{}_maskGain.nc'.format(model_name)
            cdo.gtc('0',input = toab_ymm_file, output = maskGain_file)
            maskLoss_file = diagworkdir+'/{}_maskLoss.nc'.format(model_name)
            cdo.ltc('0',input = toab_ymm_file, output = maskLoss_file)
            toabGain_file = diagworkdir+'/{}_toabGain.nc'.format(model_name)
            cdo.setrtomiss('-1000,0',input = '-mul {} {}'
                           .format(toab_ymm_file, maskGain_file), 
                           output = toabGain_file)
            toabLoss_file = diagworkdir+'/{}_toabLoss.nc'.format(model_name)
            cdo.setrtomiss('0,1000', input = '-mul {} {}'
                           .format(toab_ymm_file, maskLoss_file), 
                           output=toabLoss_file)
            teGain_file = diagworkdir + '/{}_teGain.nc'.format(model_name)
            cdo.setrtomiss('-1000,0',input = '-mul {} {}'
                           .format(te_ymm_file, maskGain_file),
                           output=teGain_file)
            teLoss_file = diagworkdir + '/{}_teLoss.nc'.format(model_name)
            cdo.setrtomiss('-1000,0',input = '-mul {} {}'
                           .format(te_ymm_file, maskLoss_file), 
                           output=teLoss_file)
            teGainm_file = diagworkdir+'/{}_teGainm.nc'.format(model_name)
            cdo.div(input = '-fldmean {} -fldmean -div {} {} '
                    .format( toabGain_file, toabGain_file, teGain_file), 
                    output = teGainm_file)
            teLossm_file = diagworkdir+'/{}_teLossm.nc'.format(model_name)
            cdo.div(input = '-fldmean {} -fldmean -div {} {} '
                    .format( toabLoss_file, toabLoss_file, teLoss_file), 
                    output = teLossm_file)
            aux_barocEff_file = diagworkdir+'/{}_aux_barocEff.nc'.format(model_name)
            cdo.sub(input = '-reci {} -reci {}'
                    .format(teLossm_file, teGainm_file), 
                    output = aux_barocEff_file)
            barocEff_file = diagworkdir+'/{}_barocEff.nc'.format(model_name)
            cdo.div(input = '{} -mulc,0.5 -add -reci {} -reci {}'
                    .format(aux_barocEff_file, teGainm_file, teLossm_file), 
                    output = barocEff_file)        
            fl = Dataset(barocEff_file)
            barocEff = fl.variables['toab'][0,0,0]
            logger.info('Baroclinic efficiency (Lucarini et al., 2011): {}\n'.format(barocEff))
            barocEff_all[i] = barocEff
        else:
        	    pass
     
        	#Water mass budget
        if wat in {'y','yes'}:
            logger.info('Computing water mass and latent energy budgets\n')
            wmassBudget_file = diagworkdir+'/{}_wmassBudget.nc'.format(model_name)
            removeif(wmassBudget_file)
            removeif(aux_file)
            cdo.sub(input = "{} {}".format(evspsbl_file, pr_file), output = aux_file)
            cdo.chname('hfls,wmass', input = aux_file,options='-b F32',
                       output=wmassBudget_file)
            wmass_gmean_file = diagworkdir + '/{}_wmass_gmean.nc'.format(model_name)
            cdo.fldmean(input = '-yearmonmean {}'.format(wmassBudget_file), 
                        output = wmass_gmean_file)        
            fl = Dataset(wmass_gmean_file)
            wmass_gmean_constant = fl.variables['wmass'][:, :, :]
            wmass_all[i,0] = np.nanmean(wmass_gmean_constant)
            wmass_all[i,1] = np.nanstd(wmass_gmean_constant)
            logger.info('Water mass budget: {}\n'.format(wmass_all[i,0]))
            #Latent energy budget
            latentEnergy_file = diagworkdir+'/{}_latentEnergy.nc'.format(model_name)
            removeif(latentEnergy_file)
            removeif(aux_file)
            cdo.sub(input = "{} -add -mulc,{} {} -mulc,{} {}"
                    .format(hfls_file, str(lcsub), prsn_file,  str(lc), prr_file),
                    output = aux_file)
            cdo.chname('hfls,latent', 
                       input = aux_file,options='-b F32',
                       output = latentEnergy_file)
            laten_gmean_file = diagworkdir+'/{}_latenergy_gmean.nc'.format(model_name)
            cdo.fldmean(input = '-yearmonmean {}'.format(latentEnergy_file), 
                        output = laten_gmean_file)        
            fl = Dataset(laten_gmean_file)
            laten_gmean_constant = fl.variables['latent'][:, :, :]
            latent_all[i,0] = np.nanmean(laten_gmean_constant)
            latent_all[i,1] = np.nanstd(laten_gmean_constant)
            logger.info('Latent energy budget: {}\n'.format(latent_all[i,0]))
            logger.info('Done\n')
        else:
        	    pass
     
        #######################
            
        if lsm in {'y','yes'}:	
            if eb in {'y','yes'}:
                logger.info('Computing energy budgets over land and oceans\n')
                toab_ocean_file = diagworkdir+'/{}_toab_ocean.nc'.format(model_name)
                cdo.mul(input='{} -eqc,0 {}'.format(toab_file, sftlf_fx[model_name]), 
                        output = toab_ocean_file)
                toab_oc_gmean_file = diagworkdir+'/{}_toab_oc_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(toab_ocean_file), 
                            output = toab_oc_gmean_file)        
                fl = Dataset(toab_oc_gmean_file)
                toab_oc_gmean_constant = fl.variables['toab'][0, 0, 0]
                logger.info('TOA energy budget over oceans: {}\n'
                          .format(toab_oc_gmean_constant))
                toab_oc_all[i] = toab_oc_gmean_constant
                toab_land_file = diagworkdir + '/{}_toab_land.nc'.format(model_name)
                cdo.sub(input = '{} {}'.format(toab_file, toab_ocean_file),
                        output = toab_land_file)
                cdo.setctomiss('0', input = toab_ocean_file, output = aux_file)
                move(aux_file, toab_ocean_file)
                cdo.setctomiss('0', input = toab_land_file, output = aux_file)
                move(aux_file, toab_land_file)
                toab_la_gmean_file = diagworkdir + '/{}_toab_la_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(toab_land_file), 
                            output = toab_la_gmean_file)        
                fl = Dataset(toab_la_gmean_file)
                toab_la_gmean_constant = fl.variables['toab'][0, 0, 0]
                logger.info('TOA energy budget over land: {}\n'.format(toab_la_gmean_constant))
                toab_la_all[i] = toab_la_gmean_constant
               
                atmb_ocean_file = diagworkdir+'/{}_atmb_ocean.nc'.format(model_name)
                cdo.mul(input = '{} -eqc,0 {}'.format(atmb_file, sftlf_fx[model_name]), 
                        output = atmb_ocean_file)
                atmb_oc_gmean_file = diagworkdir+'/{}_atmb_oc_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(atmb_ocean_file), 
                            output = atmb_oc_gmean_file)        
                fl = Dataset(atmb_oc_gmean_file)
                atmb_oc_gmean_constant = fl.variables['atmb'][0,0,0]
                logger.info('Atmospheric energy budget over oceans: {}\n'.format(atmb_oc_gmean_constant))
                atmb_oc_all[i] = atmb_oc_gmean_constant
                atmb_land_file = diagworkdir+'/{}_atmb_land.nc'.format(model_name)
                cdo.sub(input='{} {}'.format(atmb_file, atmb_ocean_file), 
                        output = atmb_land_file)
                aux_file = diagworkdir+'/{}_aux.nc'.format(model_name)
                cdo.setctomiss('0', input = atmb_ocean_file, output = aux_file)
                move(aux_file, atmb_ocean_file)
                cdo.setctomiss('0', input = atmb_land_file, output = aux_file)
                move(aux_file, atmb_land_file)
                atmb_la_gmean_file = diagworkdir + '/{}_atmb_la_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(atmb_land_file), 
                            output = atmb_la_gmean_file)        
                fl = Dataset(atmb_la_gmean_file)
                atmb_la_gmean_constant = fl.variables['atmb'][0,0,0]
                logger.info('Atmospheric energy budget over land: {}\n'.format(atmb_la_gmean_constant))
                atmb_la_all[i] = atmb_la_gmean_constant
                
                surb_ocean_file = diagworkdir+'/{}_surb_ocean.nc'.format(model_name)
                cdo.mul(input = '{} -eqc,0 {}'.format(surb_file, sftlf_fx[model_name]),
                        output = surb_ocean_file)
                surb_oc_gmean_file = diagworkdir+'/{}_surb_oc_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(surb_ocean_file), 
                            output = surb_oc_gmean_file)        
                fl = Dataset(surb_oc_gmean_file)
                surb_oc_gmean_constant = fl.variables['surb'][0,0,0]
                logger.info('Surface energy budget over oceans: {}\n'.format(surb_oc_gmean_constant))
                surb_oc_all[i] = surb_oc_gmean_constant
                surb_land_file = diagworkdir+'/{}_surb_land.nc'.format(model_name)
                cdo.sub(input = '{} {}'.format(surb_file, surb_ocean_file), 
                        output = surb_land_file)
                aux_file = diagworkdir+'/{}_aux.nc'.format(model_name)
                cdo.setctomiss('0', input = surb_ocean_file, output=aux_file)
                move(aux_file, surb_ocean_file)
                cdo.setctomiss('0', input = surb_land_file, output = aux_file)
                move(aux_file, surb_land_file)
                surb_la_gmean_file = diagworkdir+'/{}_surb_la_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(surb_land_file), 
                            output = surb_la_gmean_file)        
                fl = Dataset(surb_la_gmean_file)
                surb_la_gmean_constant = fl.variables['surb'][0, 0, 0]
                logger.info('Surface energy budget over land: {}\n'.format(surb_la_gmean_constant))
                surb_la_all[i] = surb_la_gmean_constant
                logger.info('Done\n')
            else:
                pass   
    	    
            if wat in {'y','yes'}:
                logger.info('Computing water mass and latent energy budgets over land and oceans\n')
                toab_ocean_file = diagworkdir+'/{}_toab_ocean.nc'.format(model_name)
                wmassBudget_ocean_file = diagworkdir+'/{}_wmassBudget_ocean.nc'.format(model_name)
                cdo.mul(input = '{} -eqc,0 {}'
                        .format(wmassBudget_file, sftlf_fx[model_name]), 
                        output = wmassBudget_ocean_file)
                wmassBudget_land_file = diagworkdir+'/{}_wmassBudget_land.nc'.format(model_name)
                cdo.sub(input = '{} {}'
                        .format(wmassBudget_file, wmassBudget_ocean_file), 
                        output = wmassBudget_land_file)
                cdo.setctomiss('0', input=wmassBudget_ocean_file, output = aux_file)
                move(aux_file, wmassBudget_ocean_file)
                cdo.setctomiss('0', input = wmassBudget_land_file, output = aux_file)
                move(aux_file, wmassBudget_land_file)
                wmass_oc_gmean_file = diagworkdir+'/{}_wmass_oc_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'
                            .format(wmassBudget_ocean_file), 
                            output = wmass_oc_gmean_file)        
                fl = Dataset(wmass_oc_gmean_file)
                wmass_oc_gmean_constant = fl.variables['wmass'][0, 0, 0]
                logger.info('Water mass budget over oceans: {}\n'.format(wmass_oc_gmean_constant))
                wmass_oc_all[i] = wmass_oc_gmean_constant
                wmass_la_gmean_file = diagworkdir+'/{}_wmass_la_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(wmassBudget_land_file), 
                            output=wmass_la_gmean_file)        
                fl = Dataset(wmass_la_gmean_file)
                wmass_la_gmean_constant = fl.variables['wmass'][0, 0, 0]
                logger.info('Water mass budget over land: {}\n'.format(wmass_la_gmean_constant))
                wmass_la_all[i] = wmass_la_gmean_constant
                
                latentEnergy_ocean_file = diagworkdir + '/{}_latentEnergy_ocean.nc'.format(model_name)
                cdo.mul(input = '{} -eqc,0 {}'
                        .format(latentEnergy_file, sftlf_fx[model_name]), 
                        output = latentEnergy_ocean_file)
                latentEnergy_land_file = diagworkdir+'/{}_latentEnergy_land.nc'.format(model_name)
                cdo.sub(input = '{} {}'
                        .format(latentEnergy_file, latentEnergy_ocean_file), 
                        output = latentEnergy_land_file)
                aux_file = diagworkdir+'/{}_aux.nc'.format(model_name)
                cdo.setctomiss('0', input = latentEnergy_ocean_file, 
                               output = aux_file)
                move(aux_file, latentEnergy_ocean_file)
                cdo.setctomiss('0', input = latentEnergy_land_file, 
                               output = aux_file )
                move(aux_file, latentEnergy_land_file)
                latent_oc_gmean_file = diagworkdir + '/{}_latent_oc_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(latentEnergy_ocean_file),
                            output = latent_oc_gmean_file)        
                fl = Dataset(latent_oc_gmean_file)
                latent_oc_gmean_constant = fl.variables['latent'][0,0,0]
                logger.info('Latent energy budget over oceans: {}\n'.format(latent_oc_gmean_constant))
                latent_oc_all[i] = latent_oc_gmean_constant
                latent_la_gmean_file = diagworkdir+'/{}_latent_la_gmean.nc'.format(model_name)
                cdo.timmean(input = '-fldmean {}'.format(latentEnergy_land_file),
                            output = latent_la_gmean_file)        
                fl = Dataset(latent_la_gmean_file)
                latent_la_gmean_constant = fl.variables['latent'][0, 0, 0]
                logger.info('Latent energy budget over land: {}\n'.format(latent_la_gmean_constant))
                latent_la_all[i] = latent_la_gmean_constant
                logger.info('Done\n')
            else:
                pass          
        else:
            pass
        ################################

        if lec in {'y','yes'}: 
            logger.info('Computation of the Lorenz Energy Cycle (year by year)\n')
            energy_file = diagworkdir + '/energy.nc'
            removeif(energy_file)
            cdo.merge(input = '{} {} {}'.format(ta_file, ua_file, va_file), 
                      options = '-b F32', output = energy_file)
            energy2_file = diagworkdir + '/energy_2.nc'
            removeif(energy2_file)
            cdo.merge(input = '{} {}'.format(energy_file, wap_file),
                      options = '-b F32', output = energy2_file)
            energy3_file = diagworkdir + '/energy_short.nc'
            removeif(energy3_file)
            cdo.setmisstoc('0',input = '-sellevel,10000/85000 {}'
                           .format(energy2_file),options = '-b F32', 
                           output = energy3_file)
            yrs = cdo.showyear(input = energy3_file)
            yrs = str(yrs)
            yrs2 = yrs.split()
            y = 0
            lect = np.zeros(len(yrs2))
            for yr in yrs2:
                os.chdir(auxpath)
                yr = int(filter(str.isdigit,yr))
                cdo.selyear(yr,input='-invertlat {}'.format(energy3_file), 
                            options = '-b F32', output='inputen.nc')
                fourc.fourier_coeff('inputen.nc')
                os.remove('inputen.nc')
                cdo.copy(input = 'fourier_coeff.nc', options = '-f srv', 
                         output = 'fourier_coeff.srv')
                srvfile=diagworkdir + '/{}_{}_fc.srv'.format(model_name,yr)
                ncfile=diagworkdir + '/{}_{}_fc.nc'.format(model_name,yr)
                move('fourier_coeff.srv', srvfile)
                move('fourier_coeff.nc', ncfile)
                
                diagfile = plotpath2+'/{}_{}_lec_diagram.png'.format(model_name, yr)
                logfile = plotpath2 + '/{}_{}_lec_table.txt'.format(model_name, yr)
                lect[y] = lecscript.lorenz(auxpath, diagworkdir, plotdir, model_name,
                                            yr, srvfile, ncfile, diagfile, logfile)
                caption = "Lorenz Energy Cycle for {}, year {}".format(model_name, yr)
                plot_id = "#lecdiag"
		dataIDs = "ta, ua, va, wap (time res: daily, vertical: pressure levels)"       
                #ESMValMD("both", diagfile, plot_tags, caption, plot_id, 
                #         dataIDs, diag_script, authors)
                y = y+1
            lec_all[i] = np.nanmean(lect)
            logger.info('Done\n')
            os.chdir(inputpath)
            os.remove(energy_file)
            os.remove(energy2_file)
            os.remove(energy3_file)
        else:
        	    pass
    
        if entr in {'y', 'yes'}:
            if met in {'1', '3'}:
                logger.info('Computation of the material entropy production with the indirect method\n')
                #Horizonzal material entropy production 
                horizEntropy_file = diagworkdir + '/{}_horizEntropy.nc'.format(model_name)
                removeif(horizEntropy_file)
                removeif(aux_file)
                cdo.yearmonmean(input = '-mulc,-1 -div -subc,{}  {}  {}'
                                .format(np.nanmean(toab_gmean_constant),
                                        toab_file, te_file), output = aux_file)
                cdo.chname('toab,shor',input = aux_file,options = '-b F32',
                           output = horizEntropy_file)
                horizEntropy_mean_file = diagworkdir + '/{}_horizEntropy_gmean.nc'.format(model_name)
                removeif(horizEntropy_mean_file)
                cdo.fldmean(input = horizEntropy_file,
                            options='-b F32', output = horizEntropy_mean_file)
                fl = Dataset(horizEntropy_mean_file)
                horzentr_mean = fl.variables['shor'][:, :, :]
                horzentr_all[i,0] = np.nanmean(horzentr_mean)
                horzentr_all[i,1] = np.nanstd(horzentr_mean)
                logger.info('Horizontal component of the material entropy production: {}\n'.format(horzentr_all[i,0]))
                #Vertical material entropy production
                #aux_verticalEntropy_file = diagworkdir+'/{}_aux_verticalEntropy.nc'.format(model_name)
                verticalEntropy_file = diagworkdir+'/{}_verticalEntropy.nc'.format(model_name)
                cdo.yearmonmean(input = ' -add {} -sub {} -add {} {}'
                                .format(rlds_file, rsds_file, rlus_file, rsus_file),
                                output = aux_file)
                cdo.mulc('-1', input = '-mul  -sub -yearmonmean -reci {} -yearmonmean -reci {} {}'
                         .format(ts_file, te_file, aux_file), output = verticalEntropy_file)
                cdo.chname('ts,sver',input = verticalEntropy_file,
                           options = '-b F32',output = aux_file)
                move(aux_file, verticalEntropy_file)	
                verticalEntropy_mean_file = diagworkdir+'/{}_vertEntropy_gmean.nc'.format(model_name)
                removeif(verticalEntropy_mean_file)
                cdo.fldmean(input = verticalEntropy_file,
                            options = '-b F32',output = verticalEntropy_mean_file)
                cdo.chname('ts,sver',input = verticalEntropy_mean_file,
                           options = '-b F32',output = aux_file)
                move(aux_file, verticalEntropy_mean_file)
                fl = Dataset(verticalEntropy_mean_file)
                vertentr_mean = fl.variables['sver'][:, :, :]
                vertentr_all[i,0] = np.nanmean(vertentr_mean)
                vertentr_all[i,1] = np.nanstd(vertentr_mean)
                logger.info('Vertical component of the material entropy production: {}\n'.format(vertentr_all[i,0]))
                logger.info('Done\n')
            if met in {'2','3'}:
                logger.info('Computation of the material entropy production with the direct method\n')
                logger.info('1. Sensible heat fluxes\n')
#                tasmean_file = diagworkdir + '/{}_tmean.nc'.format(model_name)
#                removeif(tasmean_file)
#                cdo.monmean(input = tas_file,options = '-b F32',
#                            output = tasmean_file)
                difftemp_file = diagworkdir+'/{}_difftemp_bl.nc'.format(model_name)
                removeif(difftemp_file)
                cdo.sub(input='-reci {}  -reci {}'.format(tabl_file, ts_file),
                        options = '-b F32',output = difftemp_file)
                sensentr_file = diagworkdir+'/{}_sens_entr.nc'.format(model_name)
                removeif(sensentr_file)
                removeif(aux_file)
                cdo.timmean(input = '-yearmonmean -monmean -mul {} {}'
                            .format(difftemp_file, hfss_file), options = '-b F32',
                            output = aux_file)
                cdo.chname('tabl,ssens', input = aux_file, options = '-b F32',
                           output = sensentr_file)
                sensentr_mean_file = diagworkdir+'/{}_sensEntropy_gmean.nc'.format(model_name)
                removeif(sensentr_mean_file)
                cdo.fldmean(input = sensentr_file,
                            options = '-b F32',output = sensentr_mean_file)
                fl = Dataset(sensentr_mean_file)
                sensentr_mean = fl.variables['ssens'][0, 0, 0]
                logger.info('Material entropy production associated with sens. heat fluxes: {}\n'.format(sensentr_mean))
                logger.info('Done\n')
                ######################################################################################################### 
                logger.info('2. Hydrological cycle\n')
                logger.info('2.1 Evaporation fluxes\n')
                evapentr_file = diagworkdir+'/{}_evap_entr.nc'.format(model_name)
                removeif(evapentr_file)
                removeif(aux_file)
                cdo.timmean(input='-yearmonmean -monmean -div {} {}'
                            .format(hfls_file, ts_file), options = '-b F32'
                            ,output = aux_file)
                cdo.chname('hfls,sevap', input = aux_file, options = '-b F32',
                           output = evapentr_file)
                evapentr_mean_file = diagworkdir+'/{}_evapEntropy_gmean.nc'.format(model_name)
                removeif(evapentr_mean_file)
                cdo.fldmean(input = evapentr_file,
                            options = '-b F32',output = evapentr_mean_file)
                fl = Dataset(evapentr_mean_file)
                evapentr_mean = fl.variables['sevap'][0, 0, 0]
                logger.info('Material entropy production associated with evaporation fluxes: {}\n'.format(evapentr_mean))
                logger.info('Done\n')
                ########################################################################################################## 
                #Masks for snowfall and rainfall
                maskrain_file = diagworkdir + '/{}_maskprecr.nc'.format(model_name)
                removeif(maskrain_file)
                cdo.gtc('1.0E-7',input = prr_file,options = ' -b F32', 
                        output = maskrain_file)
                masksnow_file = diagworkdir+'/{}_maskprecs.nc'.format(model_name)
                removeif(masksnow_file)
                cdo.gtc('1.0E-7', input = prsn_file,options = ' -b F32', 
                        output = masksnow_file)
                prrmask_file = diagworkdir + '/{}_prr_masked.nc'.format(model_name)
                removeif(prrmask_file)
                cdo.mul(input = '{} {}'.format(maskrain_file, prr_file), 
                        options = '-b F32', output = prrmask_file)
                prsnmask_file = diagworkdir +'/{}_prsn_masked.nc'.format(model_name)
                removeif(prsnmask_file)
                cdo.mul(input = '{} {}'.format(masksnow_file, prsn_file), 
                        options = '-b F32', output = prsnmask_file)
                #temperatures of the rainfall and snowfall clouds
                tliq_file = diagworkdir+'/{}_tliq.nc'.format(model_name)
                removeif(tliq_file)
                cdo.setrtomiss('-1000,0', input = '-mul {} {}'
                               .format(tlcl_file, maskrain_file),
                               options = '-b F32',output = tliq_file)
                tsol_file = diagworkdir+'/{}_tsol.nc'.format(model_name)
                removeif(tsol_file)
                cdo.setrtomiss('-1000,0',input = '-mul {} {}'
                               .format(tlcl_file, masksnow_file),
                               options = '-b F32',output = tsol_file)
                tdegl_file = diagworkdir+'/{}_tliqdeg.nc'.format(model_name)
                removeif(tdegl_file)
                cdo.subc('273.15',input = tliq_file, options = '-b F32', 
                         output = tdegl_file)
                tdegs_file = diagworkdir+'/{}_tsoldeg.nc'.format(model_name)
                removeif(tdegs_file)
                cdo.subc('273.15',input = tsol_file, options = '-b F32', 
                         output = tdegs_file)
                #Mask for ice cloud and original temperature for phase changes from ice to rain
                maskice_file = diagworkdir + '/{}_maskice.nc'.format(model_name)
                removeif(maskice_file)
                cdo.ltc('0.0',input = tdegl_file, options = '-b F32', 
                        output = maskice_file)
                ticer_file = diagworkdir + '/{}_t_icerain_file'.format(model_name)
                removeif(ticer_file)
                cdo.setrtomiss('-1000,0', input = '-mul {} {}'
                               .format(tliq_file,maskice_file), 
                               options = '-b F32', output = ticer_file)
                prrice_file = diagworkdir + '/{}_prr_ice_file.nc'.format(model_name)
                removeif(prrice_file)
                cdo.mul(input = '{} {}'.format(maskice_file, prr_file),
                        options = '-b F32', output = prrice_file)
                #Mask for water vapor cloud and original temperature for phase changes from vapor to snow
                maskvap_file = diagworkdir+'/{}_maskvap.nc'.format(model_name)
                removeif(maskvap_file)
                cdo.gtc('0.0',input = tdegs_file,options = '-b F32', 
                        output = maskvap_file)
                tvaps_file = diagworkdir + '/{}_t_vapsnow.nc'.format(model_name)
                removeif(tvaps_file)
                cdo.setrtomiss('-1000,0', input = '-mul {} {}'
                               .format(tsol_file, maskvap_file), 
                               options = '-b F32', output = tvaps_file)
                prsnvap_file = diagworkdir+'/{}_prsn_vap.nc'.format(model_name)
                removeif(prsnvap_file)
                cdo.mul(input = '{} {}'.format(maskvap_file, prsn_file), 
                        options = '-b F32', output = prsnvap_file)
                ##############################################################################################
                logger.info('2.2 Rainfall precipitation\n')
                latrain_file = diagworkdir + '/{}_latentEnergy_rain.nc'.format(model_name)
                removeif(latrain_file)
                cdo.mulc(str(lc), input = prrmask_file, options = '-b F32', 
                         output = latrain_file)
                rainentr_file = diagworkdir+'/{}_rain_entr.nc'.format(model_name)
                removeif(rainentr_file)
                removeif(aux_file)
                cdo.timmean(input = '-yearmonmean -monmean -setmisstoc,0 -div {} {}'
                            .format(latrain_file, tcloud_file), 
                            options = '-b F32', output = aux_file)
                cdo.chname('prr,srain',input = aux_file,options='-b F32',
                           output = rainentr_file)
                rainentr_mean_file = diagworkdir+'/{}_rainEntropy_gmean.nc'.format(model_name)
                removeif(rainentr_mean_file)
                cdo.fldmean(input = rainentr_file,
                            options = '-b F32',output = rainentr_mean_file)
                fl = Dataset(rainentr_mean_file)
                rainentr_mean = fl.variables['srain'][0, 0, 0]
                logger.info('Material entropy production associated with rainfall: {}\n'.format(rainentr_mean))
                logger.info('Done\n')
                ##############################################################################################
                logger.info('2.3 Snowfall precipitation\n')
                latsnow_file = diagworkdir+'/{}_latentEnergy_snow.nc'.format(model_name)
                removeif(latsnow_file)
                cdo.mulc(str(lcsub),input = prsnmask_file, options = '-b F32',
                         output = latsnow_file)
                snowentr_file = diagworkdir+'/{}_snow_entr.nc'.format(model_name)
                removeif(snowentr_file)
                removeif(aux_file)
                cdo.timmean(input = '-yearmonmean -monmean -setmisstoc,0 -div {} {}'
                            .format(latsnow_file, tcloud_file), 
                            options = '-b F32', output = aux_file)
                cdo.chname('prsn,ssnow', input = aux_file,options = '-b F32',
                           output = snowentr_file)
                snowentr_mean_file = diagworkdir+'/{}_snowEntropy_gmean.nc'.format(model_name)
                removeif(snowentr_mean_file)
                cdo.fldmean(input = snowentr_file, 
                            options = '-b F32', output = snowentr_mean_file)
                fl = Dataset(snowentr_mean_file)
                snowentr_mean = fl.variables['ssnow'][0, 0, 0]
                logger.info('Material entropy production associated with snowfall: {}\n'.format(snowentr_mean))
                logger.info('Done\n')
                ##############################################################################################
                logger.info('2.4 Phase changes from ice to rain\n')
                laticern_file = diagworkdir+'/{}_latentEnergy_phicerain.nc'.format(model_name)
                removeif(laticern_file)
                cdo.mulc(ls, input = prrice_file, options = '-b F32', 
                         output = laticern_file)
                phir_entr_file = diagworkdir+'/{}_ph_icerain_entr.nc'.format(model_name)
                removeif(phir_entr_file)
                removeif(aux_file)
                cdo.setmisstoc('0', input = '-div {} {}'.format(laticern_file,ticer_file), 
                               options = '-b F32', output = aux_file)
                cdo.chname('tlcl,slatpr', input = aux_file, options = '-b F32', 
                           output = phir_entr_file)
                phir_entr_mean_file = diagworkdir+'/{}_phir_Entropy_gmean.nc'.format(model_name)
                removeif(phir_entr_mean_file)
                cdo.timmean(input = '-fldmean {}'
                            .format(phir_entr_file), options = '-b F32', 
                            output = phir_entr_mean_file)
                fl = Dataset(phir_entr_mean_file)
                phir_entr_mean = fl.variables['slatpr'][0,0,0]
                logger.info('Material entropy production associated with phase changes from ice to rain: {}\n'.format(phir_entr_mean))
                logger.info('Done\n')
                ##############################################################################################
                logger.info('2.5 Phase changes from water vapor to snow\n')
                latvapsn_file = diagworkdir+'/{}_latentEnergy_phvapsnow.nc'.format(model_name)
                removeif(latvapsn_file)
                cdo.mulc(ls, input = prsnvap_file, options = '-b F32', 
                         output = latvapsn_file)
                phvs_entr_file = diagworkdir+'/{}_ph_vapsnow_entr.nc'.format(model_name)
                removeif(phvs_entr_file)
                removeif(aux_file)
                cdo.timmean(input = '-yearmonmean -monmean -setmisstoc,0 -div {} {}'
                            .format(latvapsn_file,tvaps_file), 
                            options = '-b F32', output = aux_file)
                cdo.chname('tlcl,slatps', input = aux_file, options = '-b F32', 
                           output = phvs_entr_file)
                phvs_entr_mean_file = diagworkdir+'/{}_phvs_Entropy_gmean.nc'.format(model_name)
                removeif(phvs_entr_mean_file)
                cdo.fldmean(input = phvs_entr_file, 
                            options='-b F32', output = phvs_entr_mean_file)
                fl = Dataset(phvs_entr_mean_file)
                phvs_entr_mean = fl.variables['slatps'][0, 0, 0]
                logger.info('Material entropy production associated with phase changes from vapor to snow: {}\n'.format(phvs_entr_mean))
                logger.info('Done\n')
                ##############################################################################################
                logger.info('2.6 Potential energy of the droplet\n')
                poten_file = diagworkdir+'/{}_potentialEnergy_droplet.nc'.format(model_name)
                removeif(poten_file)
                cdo.mulc(grav, input = '-mul {} -add {} {}'
                         .format(htop_file,prrmask_file,prsnmask_file), 
                         options = '-b F32', output = poten_file)
                potentr_file = diagworkdir+'/{}_pot_drop_entr.nc'.format(model_name)
                removeif(potentr_file)
                removeif(aux_file)
                cdo.timmean(input = '-yearmonmean -monmean -div {} {}'
                            .format(poten_file, tcolumn_file), options = '-b F32',
                            output = aux_file)
                cdo.chname('htop,spotp',input = aux_file, options = '-b F32', 
                           output = potentr_file)
                potentr_mean_file = diagworkdir+'/{}_potentialdropEnergy_gmean.nc'.format(model_name)
                removeif(potentr_mean_file)
                cdo.fldmean(input = potentr_file, 
                            options = '-b F32', output = potentr_mean_file)
                fi = Dataset(potentr_mean_file)
                potentr_mean = fi.variables['spotp'][0, 0, 0]
                logger.info('Material entropy production associated with potential energy of the droplet: {}\n'.format(potentr_mean))
                logger.info('Done\n')
                ##############################################################################################
                logger.info('3. Kinetic energy dissipation\n')
                if lec in {'y','yes'}: 
                    #minentr_file = diagworkdir+'/{}_kdiss_entr.nc'.format(model_name)
                    #removeif(minentr_file)
                    #cdo.timmean(input='-div {} {}'.format(kindiss_file,tasvert_file), options='-b F32', output=minentr_file)
                    cdo.yearmonmean(input = tasvert_file, output = aux_file)
                    fl = Dataset(aux_file)
                    tabl_mean = fl.variables['ts'][:, 0, 0]
                    minentr_mean = np.nanmean(lect / tabl_mean)
                    logger.info('Material entropy production associated with kinetic energy dissipation: {}\n'.format(minentr_mean))
                else:
                    minentr_mean = 0.010
                    logger.info('I cannot compute the material entropy production without the LEC...\n')
                    logger.info('I will assign a given value for the material entropy production attributed to LEC (0.01 W/m2*K)\n')
                ###########################################################################################
                sensentr_mean  = masktonull(sensentr_mean)
                evapentr_mean  = masktonull(evapentr_mean)
                rainentr_mean  = masktonull(rainentr_mean)
                snowentr_mean  = masktonull(snowentr_mean)
                phir_entr_mean = masktonull(phir_entr_mean)
                phvs_entr_mean = masktonull(phvs_entr_mean)
                potentr_mean   = masktonull(potentr_mean)
                minentr_mean   = masktonull(minentr_mean)
                matentr = float(sensentr_mean) - float(evapentr_mean) + float(rainentr_mean) + float(snowentr_mean)  - float(phir_entr_mean) + float(phvs_entr_mean) + float(potentr_mean) + float(minentr_mean)
                logger.info('Material entropy production with the direct method: {}\n'.format(matentr))
                matentr_all[i,0] = matentr
                if met in {'3'}:
                    diffentr = float(np.nanmean(vertentr_mean)) + float(np.nanmean(horzentr_mean)) - matentr
                    logger.info('Difference between the two methods: {}\n'.format(diffentr))
                    print('Difference of the two methods:', diffentr)
                    diffentr_all[i,0] = diffentr
                else:
                    pass
                irrevers=(matentr - float(minentr_mean)) / float(minentr_mean)
                logger.info('Degree of irreversibility of the system: {}\n'.format(irrevers))
                irrevers_all[i] = irrevers
            else:
                pass
     	else:
             pass
    
        	##########################################################################################################
        	#Graphical part starts here
        logger.info('Running the plotting module for the budgets\n')
        if eb in {'y','yes'}:
            plotsmod.balances(diagworkdir_up, plotpath2, 
                              [toab_file, atmb_file, surb_file],
                              ['toab','atmb','surb'],
                              model_name, lsm)
            oname = '{}/{}_{}_timeser.png'.format(plotpath2, model_name, 'toab')
            caption = "TOA EB annual mean time series (global; NH; SH)"
            plot_id = "#TOAEBtimeser"      
	    dataIDs = "rlut, rsdt, rsut (time res: monthly, vertical: 2D TOA)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_{}_timeser.png'.format(plotpath2, model_name, 'atmb')
	    dataIDs = "hfls, hfss, rlds, rlus, rlut, rsds, rsdt, rsus, rsut (time res: monthly, vertical: 2D TOA/surf)" 
            caption = "Atmospheric EB annual mean time series (global; NH; SH)"
            plot_id = "#AtmEBtimeser"       
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_{}_timeser.png'.format(plotpath2, model_name, 'surb')
	    dataIDs = "hfls, hfss, rlds, rlus, rsds, rsus (time res: monthly, vertical: 2D TOA/surf)" 
            caption = "Surface EB annual mean time series (global; NH; SH)"
            plot_id = "#SurEBtimeser"       
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_transp.png'.format(plotpath2, model_name)
            caption = "Annual mean northward heat transports (total,atmospheric,oceanic)"
            plot_id = "#transp"       
	    dataIDs = "hfls, hfss, rlds, rlus, rlut, rsds, rsdt, rsus, rsut (time res: monthly, vertical: 2D TOA/surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_energy_climap.png'.format(plotpath2, model_name)
            caption = "Climatological annual mean energy budgets (EB) (TOA,atmospheric,surface)"
            plot_id = "#enclimap"       
	    dataIDs = "hfls, hfss, rlds, rlus, rlut, rsds, rsdt, rsus, rsut (time res: monthly, vertical: 2D TOA/surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_scatpeak.png'.format(plotpath2,model_name)
            caption = "Meridional heat transports hemispheric peak magnitudes vs. positions (total,atmospheric,oceanic)"
            plot_id = "#scatpeak"       
	    dataIDs = "hfls, hfss, rlds, rlus, rlut, rsds, rsdt, rsus, rsut (time res: monthly, vertical: 2D TOA/surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            logger.info('Done\n')
        if wat in {'y','yes'}:
            logger.info('Running the plotting module for the water budgets\n')         
            plotsmod.balances(diagworkdir_up, plotpath2,
                              [wmassBudget_file, latentEnergy_file],
                              ['wmass','latent'], 
                              model_name, lsm)               
            oname = '{}/{}_{}_timeser.png'.format(plotpath2,model_name,'wmass')
            caption = "Water mass annual mean time series (global; NH; SH)"
            plot_id = "#wmasstimeser"       
	    dataIDs = "hfls, pr, prsn (time res: monthly, vertical: 2D surf)" 
            #ESMValMD("both",oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_{}_timeser.png'.format(plotpath2,model_name,'latent')
            caption = "Latent energy annual mean time series (global; NH; SH)"
            plot_id = "#latentimeser"       
	    dataIDs = "hfls, pr, prsn (time res: monthly, vertical: 2D surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_latent_transp.png'.format(plotpath2, model_name)
            caption = "Annual mean northward latent heat transports"
            plot_id = "#latentransp"       
	    dataIDs = "hfls, pr, prsn (time res: monthly, vertical: 2D surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_wmass_transp.png'.format(plotpath2, model_name)
            caption = "Annual mean northward water mass transports"
            plot_id = "#wmassransp"       
	    dataIDs = "hfls, pr, prsn (time res: monthly, vertical: 2D surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_wmass_climap.png'.format(plotpath2, model_name)
            caption = "Climatological annual mean water mass budgets"
            plot_id = "#wmclimap"       
	    dataIDs = "hfls, pr, prsn (time res: monthly, vertical: 2D surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            oname = '{}/{}_latent_climap.png'.format(plotpath2, model_name)
            caption = "Climatological annual mean latent energy budgets"
            plot_id = "#leclimap"       
	    dataIDs = "hfls, pr, prsn (time res: monthly, vertical: 2D surf)" 
            #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
            #         diag_script, authors)
            logger.info('Done\n')
        else:
            pass
        
        if entr in {'y','yes'}:
            if met in {'1'}:
                logger.info('Running the plotting module for the material entropy production (indirect method)\n')
                plotsmod.entropy(plotpath2, verticalEntropy_file, 'sver', 
                                 'Vertical entropy production', model_name)
                oname = '{}/{}_sver_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean vertical entropy production"
                plot_id = "#verclimap"       
	        dataIDs = "rlds, rlus, rlut, rsds, rsus, tas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, horizEntropy_file, 'shor', 
                                 'Horizontal entropy production', model_name)
                oname = '{}/{}_sver_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean horizontal entropy production"
                plot_id = "#horclimap"       
	        dataIDs = "rlut, rsdt rsut (time res: monthly, vertical: 2D TOA)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                logger.info('Done\n')
            elif met in {'2'}:
                logger.info('Running the plotting module for the material entropy production (direct method)\n')
                plotsmod.entropy(plotpath2, sensentr_file, 'ssens', 
                                 'Sensible Heat entropy production', model_name)
                oname = '{}/{}_ssens_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean entropy production associated with sensible heat fluxes"
                plot_id = "#ssensclimap"       
	        dataIDs = "hfss, hus, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, evapentr_file, 'sevap', 
                                 'Evaporation entropy production', model_name)
                oname = '{}/{}_sevap_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean entropy production associated with evaporation"
                plot_id = "#sevapclimap"      
	        dataIDs = "hfls, hfss, hus, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, rainentr_file, 'srain', 
                                 'Rainfall precipitation entropy production', model_name)
                oname = '{}/{}_srain_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean entropy production associated with rainfall"
                plot_id = "#srainclimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs,
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, snowentr_file, 'ssnow', 
                                 'Snowfall precipitation entropy production', model_name)
                oname = '{}/{}_ssnow_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean entropy production associated with snowfall"
                plot_id = "#snowclimap"       
	        dataIDs = "hfss, hus, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, phir_entr_file, 'slatpr', 
                                 'Phase changes ice -> rain entropy production', 
                                 model_name)    
                oname = '{}/{}_slatpr_climap.png'.format(plotpath2,model_name)
                caption = "Climatological annual mean entropy production associated with phase changes ice->rain"
                plot_id = "#slatprclimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, phvs_entr_file, 'slatps', 
                                 'Phase changes vapor -> snow entropy production', 
                                 model_name)
                oname = '{}/{}_slatps_climap.png'.format(plotpath2,model_name)
                caption = "Climatological annual mean entropy production associated with phase chanes vapour->snow"
                plot_id = "#slatpslimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, potentr_file, 'spotp', 
                                 'Potential energy entropy production', 
                                 model_name)
                oname = '{}/{}_spotp_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean entropy production associated with potential energy"
                plot_id = "#spotpclimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                logger.info('Done\n')
            elif met in {'3'}:
                logger.info('Running the plotting module for the material entropy production (indirect method)\n')
                plotsmod.entropy(plotpath2, verticalEntropy_file, 'sver', 
                                 'Vertical entropy production', model_name)
                oname = '{}/{}_sver_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean vertical entropy production"
                plot_id = "#verclimap"       
	        dataIDs = "rlds, rlus, rlut, rsds, rsus, tas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2,horizEntropy_file,'shor','Horizontal entropy production',model_name)
                oname = '{}/{}_sver_climap.png'.format(plotpath2, model_name)
                caption = "Climatological annual mean horizontal entropy production"
                plot_id = "#hprclimap"       
	        dataIDs = "rlut, rsdt rsut (time res: monthly, vertical: 2D TOA)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                logger.info('Done\n')
                logger.info('Running the plotting module for the material entropy production (direct method)\n')
                plotsmod.entropy(plotpath2, sensentr_file, 'ssens', 
                                 'Sensible Heat entropy production', model_name)
                oname = '{}/{}_ssens_climap.png'.format(plotpath2, model_name)
                caption = "Sensible heat fluxes"
                plot_id = "#ssensclimap"       
	        dataIDs = "hfss, hus, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, evapentr_file, 'sevap', 
                                 'Evaporation entropy production', model_name)
                oname = '{}/{}_sevap_climap.png'.format(plotpath2, model_name)
                caption = "Evaporation"
                plot_id = "#sevapclimap"      
	        dataIDs = "hfls, hfss, hus, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, rainentr_file, 'srain', 
                                 'Rainfall precipitation entropy production',
                                 model_name)
                oname = '{}/{}_srain_climap.png'.format(plotpath2, model_name)
                caption = "Rainfall precipitation"
                plot_id = "#srainclimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, snowentr_file, 'ssnow', 
                                 'Snowfall precipitation entropy production', 
                                 model_name)
                oname = '{}/{}_ssnow_climap.png'.format(plotpath2, model_name)
                caption = "Snowfall precipitation"
                plot_id = "#snowclimap"       
	        dataIDs = "hfss, hus, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, phir_entr_file, 'slatpr', 
                                 'Phase changes ice -> rain entropy production', 
                                 model_name)    
                oname = '{}/{}_slatpr_climap.png'.format(plotpath2, model_name)
                caption = "Phase changes ice->rain"
                plot_id = "#slatprclimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, phvs_entr_file, 'slatps', 
                                 'Phase changes vapor -> snow entropy production', 
                                 model_name)
                oname = '{}/{}_slatps_climap.png'.format(plotpath2, model_name)
                caption = "Phase chanes water vapour->snow"
                plot_id = "#slatpslimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                plotsmod.entropy(plotpath2, potentr_file, 'spotp', 
                                 'Potential energy entropy production', model_name)
                oname = '{}/{}_spotp_climap.png'.format(plotpath2, model_name)
                caption = "Potential energy of the droplet"
                plot_id = "#spotpclimap"       
	        dataIDs = "hfss, hus, pr, prsn, ps, rlut, tas, ts, uas, vas (time res: monthly, vertical: 2D TOA/surf)" 
                #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
                #         diag_script, authors)
                logger.info('Done\n')
            else:
                pass
        else:
            pass
        logger.info('Done for model: {} \n'.format(model_name))
        i = i+1
    
    logger.info('I will now start multi-model plots')
    if eb in {'y', 'yes'}:   
        aux_file        = diagworkdir_up + '/aux.nc'
        tot_transp_file = diagworkdir_up + '/total_transp_mean.nc'
        listf = filter(lambda x: x.startswith('total_transp_mean_'),
                      os.listdir(diagworkdir_up))
        if modnum>1: 
            cdo.merge(input = '{}/{} {}/{}'
                      .format(diagworkdir_up, str(listf[0]), diagworkdir_up,
                              str(listf[1])), output = tot_transp_file)
            for j in np.arange(2,modnum):
                removeif(aux_file)
                cdo.merge(input = '{} {}/{}'
                          .format(tot_transp_file, diagworkdir_up, listf[j]), 
                          options = '-b F32', output = aux_file)
                move(aux_file,tot_transp_file)
        else:
            inp=diagworkdir_up+'/'+listf[0]
            move(inp,tot_transp_file)
        atm_transp_file = diagworkdir_up+'/atmos_transp_mean.nc'.format(model_name)            
        listf = filter(lambda x: x.startswith('atmos_transp_mean_'),
                       os.listdir(diagworkdir_up))
        if modnum>1: 
            cdo.merge(input = '{}/{} {}/{}'
                      .format(diagworkdir_up, str(listf[0]), diagworkdir_up,
                              str(listf[1])), output = atm_transp_file)
            for j in np.arange(2, modnum):
                removeif(aux_file)
                cdo.merge(input = '{} {}/{}'
                          .format(atm_transp_file, diagworkdir_up, listf[j]),
                          options = '-b F32', output = aux_file)
                move(aux_file, atm_transp_file)
        else:
            inp=diagworkdir_up+'/'+listf[0]
            move(inp,atm_transp_file)
        oce_transp_file = diagworkdir_up + '/ocean_transp_mean.nc'.format(model_name)
        listf=filter(lambda x: x.startswith('ocean_transp_mean_'),
                     os.listdir(diagworkdir_up))
        if modnum>1: 
            cdo.merge(input = '{}/{} {}/{}'
                      .format(diagworkdir_up, str(listf[0]), diagworkdir_up,
                              str(listf[1])), output = oce_transp_file)
            for j in np.arange(2, modnum):
                removeif(aux_file)
                cdo.merge(input = '{} {}/{}'
                          .format(oce_transp_file, diagworkdir_up, listf[j]),
                          options = '-b F32', output = aux_file)
                move(aux_file,oce_transp_file)
        else:
            inp=diagworkdir_up+'/'+listf[0]
            move(inp,oce_transp_file)
        logger.info('Meridional heat transports')
        
        fig = plt.figure()
        fig.set_size_inches(12, 22)
        ax = plt.subplot(311)
        ax.set_figsize = (50,50)
        dataset = Dataset(tot_transp_file)
        toat = dataset.variables['total'][:]
        lats = dataset.variables['lat'][:]
        rank = toat.ndim
        if rank == 1:
            plt.plot(np.array(lats), np.array(toat), color = 'black', 
                     linewidth = 1.)
        else:
            for i in np.arange(modnum):
                plt.plot(np.array(lats), np.array(toat[i,:]),color = 'black', 
                         linewidth = 1.)
        #for j in np.arange(1, modnum-1):
        #    name='total_{}'.format(str(j + 1))
        #    toat =  dataset.variables[name][:]
        #    lats =  np.linspace(-90, 90, len(toat))
        #    plt.plot(lats, toat, color = 'black', linewidth = 1.)
        plt.title('Total heat transports', fontsize = 10)
        plt.xlabel('Latitude [deg]', fontsize = 10)
        plt.ylabel('[W]', fontsize = 10)
        plt.tight_layout()
        plt.ylim([-6.25E15, 6.25E15])
        plt.xlim(-90, 90)
        plt.grid()
        ax=plt.subplot(312)
        ax.set_figsize=(50, 50)
        dataset = Dataset(atm_transp_file)
        atmt=dataset.variables['atmos'][:]
        lats=dataset.variables['lat'][:]
        rank=toat.ndim
        if rank == 1:
            plt.plot(np.array(lats), np.array(atmt),color = 'black', 
                     linewidth = 1.)
        else:
            for i in np.arange(modnum):
                plt.plot(np.array(lats), np.array(atmt[i,:]), color = 'black', 
                         linewidth = 1.)
        #for j in np.arange(1, modnum-1):
        #    name='atmos_{}'.format(str(j + 1))
        #    atmt =  dataset.variables[name][:]
        #    lats =  np.linspace(-90, 90, len(atmt))
        #    plt.plot(lats, atmt, color = 'black', linewidth = 1.)
        plt.title('Atmospheric heat transports', fontsize = 10)
        plt.xlabel('Latitude [deg]', fontsize = 10)
        plt.ylabel('[W]', fontsize = 10)
        plt.tight_layout()
        plt.ylim([-6.25E15, 6.25E15])
        plt.xlim(-90, 90)
        plt.grid()
        ax=plt.subplot(313)
        ax.set_figsize=(50, 50)
        dataset = Dataset(oce_transp_file)
        surt=dataset.variables['ocean'][:]
        lats=dataset.variables['lat'][:]
        rank=toat.ndim
        if rank == 1:
            plt.plot(np.array(lats), np.array(surt),color = 'black', 
                     linewidth = 1.)
        else:
            for i in np.arange(modnum):
                plt.plot(np.array(lats), np.array(surt[i,:]), color = 'black', 
                         linewidth = 1.)
        #for j in np.arange(1, modnum-1):
        #    name = 'ocean_{}'.format(str(j + 1))
        #    surt =  dataset.variables[name][:]
        #    lats =  np.linspace(-90, 90, len(surt))
        #    plt.plot(lats, surt, color = 'black',linewidth = 1.)
        plt.title('Oceanic heat transports',fontsize = 10)
        plt.xlabel('Latitude [deg]', fontsize = 10)
        plt.ylabel('[W]', fontsize = 10)
        plt.tight_layout()
        plt.ylim([-3E15, 3E15])
        plt.xlim(-90, 90)
        plt.grid()
        oname=plotdir + 'meridional_transp.png'
        plt.savefig(oname)
        plt.show(fig)
        plt.close(fig)
        caption = "Annual mean meridional northward heat transports"
        plot_id = "#Mertransp"       
        dataIDs = "hfss, hfls, rlds, rlus, rlut, rsds, rsdt, rsus, rsut" 
#        ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, 
#                 diag_script, authors)
    else:
        pass
    
    logger.info('Scatter plots')

    fig=plt.figure()
    fig.set_size_inches(12, 22)
    colors = (0, 0, 0)
    ax=plt.subplot(321)
    ax.set_figsize=(50, 50)
    plt.scatter(toab_all[:,0], atmb_all[:,0], c = colors, alpha = 1)
    plt.scatter(np.nanmean(toab_all[:,0]), np.nanmean(atmb_all[:,0]), c = 'red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(toab_all[:, 0]),
                          semimin = np.nanstd(atmb_all[:, 0]),
                          phi = 0, x_cent = np.nanmean(toab_all[:, 0]),
                          y_cent=np.nanmean(atmb_all[:,0]),ax=ax)
    plt.title('TOA vs. atmospheric energy budget', fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('R_t [W m-2]',fontsize = 10)
    plt.ylabel('F_a [W m-2]',fontsize = 10)
    dx = 0.01 * (max(toab_all[:,0]) - min(toab_all[:,0]))
    dy = 0.01 * (max(atmb_all[:,0]) - min(atmb_all[:,0]))
    for i in np.arange(modnum):
        ax.annotate(str(i + 1), (toab_all[i,0], atmb_all[i,0]), 
                    xytext = (toab_all[i,0] + dx, atmb_all[i,0] + dy),
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    ax = plt.subplot(322)
    ax.set_figsize=(50, 50)
    plt.scatter(barocEff_all, irrevers_all, c = colors, alpha = 1)
    plt.scatter(np.nanmean(barocEff_all), np.nanmean(irrevers_all), c='red')
    plotsmod.plot_ellipse(semimaj = np.std(barocEff_all), 
                          semimin = np.std(irrevers_all),
                          phi = 0, x_cent = np.nanmean(barocEff_all),
                          y_cent = np.nanmean(irrevers_all), ax = ax)
    plt.title('Baroclinic efficiency vs. Irreversibility', fontsize = 10)
    rcParams['axes.titlepad'] = 1.
    rcParams['axes.labelpad'] = 1
    plt.xlabel('Eta', fontsize = 10)
    plt.ylabel('Alpha', fontsize = 10)
    dx = 0.01 * (max(barocEff_all) - min(barocEff_all))
    dy = 0.01 * (max(irrevers_all) - min(irrevers_all))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (barocEff_all[i], irrevers_all[i]), 
                    xytext = (barocEff_all[i] + dx, irrevers_all[i] + dy), 
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    ax = plt.subplot(323)
    ax.set_figsize = (50,50)
    plt.scatter(horzentr_all[:,0], vertentr_all[:,0], c = colors, alpha = 1)
    plt.scatter(np.nanmean(horzentr_all[:,0]), np.nanmean(vertentr_all[:,0]), 
                c = 'red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(horzentr_all[:,0]), 
                          semimin = np.nanstd(vertentr_all[:,0]),
                          phi = 0, x_cent = np.nanmean(horzentr_all[:,0]),
                          y_cent = np.nanmean(vertentr_all[:,0]), ax = ax)
    xrang=abs(max(horzentr_all[:,0])-min(horzentr_all[:,0]))
    yrang=abs(max(vertentr_all[:,0])-min(vertentr_all[:,0]))
    plt.xlim(min(horzentr_all[:,0])-0.1*xrang, max(horzentr_all[:,0])+0.1*xrang)
    plt.ylim(min(vertentr_all[:,0])-0.1*yrang, max(vertentr_all[:,0])+0.1*yrang)
    xx = np.linspace(min(horzentr_all[:,0])-0.1*xrang, max(horzentr_all[:,0])+0.1*xrang, 10)
    yy = np.linspace(min(vertentr_all[:,0])-0.1*yrang, max(vertentr_all[:,0])+0.1*yrang, 10)
    X, Y = np.meshgrid(xx, yy)
    Z = X + Y
    cp = plt.contour(X, Y, Z, colors = 'black', linestyles='dashed',
                     linewidths = 1.)    
    plt.clabel(cp, inline = True, inline_spacing = -4, fontsize = 8)
    plt.title('Vertical vs. horizontal material entropy production', fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('S_hor [W m-2 K-1]', fontsize = 10)
    plt.ylabel('S_ver [W m-2 K-1]', fontsize = 10)
    dx=0.01 * (max(horzentr_all[:,0]) - min(horzentr_all[:,0]))
    dy=0.01 * (max(vertentr_all[:,0]) - min(vertentr_all[:,0]))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (horzentr_all[i,0], vertentr_all[i,0]), 
                    xytext = (horzentr_all[i,0] + dx, vertentr_all[i,0] + dy),
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    indentr_all = horzentr_all[:,0] + vertentr_all[:,0]
    
    ax = plt.subplot(324)
    ax.set_figsize = (50,50)
    plt.scatter(indentr_all, diffentr_all[:,0], c = colors, alpha = 1)
    plt.scatter(np.nanmean(indentr_all), np.nanmean(diffentr_all[:,0]), c='red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(indentr_all), 
                          semimin = np.nanstd(diffentr_all[:,0]),
                          phi = 0, x_cent = np.nanmean(indentr_all),
                          y_cent = np.nanmean(diffentr_all[:,0]), ax = ax)
    plt.title('Indirect material entropy production vs. methods difference', 
              fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('S_mat [W m-2 K-1]', fontsize = 10)
    plt.ylabel('Delta S_mat [W m-2 K-1]', fontsize = 10)
    xrang=abs(max(indentr_all)-min(indentr_all))
    yrang=abs(max(diffentr_all[:,0])-min(diffentr_all[:,0]))
    plt.xlim(min(indentr_all)-0.1*xrang, max(indentr_all)+0.1*yrang)
    plt.ylim(min(diffentr_all[:,0])-0.1*yrang, max(diffentr_all[:,0])+0.1*yrang)
    dx = 0.01 * (max(indentr_all) - min(indentr_all))
    dy = 0.01 * (max(diffentr_all[:,0]) - min(diffentr_all[:,0]))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (indentr_all[i], diffentr_all[i,0]), 
                    xytext = (indentr_all[i] + dx, diffentr_all[i,0] + dy),
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    ax = plt.subplot(325)
    ax.set_figsize=(50, 50)
    plt.scatter(te_all, indentr_all, c = colors, alpha = 1)
    plt.scatter(np.nanmean(te_all), np.nanmean(indentr_all), c ='red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(te_all), 
                          semimin = np.nanstd(indentr_all),
                          phi = 0, x_cent = np.nanmean(te_all), 
                          y_cent = np.nanmean(indentr_all), ax = ax)
    plt.title('Indirect material entropy production vs. emission temperature', 
              fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('T_E [K]',fontsize = 10)
    plt.ylabel('S_mat [W m-2 K-1]', fontsize = 10)
    dx= 0.01 * (max(te_all) - min(te_all))
    dy= 0.01 * (max(indentr_all) - min(indentr_all))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (te_all[i], indentr_all[i]), 
                    xytext=(te_all[i] + dx, (indentr_all[i]) + dy), fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    ax=plt.subplot(326)
    ax.set_figsize=(50, 50)
    plt.scatter(te_all, barocEff_all, c = colors, alpha = 1)
    plt.scatter(np.nanmean(te_all), np.nanmean(barocEff_all), c = 'red')
    plotsmod.plot_ellipse(semimaj = np.std(te_all), 
                          semimin = np.std(barocEff_all),
                          phi = 0, x_cent = np.nanmean(te_all), 
                          y_cent = np.nanmean(barocEff_all), ax = ax)
    plt.title('Baroclinic efficiency vs. emission temperature',fontsize=10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('T_E [K]', fontsize = 10)
    plt.ylabel('Eta', fontsize = 10)
    dx = 0.01 * (max(te_all) - min(te_all))
    dy = 0.01 * (max(barocEff_all) - min(barocEff_all))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (te_all[i], barocEff_all[i]),
                    xytext = (te_all[i] + dx,barocEff_all[i] + dy),
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.grid()
    
    oname = plotdir + 'scatters_summary.png'
    plt.savefig(oname)
    plt.subplots_adjust(hspace = .3)
    caption = "Summary scatters from thermodynamics"
    plot_id = "#scatsum"
    dataIDs = "hfss, hfls, hus, pr, prsn, ps, rlds, rlus, rlut, rsds, rsdt, rsus, rsut, ta, tas, ts, ua, uas, va, vas, wap" 
#    ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, diag_script, authors)
    
    #Scatter plot of climatological mean values vs. interannual variability for each model
    logger.info('Scatter plots for inter-annual variability of some quantities')
    
    fig = plt.figure()
    fig.set_size_inches(12, 22)
    colors = (0, 0, 0)
    ax = plt.subplot(221)
    ax.set_figsize = (50, 50)
    plt.scatter(toab_all[:, 0], toab_all[:, 1], c = colors, alpha = 1)
    plt.scatter(np.nanmean(toab_all[:,0]), np.nanmean(toab_all[:,1]), c = 'red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(toab_all[:,0]), 
                          semimin = np.nanstd(toab_all[:,1]),
                          phi = 0, x_cent = np.nanmean(toab_all[:,0]), 
                          y_cent = np.nanmean(toab_all[:,1]), ax = ax)
    plt.title('TOA energy budget', fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('R_t [W m-2]', fontsize = 10)
    plt.ylabel('Sigma (R_t) [W m-2]', fontsize = 10)
    dx = 0.01 * (max(toab_all[:,0]) - min(toab_all[:,0]))
    dy = 0.01*(max(toab_all[:,1]) - min(toab_all[:,1]))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (toab_all[i,0], toab_all[i,1]), 
                    xytext = (toab_all[i,0] + dx, toab_all[i,1] + dy),
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    ax = plt.subplot(222)
    ax.set_figsize = (50,50)
    plt.scatter(atmb_all[:,0], atmb_all[:,1], c = colors, alpha = 1)
    plt.scatter(np.nanmean(atmb_all[:,0]), np.nanmean(atmb_all[:,1]), c = 'red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(atmb_all[:,0]), 
                          semimin = np.nanstd(atmb_all[:,1]),
                          phi = 0, x_cent = np.nanmean(atmb_all[:,0]), 
                          y_cent = np.nanmean(atmb_all[:,1]), ax = ax)
    plt.title('Atmospheric energy budget', fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('F_a [W m-2]', fontsize = 10)
    plt.ylabel('Sigma (F_a) [W m-2]', fontsize = 10)
    dx = 0.01 * (max(atmb_all[:,0]) - min(atmb_all[:,0]))
    dy = 0.01 * (max(atmb_all[:,1]) - min(atmb_all[:,1]))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (atmb_all[i,0], atmb_all[i,1]), 
                    xytext = (atmb_all[i,0] + dx,atmb_all[i,1] + dy),
                    fontsize = 10)
    ax.tick_params(axis='both', which='major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    ax = plt.subplot(223)
    ax.set_figsize = (50, 50)
    plt.scatter(surb_all[:,0], surb_all[:,1], c = colors, alpha = 1)
    plt.scatter(np.nanmean(surb_all[:,0]), np.nanmean(surb_all[:,1]), c = 'red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(surb_all[:,0]),
                          semimin = np.nanstd(surb_all[:,1]),
                          phi = 0,x_cent = np.nanmean(surb_all[:,0]),
                          y_cent = np.nanmean(surb_all[:,1]), ax = ax)
    plt.title('Surface energy budget', fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    plt.xlabel('F_s [W m-2]',fontsize = 10)
    plt.ylabel('Sigma (F_s) [W m-2]', fontsize = 10)
    dx = 0.01 * (max(surb_all[:,0]) - min(surb_all[:,0]))
    dy = 0.01 * (max(surb_all[:,1]) - min(surb_all[:,1]))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (surb_all[i,0], surb_all[i,1]), 
                    xytext=(surb_all[i,0] + dx,surb_all[i,1] + dy), fontsize=10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    
    entr_all = horzentr_all + vertentr_all
    
    ax = plt.subplot(224)
    ax.set_figsize = (50, 50)
    plt.scatter(entr_all[:,0], entr_all[:,1], c = colors, alpha = 1)
    plt.scatter(np.nanmean(entr_all[:,0]), np.nanmean(entr_all[:,1]), c = 'red')
    plotsmod.plot_ellipse(semimaj = np.nanstd(entr_all[:,0]),
                          semimin = np.nanstd(entr_all[:,1]),
                          phi = 0, x_cent = np.nanmean(entr_all[:,0]),
                          y_cent = np.nanmean(entr_all[:,1]), ax = ax)
    plt.title('Indirect material entropy production',fontsize = 10)
    rcParams['axes.titlepad'] = 1
    rcParams['axes.labelpad'] = 1
    xrang=abs(max(entr_all[:,0])-min(entr_all[:,0]))
    yrang=abs(max(entr_all[:,1])-min(entr_all[:,1]))
    plt.xlim(min(entr_all[:,0])-0.1*xrang, max(entr_all[:,0])+0.1*xrang)
    plt.ylim(min(entr_all[:,1])-0.1*yrang, max(entr_all[:,1])+0.1*yrang)
    plt.xlabel('S_i [W m-2 K-1]',fontsize = 10)
    plt.ylabel('Sigma (F_s) [W m-2 K-1]',fontsize = 10)
    dx = 0.01 * (max(entr_all[:,0]) - min(entr_all[:,0]))
    dy = 0.01 * (max(entr_all[:,1]) - min(entr_all[:,1]))
    for i in np.arange(modnum):
        ax.annotate(str(i+1), (entr_all[i,0], entr_all[i,1]), 
                    xytext = (entr_all[i,0] + dx, entr_all[i,1] + dy), 
                    fontsize = 10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 8)
    plt.subplots_adjust(hspace = .3)
    plt.grid()
    plt.savefig(plotdir + 'scatters_variability.png')
    plt.show(fig)
    plt.close(fig)
    caption = "Inter-annual variability vs. annual mean EBs and entropy"
    plot_id = "#scatsum"
    dataIDs = "hfss, hfls, rlds, rlus, rlut, rsds, rsdt, rsus, rsut, ta, tas, va, wap" 
    #ESMValMD("both", oname, plot_tags, caption, plot_id, dataIDs, diag_script, authors)
    
    logger.info("The diagnostic has finished. Now closing...\n")
        
#    e.write_references( diag_script,["A_lemb_va"],
#                                    ["A_kold_ni"], 	       # contributors
#                                    ["D0001","D0002"],      # diag_references
#                                    [""],                   # obs_references
#                                    ["P_trr181"],           # proj_references
#                                    project_info,
#                                    verbosity,
#                                    False)
    
    #write_references(diag_script, ["A_lemb_va","A_kold_ni"], ["D_0001", "D_0002"],"P_trr181")
    
    
    #return(project_info)


def removeif(filename):
    try:
        os.remove(filename)
    except OSError:
	pass

def masktonull(value):
    try:
    	value = float(value)
    except Warning:
	value = 0
    return value


if __name__ == '__main__':
	
   with run_diagnostic() as config:
       main(config)
