"""INTERNAL COMPUTATIONS.

Module containing all the core computations.

This module contains all the basic computations needed by the thermodynamics
diagnostic tool.

The functions that are here contained are:
- baroceff: function for the baroclinic efficiency;
- budgets: function for the energy budgets (TOA, atmospheric, surface);
- direntr: function for the material entropy production (direct method);
- entr: function for the computation of entropy as energy/temperature;
- evapentr: function for the evaporation related material entropy production;
- indentr: function for material entropy production (indirect method);
- kinentr: function for the kin. en. diss. related material entropy production;
- landoc_budg: function for budget computations over land and oceans;
- mask_precip: function for masking rainfall and snowfall regions;
- masktonull: function for masking nan values to null;
- meltentr: function for the entropy production from ground snow melting;
- potentr: function for the entropy production from pot. en. of the droplet;
- rainentr: function for the entropy production from rainfall precipitation;
- removeif: function for conditional file deleting;
- sensentr: function for the entropy production from sensible heat fluxes;
- snowentr: function for the entropy production from snowfall precipitation;
- wmbudg: function for water mass and latent energy budgets;
- write_eb: function for writing global mean energy budgets to file;

@author: valerio.lembo@uni-hamburg.de, Valerio Lembo, Hamburg University, 2019.
"""

import os
from shutil import move

import numpy as np
from cdo import Cdo
from netCDF4 import Dataset

import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.thermodyn_diagtool import mkthe

L_C = 2501000  # latent heat of condensation
LC_SUB = 2835000  # latent heat of sublimation
L_S = 334000  # latent heat of solidification
GRAV = 9.81  # gravity acceleration


def baroceff(model, wdir, aux_file, toab_file, te_file):
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
    """
    cdo = Cdo()
    removeif(aux_file)
    gain_file = wdir + '/{}_maskGain.nc'.format(model)
    cdo.gtc('0', input=toab_file, output=gain_file)
    loss_file = wdir + '/{}_maskLoss.nc'.format(model)
    cdo.ltc('0', input=toab_file, output=loss_file)
    toabgain_file = wdir + '/{}_toabGain.nc'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(toab_file, gain_file),
                   output=toabgain_file)
    toabloss_file = wdir + '/{}_toabLoss.nc'.format(model)
    cdo.setrtomiss('0,1000',
                   input='-mul {} {}'.format(toab_file, loss_file),
                   output=toabloss_file)
    tegain_file = wdir + '/{}_teGain.nc'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(te_file, gain_file),
                   output=tegain_file)
    teloss_file = wdir + '/{}_teLoss.nc'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(te_file, loss_file),
                   output=teloss_file)
    tegainm_file = wdir + '/{}_teGainm.nc'.format(model)
    cdo.div(input='-fldmean {0} -fldmean -div {0} {1} '.format(
        toabgain_file, tegain_file),
        output=tegainm_file)
    telossm_file = wdir + '/{}_teLossm.nc'.format(model)
    cdo.div(input='-fldmean {0} -fldmean -div {0} {1} '.format(
        toabloss_file, teloss_file),
        output=telossm_file)
    aux_baroceff_file = (wdir + '/{}_aux_barocEff.nc'.format(model))
    cdo.sub(input='-reci {} -reci {}'.format(telossm_file, tegainm_file),
            output=aux_baroceff_file)
    baroceff_file = wdir + '/{}_barocEff.nc'.format(model)
    cdo.div(input='{} -mulc,0.5 -add -reci {} -reci {}'.format(
        aux_baroceff_file, tegainm_file, telossm_file),
        output=baroceff_file)
    with Dataset(baroceff_file) as f_l:
        baroc = f_l.variables['toab'][0, 0, 0]
    remove_files = [
        gain_file, loss_file, toabgain_file, toabloss_file, tegain_file,
        teloss_file, tegainm_file, telossm_file, aux_baroceff_file
    ]
    for filen in remove_files:
        os.remove(filen)
    return baroc


def budgets(model, wdir, aux_file, input_data):
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
    """
    cdo = Cdo()
    hfls_file = e.select_metadata(input_data, short_name='hfls',
                                  dataset=model)[0]['filename']
    hfss_file = e.select_metadata(input_data, short_name='hfss',
                                  dataset=model)[0]['filename']
    rlds_file = e.select_metadata(input_data, short_name='rlds',
                                  dataset=model)[0]['filename']
    rlus_file = e.select_metadata(input_data, short_name='rlus',
                                  dataset=model)[0]['filename']
    rlut_file = e.select_metadata(input_data, short_name='rlut',
                                  dataset=model)[0]['filename']
    rsds_file = e.select_metadata(input_data, short_name='rsds',
                                  dataset=model)[0]['filename']
    rsus_file = e.select_metadata(input_data, short_name='rsus',
                                  dataset=model)[0]['filename']
    rsdt_file = e.select_metadata(input_data, short_name='rsdt',
                                  dataset=model)[0]['filename']
    rsut_file = e.select_metadata(input_data, short_name='rsut',
                                  dataset=model)[0]['filename']
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
    toab_ymm_file = wdir + '/{}_toab_ymm.nc'.format(model)
    cdo.yearmonmean(input=toab_file, output=toab_ymm_file)
    # Surface energy budget
    removeif(aux_file)
    cdo.add(input=" {} {}".format(rsds_file, rlds_file), output=aux_surb_file)
    cdo.sub(input="-sub -sub -sub {} {} {} {} {}".format(
        aux_surb_file, rsus_file, rlus_file, hfls_file, hfss_file),
        output=aux_file)
    surb_gmean = write_eb('rsds', 'surb', aux_file, surb_file, surb_gmean_file)
    # Atmospheric energy budget
    removeif(aux_file)
    cdo.sub(input="{} {}".format(toab_file, surb_file), output=aux_file)
    atmb_gmean = write_eb('toab', 'atmb', aux_file, atmb_file, atmb_gmean_file)
    eb_gmean = [toab_gmean, atmb_gmean, surb_gmean]
    eb_file = [toab_file, atmb_file, surb_file]
    # Delete files
    filenames = [
        aux_surb_file, toab_gmean_file, atmb_gmean_file, surb_gmean_file
    ]
    for filen in filenames:
        os.remove(filen)
    input_list = [
        hfls_file, hfss_file, rlds_file, rlus_file, rlut_file, rsds_file,
        rsdt_file, rsus_file, rsut_file
    ]
    return input_list, eb_gmean, eb_file, toab_ymm_file


def direntr(logger, model, wdir, input_data, aux_file, te_file, lect, flags):
    """Compute the material entropy production with the direct method.

    The function computes the material entropy production with the direct
    method, explicitly retrieving the components related to evaporation,
    rainfall and snowfall precipitation, snow melting at the ground, potential
    energy of the droplet, sensible heat fluxes and kinetic energy dissipation
    (from Lorenz Energy Cycle, LEC). The outputs are stored as NC files in
    terms of global mean time series, and in terms of annual mean
    (time,lat,lon) fields.

    Arguments:
    - logger: the log file where the global mean values are printed out;
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - filelist: the list containing all the input files;
    - aux_file: the name of a dummy aux. file to be used for computations;
    - lect: the annual mean value of the LEC strength;
    - lec: a flag having y (yes) value if the LEC is computed, n (no) if not.
    In the latter case, a reference value of 0.010 W*m-2*K-1 is given for the
    material entropy production related to the kinetic energy dissipation;
    - flags: a list of flags containing information on whether the water mass
    and energy budgets are computed, if the material entropy production has to
    be computed, if using the indirect, the direct method, or both methods;
    """
    lec = flags[1]
    aux_files = mkthe.init_mkthe_direntr(model, wdir, input_data, te_file,
                                         flags)
    htop_file = aux_files[1]
    prr_file = aux_files[2]
    tabl_file = aux_files[3]
    tasvert_file = aux_files[4]
    tcloud_file = aux_files[5]
    tcolumn_file = aux_files[6]
    tlcl_file = aux_files[7]
    hfls_file = e.select_metadata(input_data, short_name='hfls',
                                  dataset=model)[0]['filename']
    hfss_file = e.select_metadata(input_data, short_name='hfss',
                                  dataset=model)[0]['filename']
    prsn_file = e.select_metadata(input_data, short_name='prsn',
                                  dataset=model)[0]['filename']
    ts_file = e.select_metadata(input_data, short_name='ts',
                                dataset=model)[0]['filename']
    logger.info('Computation of the material entropy '
                'production with the direct method\n')
    logger.info('1. Sensible heat fluxes\n')
    infile_list = [hfss_file, tabl_file, ts_file]
    ssens, sensentr_file = sensentr(model, wdir, infile_list, aux_file)
    logger.info(
        'Material entropy production associated with '
        'sens. heat fluxes: %s\n', ssens)
    logger.info('2. Hydrological cycle\n')
    logger.info('2.1 Evaporation fluxes\n')
    infile_list = [hfls_file, ts_file]
    sevap, evapentr_file = evapentr(model, wdir, infile_list, aux_file)
    logger.info(
        'Material entropy production associated with '
        'evaporation fluxes: %s\n', sevap)
    infile_mask = [prr_file, prsn_file, tlcl_file]
    prrmask_file, prsnmask_file = mask_precip(model, wdir, infile_mask)
    logger.info('2.2 Rainfall precipitation\n')
    infile_rain = [prrmask_file, tcloud_file]
    srain, rainentr_file = rainentr(model, wdir, infile_rain, aux_file)
    logger.info(
        'Material entropy production associated with '
        'rainfall: %s\n', srain)
    logger.info('2.3 Snowfall precipitation\n')
    infile_snow = [prsnmask_file, tcloud_file]
    ssnow, latsnow_file, snowentr_file = snowentr(model, wdir, infile_snow,
                                                  aux_file)
    logger.info(
        'Material entropy production associated with '
        'snowfall: %s\n', ssnow)
    logger.info('2.4 Melting of snow at the surface \n')
    smelt, meltentr_file = meltentr(model, wdir, latsnow_file, aux_file)
    logger.info(
        'Material entropy production associated with snow '
        'melting: %s\n', smelt)
    logger.info('2.5 Potential energy of the droplet\n')
    infile_pot = [htop_file, prrmask_file, prsnmask_file, tcolumn_file]
    spot, potentr_file = potentr(model, wdir, infile_pot, aux_file)
    logger.info(
        'Material entropy production associated with '
        'potential energy of the droplet: %s\n', spot)
    os.remove(prrmask_file)
    os.remove(prsnmask_file)
    logger.info('3. Kinetic energy dissipation\n')
    skin = kinentr(logger, aux_file, tasvert_file, lect, lec)
    matentr = (float(ssens) - float(sevap) + float(srain) + float(ssnow) +
               float(spot) + float(skin) - float(smelt))
    logger.info('Material entropy production with '
                'the direct method: %s\n', matentr)
    irrevers = ((matentr - float(skin)) / float(skin))
    for filen in aux_files:
        os.remove(filen)
    entr_list = [
        sensentr_file, evapentr_file, rainentr_file, snowentr_file,
        meltentr_file, potentr_file
    ]
    return matentr, irrevers, entr_list


def entr(filelist, nin, nout, entr_file, entr_mean_file):
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
    """
    cdo = Cdo()
    en_file = filelist[0]
    tem_file = filelist[1]
    aux_file = filelist[2]
    removeif(aux_file)
    cdo.timmean(input='-yearmonmean -monmean -div {} {}'.format(
        en_file, tem_file),
        options='-b F32',
        output=aux_file)
    entr_gmean = write_eb(nin, nout, aux_file, entr_file, entr_mean_file)
    return entr_gmean


def evapentr(model, wdir, infile, aux_file):
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
    """
    evapentr_file = wdir + '/{}_evap_entr.nc'.format(model)
    evapentr_mean_file = wdir + '/{}_evapEntropy_gmean.nc'.format(model)
    flist = [infile[0], infile[1], aux_file]
    evapentr_gmean = entr(flist, 'hfls', 'sevap', evapentr_file,
                          evapentr_mean_file)
    evapentr_gmean = masktonull(evapentr_gmean)
    os.remove(evapentr_mean_file)
    return evapentr_gmean, evapentr_file


def indentr(model, wdir, infile, input_data, aux_file, toab_gmean):
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
    """
    cdo = Cdo()
    rlds_file = e.select_metadata(input_data, short_name='rlds',
                                  dataset=model)[0]['filename']
    rlus_file = e.select_metadata(input_data, short_name='rlus',
                                  dataset=model)[0]['filename']
    rsds_file = e.select_metadata(input_data, short_name='rsds',
                                  dataset=model)[0]['filename']
    rsus_file = e.select_metadata(input_data, short_name='rsus',
                                  dataset=model)[0]['filename']
    ts_file = e.select_metadata(input_data, short_name='ts',
                                dataset=model)[0]['filename']
    horzentropy_file = wdir + '/{}_horizEntropy.nc'.format(model)
    vertenergy_file = wdir + '/{}_verticalEnergy.nc'.format(model)
    vertentropy_file = wdir + '/{}_verticalEntropy.nc'.format(model)
    vertentropy_mean_file = wdir + '/{}_vertEntropy_gmean.nc'.format(model)
    horzentropy_mean_file = wdir + '/{}_horizEntropy_gmean.nc'.format(model)
    removeif(aux_file)
    cdo.yearmonmean(input='-mulc,-1 -div -subc,{}  {}  {}'.format(
        np.nanmean(toab_gmean), infile[1], infile[0]),
        output=aux_file)
    horzentr_mean = write_eb('toab', 'shor', aux_file, horzentropy_file,
                             horzentropy_mean_file)
    cdo.yearmonmean(input=' -add {} -sub {} -add {} {}'.format(
        rlds_file, rsds_file, rlus_file, rsus_file),
        output=vertenergy_file)
    cdo.mul(input='{} -sub -yearmonmean -reci {} -yearmonmean -reci {}'.format(
        vertenergy_file, infile[0], ts_file),
            output=aux_file)
    vertentr_mean = write_eb('rlds', 'sver', aux_file, vertentropy_file,
                             vertentropy_mean_file)
    remove_files = [
        horzentropy_mean_file, vertenergy_file, vertentropy_mean_file
    ]
    for filen in remove_files:
        os.remove(filen)
    return horzentr_mean, vertentr_mean, horzentropy_file, vertentropy_file


def kinentr(logger, aux_file, tasvert_file, lect, lec):
    """Compute the material entropy production from kin. energy dissipation.

    The function computes the material entropy production associated with the
    kinetic energy dissipation, through the intensity of the LEC.

    Arguments:
    - aux_file: the name of a dummy aux. file to be used for computations;
    - tasvert_file: a file containing the vertically integrated boundary layer
      temperature;
    - lect: an array containing the annual mean LEC intensity;
    - lec: a flag marking whether the LEC has been previously computed or not
    """
    cdo = Cdo()
    removeif(aux_file)
    if lec is True:
        cdo.yearmonmean(input=tasvert_file, output=aux_file)
        with Dataset(aux_file) as f_l:
            tabl_mean = f_l.variables['ts'][:, 0, 0]
        minentr_mean = np.nanmean(lect / tabl_mean)
        logger.info(
            'Material entropy production associated with '
            'kinetic energy dissipation: %s\n', minentr_mean)
        minentr_mean = masktonull(minentr_mean)
    else:
        minentr_mean = 0.010
        logger.info('I cannot compute the material entropy '
                    'production without the LEC...\n')
        logger.info('I will assign a given value for the material '
                    'entropy production attributed to LEC '
                    '(0.01 W/m2*K)\n')
    return minentr_mean


def landoc_budg(model, wdir, infile, mask, name):
    """Compute budgets separately on land and oceans.

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: the file containing the original budget field as (time,lat,lon);
    - mask: the file containing the land-sea mask;
    - name: the variable name as in the input file;
    """
    cdo = Cdo()
    ocean_file = wdir + '/{}_{}_ocean.nc'.format(model, name)
    oc_gmean_file = wdir + '/{}_{}_oc_gmean.nc'.format(model, name)
    land_file = wdir + '/{}_{}_land.nc'.format(model, name)
    la_gmean_file = wdir + '/{}_{}_la_gmean.nc'.format(model, name)
    aux_file = wdir + '/aux.nc'
    removeif(aux_file)
    cdo.mul(input='{} -eqc,0 {}'.format(infile, mask), output=ocean_file)
    cdo.timmean(input='-fldmean {}'.format(ocean_file), output=oc_gmean_file)
    with Dataset(oc_gmean_file) as f_l:
        oc_gmean = f_l.variables[name][0, 0, 0]
    cdo.sub(input='{} {}'.format(infile, ocean_file), output=land_file)
    cdo.setctomiss('0', input=ocean_file, output=aux_file)
    move(aux_file, ocean_file)
    cdo.setctomiss('0', input=land_file, output=aux_file)
    move(aux_file, land_file)
    cdo.timmean(input='-fldmean {}'.format(land_file), output=la_gmean_file)
    with Dataset(la_gmean_file) as f_l:
        la_gmean = f_l.variables[name][0, 0, 0]
    remove_files = [ocean_file, oc_gmean_file, land_file, la_gmean_file]
    for filen in remove_files:
        os.remove(filen)
    return oc_gmean, la_gmean


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
    cdo.mul(input='{} {}'.format(maskrain_file, prr_file),
            options='-b F32',
            output=prrmask_file)
    prsnmask_file = wdir + '/{}_prsn_masked.nc'.format(model)
    cdo.mul(input='{} {}'.format(masksnow_file, prsn_file),
            options='-b F32',
            output=prsnmask_file)
    # Temperatures of the rainfall and snowfall clouds
    tliq_file = wdir + '/{}_tliq.nc'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(tlcl_file, maskrain_file),
                   options='-b F32',
                   output=tliq_file)
    tsol_file = wdir + '/{}_tsol.nc'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(tlcl_file, masksnow_file),
                   options='-b F32',
                   output=tsol_file)
    tdegl_file = wdir + '/{}_tliqdeg.nc'.format(model)
    cdo.subc('273.15', input=tliq_file, options='-b F32', output=tdegl_file)
    tdegs_file = wdir + '/{}_tsoldeg.nc'.format(model)
    cdo.subc('273.15', input=tsol_file, options='-b F32', output=tdegs_file)
    # Mask for ice cloud and temperature for phase changes from ice to rain
    maskice_file = wdir + '/{}_maskice.nc'.format(model)
    cdo.ltc('0.0', input=tdegl_file, options='-b F32', output=maskice_file)
    ticer_file = wdir + '/{}_t_icerain_file'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(tliq_file, maskice_file),
                   options='-b F32',
                   output=ticer_file)
    prrice_file = wdir + '/{}_prr_ice_file.nc'.format(model)
    cdo.mul(input='{} {}'.format(maskice_file, prr_file),
            options='-b F32',
            output=prrice_file)
    # Mask for vapor cloud and temperature for phase changes from vapor to snow
    maskvap_file = wdir + '/{}_maskvap.nc'.format(model)
    cdo.gtc('0.0', input=tdegs_file, options='-b F32', output=maskvap_file)
    tvaps_file = wdir + '/{}_t_vapsnow.nc'.format(model)
    cdo.setrtomiss('-1000,0',
                   input='-mul {} {}'.format(tsol_file, maskvap_file),
                   options='-b F32',
                   output=tvaps_file)
    prsnvap_file = wdir + '/{}_prsn_vap.nc'.format(model)
    cdo.mul(input='{} {}'.format(maskvap_file, prsn_file),
            options='-b F32',
            output=prsnvap_file)
    remove_files = [
        maskrain_file, masksnow_file, tliq_file, tsol_file, tdegl_file,
        tdegs_file, maskice_file, ticer_file, prrice_file, maskvap_file,
        tvaps_file, prsnvap_file
    ]
    for filen in remove_files:
        os.remove(filen)
    return prrmask_file, prsnmask_file


def masktonull(value):
    """Replace missing values with zeros."""
    try:
        value = float(value)
    except Warning:
        value = 0
    return value


def meltentr(model, wdir, latsnow_file, aux_file):
    """Compute entropy production related to snow melting at the ground.

    The function computes the material entropy production related to snow
    melting at the ground, as part of the material entropy production
    obtained with the direct method (after Lucarini et al., 2011).

    Arguments:
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: the latent energy associated with snowfall precipitation;
    - aux_file: the name of a dummy aux. file to be used for computations;
    """
    cdo = Cdo()
    removeif(aux_file)
    latmelt_file = (wdir + '/{}_latentEnergy_snowmelt.nc'.format(model))
    meltentr_file = (wdir + '/{}_snowmelt_entr.nc'.format(model))
    meltentr_mean_file = wdir + '/{}_snowmeltEntropy_gmean.nc'.format(model)
    cdo.mulc(str(L_S),
             input='-divc,{} {}'.format(str(LC_SUB), latsnow_file),
             options='-b F32',
             output=latmelt_file)
    cdo.timmean(
        input='-yearmonmean -monmean -setmisstoc,0 -divc,273.15 {}'.format(
            latmelt_file),
        options='-b F32',
        output=aux_file)
    cdo.chname('prsn,smelt',
               input=aux_file,
               options='-b F32',
               output=meltentr_file)
    cdo.fldmean(input=meltentr_file,
                options='-b F32',
                output=meltentr_mean_file)
    with Dataset(meltentr_mean_file) as f_l:
        meltentr_gmean = f_l.variables['smelt'][0, 0, 0]
    meltentr_gmean = masktonull(meltentr_gmean)
    remove_files = [latmelt_file, meltentr_mean_file]
    for filen in remove_files:
        os.remove(filen)
    os.remove(latsnow_file)
    return meltentr_gmean, meltentr_file


def potentr(model, wdir, infile, aux_file):
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
    """
    cdo = Cdo()
    removeif(aux_file)
    htop_file = infile[0]
    prrmask_file = infile[1]
    prsnmask_file = infile[2]
    tcolumn_file = infile[3]
    poten_file = wdir + '/{}_potEnergy_drop.nc'.format(model)
    potentr_file = wdir + '/{}_pot_drop_entr.nc'.format(model)
    potentr_mean_file = wdir + '/{}_potEnergy_drop_gmean.nc'.format(model)
    cdo.mulc(GRAV,
             input='-mul {} -add {} {}'.format(htop_file, prrmask_file,
                                               prsnmask_file),
             options='-b F32',
             output=poten_file)
    flist = [poten_file, tcolumn_file, aux_file]
    potentr_gmean = entr(flist, 'htop', 'spotp', potentr_file,
                         potentr_mean_file)
    potentr_gmean = masktonull(potentr_gmean)
    remove_files = [poten_file, potentr_mean_file]
    for filen in remove_files:
        os.remove(filen)
    return potentr_gmean, potentr_file


def rainentr(model, wdir, infile, aux_file):
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
    """
    cdo = Cdo()
    prrmask_file = infile[0]
    removeif(aux_file)
    latrain_file = wdir + '/{}_latentEnergy_rain.nc'.format(model)
    rainentr_file = wdir + '/{}_rain_entr.nc'.format(model)
    rainentr_mean_file = wdir + '/{}_rainEntropy_gmean.nc'.format(model)
    cdo.mulc(str(L_C),
             input='-setmisstoc,0 {}'.format(prrmask_file),
             options='-b F32',
             output=latrain_file)
    flist = [latrain_file, infile[1], aux_file]
    rainentr_gmean = entr(flist, 'prr', 'srain', rainentr_file,
                          rainentr_mean_file)
    rainentr_gmean = masktonull(rainentr_gmean)
    remove_files = [latrain_file, rainentr_mean_file]
    for filen in remove_files:
        os.remove(filen)
    return rainentr_gmean, rainentr_file


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def sensentr(model, wdir, infile, aux_file):
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
    """
    cdo = Cdo()
    difftemp_file = wdir + '/{}_difftemp_bl.nc'.format(model)
    sensentr_file = (wdir + '/{}_sens_entr.nc'.format(model))
    sensentr_mean_file = wdir + '/{}_sensEntropy_gmean.nc'.format(model)
    cdo.reci(input='-sub -reci {}  -reci {}'.format(infile[1], infile[2]),
             options='-b F32',
             output=difftemp_file)
    flist = [infile[0], difftemp_file, aux_file]
    sensentr_gmean = entr(flist, 'hfss', 'ssens', sensentr_file,
                          sensentr_mean_file)
    sensentr_gmean = masktonull(sensentr_gmean)
    remove_files = [difftemp_file, sensentr_mean_file]
    for filen in remove_files:
        os.remove(filen)
    return sensentr_gmean, sensentr_file


def snowentr(model, wdir, infile, aux_file):
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
    """
    cdo = Cdo()
    prsnmask_file = infile[0]
    removeif(aux_file)
    latsnow_file = wdir + '/{}_latentEnergy_snow.nc'.format(model)
    snowentr_file = wdir + '/{}_snow_entr.nc'.format(model)
    snowentr_mean_file = wdir + '/{}_snowEntropy_gmean.nc'.format(model)
    cdo.mulc(str(LC_SUB),
             input='-setmisstoc,0 {}'.format(prsnmask_file),
             options='-b F32',
             output=latsnow_file)
    flist = [latsnow_file, infile[1], aux_file]
    snowentr_gmean = entr(flist, 'prsn', 'ssnow', snowentr_file,
                          snowentr_mean_file)
    snowentr_gmean = masktonull(snowentr_gmean)
    os.remove(snowentr_mean_file)
    return snowentr_gmean, latsnow_file, snowentr_file


def wmbudg(model, wdir, aux_file, input_data, auxlist):
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
    - auxlist: a list of auxiliary files;
    """
    cdo = Cdo()
    hfls_file = e.select_metadata(input_data, short_name='hfls',
                                  dataset=model)[0]['filename']
    pr_file = e.select_metadata(input_data, short_name='pr',
                                dataset=model)[0]['filename']
    prsn_file = e.select_metadata(input_data, short_name='prsn',
                                  dataset=model)[0]['filename']
    wmbudg_file = wdir + '/{}_wmb.nc'.format(model)
    wm_gmean_file = wdir + '/{}_wmb_gmean.nc'.format(model)
    latene_file = wdir + '/{}_latent.nc'.format(model)
    latene_gmean_file = wdir + '/{}_latent_gmean.nc'.format(model)
    removeif(aux_file)
    cdo.sub(input="{} {}".format(auxlist[0], pr_file), output=aux_file)
    wmass_gmean = write_eb('hfls', 'wmb', aux_file, wmbudg_file, wm_gmean_file)
    removeif(aux_file)
    cdo.sub(input="{} -add -mulc,{} {} -mulc,{} {}".format(
        hfls_file, str(LC_SUB), prsn_file, str(L_C), auxlist[1]),
            output=aux_file)
    latent_gmean = write_eb('hfls', 'latent', aux_file, latene_file,
                            latene_gmean_file)
    varlist = [wmass_gmean, latent_gmean]
    fileout = [wmbudg_file, latene_file]
    remove_files = [wm_gmean_file, latene_gmean_file]
    for filen in remove_files:
        os.remove(filen)
    return varlist, fileout


def write_eb(namein, nameout, aux_file, d3_file, gmean_file):
    """Change variable name in the NetCDF file and compute averages.

    Arguments:
    - namein: initial name of the variable;
    - nameout: final name of the variable;
    - aux_file: the name of an auxiliary file;
    - d3_file: the file containing (time,lat,lon) fields;
    - gmean_file: the name of a file where to put the annual and globally
      averaged fields;
    """
    cdo = Cdo()
    ch_name = '{},{}'.format(namein, nameout)
    cdo.chname(ch_name, input=aux_file, options='-b F32', output=d3_file)
    cdo.fldmean(input='-yearmonmean {}'.format(d3_file), output=gmean_file)
    with Dataset(gmean_file) as f_l:
        constant = f_l.variables[nameout][:]
    return constant
