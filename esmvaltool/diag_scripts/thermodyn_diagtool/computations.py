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
    ----------
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - aux_file: the name of a dummy aux. file to be used for computations;
    - toab_file: a file containing the annual mean TOA energy budgets
      (time,lon,lat);
    - te_file: a file containing the annual mean emission temperature
      (time,lon,lat);

    Returns
    -------
    The annual mean baroclinic efficiency (after Lucarini et al. 2011).

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    removeif(aux_file)
    gain_file = wdir + f"/{model}_maskGain.nc"
    cdo.gtc("0", input=toab_file, output=gain_file)
    loss_file = wdir + f"/{model}_maskLoss.nc"
    cdo.ltc("0", input=toab_file, output=loss_file)
    toabgain_file = wdir + f"/{model}_toabGain.nc"
    cdo.setrtomiss(
        "-1000,0", input=f"-mul {toab_file} {gain_file}", output=toabgain_file
    )
    toabloss_file = wdir + f"/{model}_toabLoss.nc"
    cdo.setrtomiss(
        "0,1000", input=f"-mul {toab_file} {loss_file}", output=toabloss_file
    )
    tegain_file = wdir + f"/{model}_teGain.nc"
    cdo.setrtomiss(
        "-1000,0", input=f"-mul {te_file} {gain_file}", output=tegain_file
    )
    teloss_file = wdir + f"/{model}_teLoss.nc"
    cdo.setrtomiss(
        "-1000,0", input=f"-mul {te_file} {loss_file}", output=teloss_file
    )
    tegainm_file = wdir + f"/{model}_teGainm.nc"
    cdo.div(
        input=f"-fldmean {toabgain_file} -fldmean -div {toabgain_file} {tegain_file} ",
        output=tegainm_file,
    )
    telossm_file = wdir + f"/{model}_teLossm.nc"
    cdo.div(
        input=f"-fldmean {toabloss_file} -fldmean -div {toabloss_file} {teloss_file} ",
        output=telossm_file,
    )
    aux_baroceff_file = wdir + f"/{model}_aux_barocEff.nc"
    cdo.sub(
        input=f"-reci {telossm_file} -reci {tegainm_file}",
        output=aux_baroceff_file,
    )
    baroceff_file = wdir + f"/{model}_barocEff.nc"
    cdo.div(
        input=f"{aux_baroceff_file} -mulc,0.5 -add -reci {tegainm_file} -reci {telossm_file}",
        output=baroceff_file,
    )
    with Dataset(baroceff_file) as f_l:
        baroc = f_l.variables["toab"][0, 0, 0]
    remove_files = [
        gain_file,
        loss_file,
        toabgain_file,
        toabloss_file,
        tegain_file,
        teloss_file,
        tegainm_file,
        telossm_file,
        aux_baroceff_file,
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
    ----------
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - aux_file: the name of a dummy aux. file to be used for computations;
    - filelist: a list of file names containing the input fields;

    Returns
    -------
    The list of input files, the global mean budget time series, a file
    containing the budget fields, a file containing the annual mean TOA budget
    value;

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    hfls_file = e.select_metadata(
        input_data, short_name="hfls", dataset=model
    )[0]["filename"]
    hfss_file = e.select_metadata(
        input_data, short_name="hfss", dataset=model
    )[0]["filename"]
    rlds_file = e.select_metadata(
        input_data, short_name="rlds", dataset=model
    )[0]["filename"]
    rlus_file = e.select_metadata(
        input_data, short_name="rlus", dataset=model
    )[0]["filename"]
    rlut_file = e.select_metadata(
        input_data, short_name="rlut", dataset=model
    )[0]["filename"]
    rsds_file = e.select_metadata(
        input_data, short_name="rsds", dataset=model
    )[0]["filename"]
    rsus_file = e.select_metadata(
        input_data, short_name="rsus", dataset=model
    )[0]["filename"]
    rsdt_file = e.select_metadata(
        input_data, short_name="rsdt", dataset=model
    )[0]["filename"]
    rsut_file = e.select_metadata(
        input_data, short_name="rsut", dataset=model
    )[0]["filename"]
    toab_file = wdir + f"/{model}_toab.nc"
    toab_gmean_file = wdir + f"/{model}_toab_gmean.nc"
    surb_file = wdir + f"/{model}_surb.nc"
    aux_surb_file = wdir + f"/{model}_aux_surb.nc"
    surb_gmean_file = wdir + f"/{model}_surb_gmean.nc"
    atmb_file = wdir + f"/{model}_atmb.nc"
    atmb_gmean_file = wdir + f"/{model}_atmb_gmean.nc"
    removeif(aux_file)
    cdo.sub(input=f"-sub {rsdt_file} {rsut_file} {rlut_file}", output=aux_file)
    toab_gmean = write_eb("rsdt", "toab", aux_file, toab_file, toab_gmean_file)
    toab_ymm_file = wdir + f"/{model}_toab_ymm.nc"
    cdo.yearmonmean(input=toab_file, output=toab_ymm_file)
    # Surface energy budget
    removeif(aux_file)
    cdo.add(input=f" {rsds_file} {rlds_file}", output=aux_surb_file)
    cdo.sub(
        input=f"-sub -sub -sub {aux_surb_file} {rsus_file} {rlus_file} {hfls_file} {hfss_file}",
        output=aux_file,
    )
    surb_gmean = write_eb("rsds", "surb", aux_file, surb_file, surb_gmean_file)
    # Atmospheric energy budget
    removeif(aux_file)
    cdo.sub(input=f"{toab_file} {surb_file}", output=aux_file)
    atmb_gmean = write_eb("toab", "atmb", aux_file, atmb_file, atmb_gmean_file)
    eb_gmean = [toab_gmean, atmb_gmean, surb_gmean]
    eb_file = [toab_file, atmb_file, surb_file]
    # Delete files
    filenames = [
        aux_surb_file,
        toab_gmean_file,
        atmb_gmean_file,
        surb_gmean_file,
    ]
    for filen in filenames:
        os.remove(filen)
    input_list = [
        hfls_file,
        hfss_file,
        rlds_file,
        rlus_file,
        rlut_file,
        rsds_file,
        rsdt_file,
        rsus_file,
        rsut_file,
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
    ----------
    logger: the log file where the global mean values are printed out;
    model: the model name;
    wdir: the working directory where the outputs are stored;
    filelist: the list containing all the input files;
    aux_file: the name of a dummy aux. file to be used for computations;
    lect: the annual mean value of the LEC strength;
    flags: a list of flags containing information on whether the water mass
           and energy budgets are computed, if the material entropy production
           has to be computed, if using the indirect, the direct method, or
           both methods;

    Returns
    -------
    The annual mean entropy production with the direct method, the degree of
    irreversibility, the list of input files for the computation.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    lec = flags[1]
    aux_files = mkthe.init_mkthe_direntr(
        model, wdir, input_data, te_file, flags
    )
    htop_file = aux_files[1]
    prr_file = aux_files[2]
    tabl_file = aux_files[3]
    tasvert_file = aux_files[4]
    tcloud_file = aux_files[5]
    tcolumn_file = aux_files[6]
    tlcl_file = aux_files[7]
    hfls_file = e.select_metadata(
        input_data, short_name="hfls", dataset=model
    )[0]["filename"]
    hfss_file = e.select_metadata(
        input_data, short_name="hfss", dataset=model
    )[0]["filename"]
    prsn_file = e.select_metadata(
        input_data, short_name="prsn", dataset=model
    )[0]["filename"]
    ts_file = e.select_metadata(input_data, short_name="ts", dataset=model)[0][
        "filename"
    ]
    logger.info(
        "Computation of the material entropy "
        "production with the direct method\n"
    )
    logger.info("1. Sensible heat fluxes\n")
    infile_list = [hfss_file, tabl_file, ts_file]
    ssens, sensentr_file = sensentr(model, wdir, infile_list, aux_file)
    logger.info(
        "Material entropy production associated with sens. heat fluxes: %s\n",
        ssens,
    )
    logger.info("2. Hydrological cycle\n")
    logger.info("2.1 Evaporation fluxes\n")
    infile_list = [hfls_file, ts_file]
    sevap, evapentr_file = evapentr(model, wdir, infile_list, aux_file)
    logger.info(
        "Material entropy production associated with evaporation fluxes: %s\n",
        sevap,
    )
    infile_mask = [prr_file, prsn_file, tlcl_file]
    prrmask_file, prsnmask_file = mask_precip(model, wdir, infile_mask)
    logger.info("2.2 Rainfall precipitation\n")
    infile_rain = [prrmask_file, tcloud_file]
    srain, rainentr_file = rainentr(model, wdir, infile_rain, aux_file)
    logger.info(
        "Material entropy production associated with rainfall: %s\n", srain
    )
    logger.info("2.3 Snowfall precipitation\n")
    infile_snow = [prsnmask_file, tcloud_file]
    ssnow, latsnow_file, snowentr_file = snowentr(
        model, wdir, infile_snow, aux_file
    )
    logger.info(
        "Material entropy production associated with snowfall: %s\n", ssnow
    )
    logger.info("2.4 Melting of snow at the surface \n")
    smelt, meltentr_file = meltentr(model, wdir, latsnow_file, aux_file)
    logger.info(
        "Material entropy production associated with snow melting: %s\n", smelt
    )
    logger.info("2.5 Potential energy of the droplet\n")
    infile_pot = [htop_file, prrmask_file, prsnmask_file, tcolumn_file]
    spot, potentr_file = potentr(model, wdir, infile_pot, aux_file)
    logger.info(
        "Material entropy production associated with "
        "potential energy of the droplet: %s\n",
        spot,
    )
    os.remove(prrmask_file)
    os.remove(prsnmask_file)
    logger.info("3. Kinetic energy dissipation\n")
    skin = kinentr(logger, aux_file, tasvert_file, lect, lec)
    matentr = (
        float(ssens)
        - float(sevap)
        + float(srain)
        + float(ssnow)
        + float(spot)
        + float(skin)
        - float(smelt)
    )
    logger.info(
        "Material entropy production with the direct method: %s\n", matentr
    )
    irrevers = (matentr - float(skin)) / float(skin)
    for filen in aux_files:
        os.remove(filen)
    entr_list = [
        sensentr_file,
        evapentr_file,
        rainentr_file,
        snowentr_file,
        meltentr_file,
        potentr_file,
    ]
    return matentr, irrevers, entr_list


def entr(filelist, nin, nout, entr_file, entr_mean_file):
    """Obtain the entropy dividing some energy by some working temperature.

    This function ingests an energy and a related temperature, then writes
    (time,lat,lon) entropy fluxes and entropy flux annual mean values to NC
    files.

    Arguments:
    ----------
    filelist: a list of file containing the name of the energy file, of the
              temperature file and of an auxiliary file needed for computation;
    nin: the variable name of the input energy fields;
    nout: the variable name to attribute to the entropy flux in the NC file;
    entr_file: the name of the file containing the 3D entropy fluxes;
    entr_mean_file: the name of the file containing the global annual mean
                    entropy value;

    Returns
    -------
    The annual global mean value of entropy.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    en_file = filelist[0]
    tem_file = filelist[1]
    aux_file = filelist[2]
    removeif(aux_file)
    cdo.timmean(
        input=f"-yearmonmean -monmean -div {en_file} {tem_file}",
        options="-b F32",
        output=aux_file,
    )
    entr_gmean = write_eb(nin, nout, aux_file, entr_file, entr_mean_file)
    return entr_gmean


def evapentr(model, wdir, infile, aux_file):
    """Compute entropy production related to evaporation fluxes.

    The function computes the material entropy production related to
    evaporation fluxes, as part of the material entropy production
    obtained with the direct method (after Lucarini et al., 2011).

    Arguments:
    ----------
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: a list of file containing hfls and ts, respectively
      (with dimensions (time,lat,lon);
    - aux_file: the name of a dummy aux. file to be used for computations;

    Returns
    -------
    The global annual mean entropy production related to evaporation, the
    file containing it.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    evapentr_file = wdir + f"/{model}_evap_entr.nc"
    evapentr_mean_file = wdir + f"/{model}_evapEntropy_gmean.nc"
    flist = [infile[0], infile[1], aux_file]
    evapentr_gmean = entr(
        flist, "hfls", "sevap", evapentr_file, evapentr_mean_file
    )
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: a list of files, containing each the fields rlds, rlus, rsds,
            rsus, emission temperature (te), TOA energy budget (toab) and ts;
    toab_file: a file containing the annual mean TOA energy budgets
              (time,lon,lat);
    aux_file: the name of a dummy aux. file to be used for computations;
    toab_gmean: the climatological annaul mean TOA energy budget;

    Returns
    -------
    The annual mean vertical and horizontal components of the entropy
    production with the indirect method, the file containing them.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    rlds_file = e.select_metadata(
        input_data, short_name="rlds", dataset=model
    )[0]["filename"]
    rlus_file = e.select_metadata(
        input_data, short_name="rlus", dataset=model
    )[0]["filename"]
    rsds_file = e.select_metadata(
        input_data, short_name="rsds", dataset=model
    )[0]["filename"]
    rsus_file = e.select_metadata(
        input_data, short_name="rsus", dataset=model
    )[0]["filename"]
    ts_file = e.select_metadata(input_data, short_name="ts", dataset=model)[0][
        "filename"
    ]
    horzentropy_file = wdir + f"/{model}_horizEntropy.nc"
    vertenergy_file = wdir + f"/{model}_verticalEnergy.nc"
    vertentropy_file = wdir + f"/{model}_verticalEntropy.nc"
    vertentropy_mean_file = wdir + f"/{model}_vertEntropy_gmean.nc"
    horzentropy_mean_file = wdir + f"/{model}_horizEntropy_gmean.nc"
    removeif(aux_file)
    cdo.yearmonmean(
        input=f"-mulc,-1 -div -subc,{np.nanmean(toab_gmean)}  {infile[1]}  {infile[0]}",
        output=aux_file,
    )
    horzentr_mean = write_eb(
        "toab", "shor", aux_file, horzentropy_file, horzentropy_mean_file
    )
    cdo.yearmonmean(
        input=f" -add {rlds_file} -sub {rsds_file} -add {rlus_file} {rsus_file}",
        output=vertenergy_file,
    )
    cdo.mul(
        input=f"{vertenergy_file} -sub -yearmonmean -reci {infile[0]} -yearmonmean -reci {ts_file}",
        output=aux_file,
    )
    vertentr_mean = write_eb(
        "rlds", "sver", aux_file, vertentropy_file, vertentropy_mean_file
    )
    remove_files = [
        horzentropy_mean_file,
        vertenergy_file,
        vertentropy_mean_file,
    ]
    for filen in remove_files:
        os.remove(filen)
    return horzentr_mean, vertentr_mean, horzentropy_file, vertentropy_file


def kinentr(logger, aux_file, tasvert_file, lect, lec):
    """Compute the material entropy production from kin. energy dissipation.

    The function computes the material entropy production associated with the
    kinetic energy dissipation, through the intensity of the LEC.

    Arguments:
    ----------
    aux_file: the name of a dummy aux. file to be used for computations;
    tasvert_file: a file containing the vertically integrated boundary layer
                  temperature;
    lect: an array containing the annual mean LEC intensity;
    lec: a flag marking whether the LEC has been previously computed or not

    Returns
    -------
    The global annual mean entropy production related to kinetic energy
    dissipation.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    removeif(aux_file)
    if lec == "True":
        cdo.yearmonmean(input=tasvert_file, output=aux_file)
        with Dataset(aux_file) as f_l:
            tabl_mean = f_l.variables["ts"][:, 0, 0]
        minentr_mean = np.nanmean(lect / tabl_mean)
        logger.info(
            "Material entropy production associated with "
            "kinetic energy dissipation: %s\n",
            minentr_mean,
        )
        minentr_mean = masktonull(minentr_mean)
    else:
        minentr_mean = 0.010
        logger.info(
            "I cannot compute the material entropy "
            "production without the LEC...\n"
        )
        logger.info(
            "I will assign a given value for the material "
            "entropy production attributed to LEC "
            "(0.01 W/m2*K)\n"
        )
    return minentr_mean


def landoc_budg(model, wdir, infile, mask, name):
    """Compute budgets separately on land and oceans.

    Arguments:
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: the file containing the original budget field as (time,lat,lon);
    mask: the file containing the land-sea mask;
    name: the variable name as in the input file;

    Returns
    -------
    The mean budgets over land and over oceans.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    ocean_file = wdir + f"/{model}_{name}_ocean.nc"
    oc_gmean_file = wdir + f"/{model}_{name}_oc_gmean.nc"
    land_file = wdir + f"/{model}_{name}_land.nc"
    la_gmean_file = wdir + f"/{model}_{name}_la_gmean.nc"
    aux_file = wdir + "/aux.nc"
    removeif(aux_file)
    cdo.mul(input=f"{infile} -eqc,0 {mask}", output=ocean_file)
    cdo.timmean(input=f"-fldmean {ocean_file}", output=oc_gmean_file)
    with Dataset(oc_gmean_file) as f_l:
        oc_gmean = f_l.variables[name][0, 0, 0]
    cdo.sub(input=f"{infile} {ocean_file}", output=land_file)
    cdo.setctomiss("0", input=ocean_file, output=aux_file)
    move(aux_file, ocean_file)
    cdo.setctomiss("0", input=land_file, output=aux_file)
    move(aux_file, land_file)
    cdo.timmean(input=f"-fldmean {land_file}", output=la_gmean_file)
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: a list of input file, containing rainfall precipitation (prr) and
            prsn, respectively (dimensions (time,lat,lon));

    Returns
    -------
    The files containing masked rainfall and snowfall precipitation fields,
    respectively.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    prr_file = infile[0]
    prsn_file = infile[1]
    tlcl_file = infile[2]
    # Prepare masks for snowfall and rainfall
    maskrain_file = wdir + f"/{model}_maskprecr.nc"
    cdo.gtc("1.0E-7", input=prr_file, options=" -b F32", output=maskrain_file)
    masksnow_file = wdir + f"/{model}_maskprecs.nc"
    cdo.gtc("1.0E-7", input=prsn_file, options=" -b F32", output=masksnow_file)
    prrmask_file = wdir + f"/{model}_prr_masked.nc"
    cdo.mul(
        input=f"{maskrain_file} {prr_file}",
        options="-b F32",
        output=prrmask_file,
    )
    prsnmask_file = wdir + f"/{model}_prsn_masked.nc"
    cdo.mul(
        input=f"{masksnow_file} {prsn_file}",
        options="-b F32",
        output=prsnmask_file,
    )
    # Temperatures of the rainfall and snowfall clouds
    tliq_file = wdir + f"/{model}_tliq.nc"
    cdo.setrtomiss(
        "-1000,0",
        input=f"-mul {tlcl_file} {maskrain_file}",
        options="-b F32",
        output=tliq_file,
    )
    tsol_file = wdir + f"/{model}_tsol.nc"
    cdo.setrtomiss(
        "-1000,0",
        input=f"-mul {tlcl_file} {masksnow_file}",
        options="-b F32",
        output=tsol_file,
    )
    tdegl_file = wdir + f"/{model}_tliqdeg.nc"
    cdo.subc("273.15", input=tliq_file, options="-b F32", output=tdegl_file)
    tdegs_file = wdir + f"/{model}_tsoldeg.nc"
    cdo.subc("273.15", input=tsol_file, options="-b F32", output=tdegs_file)
    # Mask for ice cloud and temperature for phase changes from ice to rain
    maskice_file = wdir + f"/{model}_maskice.nc"
    cdo.ltc("0.0", input=tdegl_file, options="-b F32", output=maskice_file)
    ticer_file = wdir + f"/{model}_t_icerain_file"
    cdo.setrtomiss(
        "-1000,0",
        input=f"-mul {tliq_file} {maskice_file}",
        options="-b F32",
        output=ticer_file,
    )
    prrice_file = wdir + f"/{model}_prr_ice_file.nc"
    cdo.mul(
        input=f"{maskice_file} {prr_file}",
        options="-b F32",
        output=prrice_file,
    )
    # Mask for vapor cloud and temperature for phase changes from vapor to snow
    maskvap_file = wdir + f"/{model}_maskvap.nc"
    cdo.gtc("0.0", input=tdegs_file, options="-b F32", output=maskvap_file)
    tvaps_file = wdir + f"/{model}_t_vapsnow.nc"
    cdo.setrtomiss(
        "-1000,0",
        input=f"-mul {tsol_file} {maskvap_file}",
        options="-b F32",
        output=tvaps_file,
    )
    prsnvap_file = wdir + f"/{model}_prsn_vap.nc"
    cdo.mul(
        input=f"{maskvap_file} {prsn_file}",
        options="-b F32",
        output=prsnvap_file,
    )
    remove_files = [
        maskrain_file,
        masksnow_file,
        tliq_file,
        tsol_file,
        tdegl_file,
        tdegs_file,
        maskice_file,
        ticer_file,
        prrice_file,
        maskvap_file,
        tvaps_file,
        prsnvap_file,
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
    ----------
    - model: the model name;
    - wdir: the working directory where the outputs are stored;
    - infile: the latent energy associated with snowfall precipitation;
    - aux_file: the name of a dummy aux. file to be used for computations;

    Returns
    -------
    The global annual mean entropy production related to evaporation, the
    file containing it.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    removeif(aux_file)
    latmelt_file = wdir + f"/{model}_latentEnergy_snowmelt.nc"
    meltentr_file = wdir + f"/{model}_snowmelt_entr.nc"
    meltentr_mean_file = wdir + f"/{model}_snowmeltEntropy_gmean.nc"
    cdo.mulc(
        str(L_S),
        input=f"-divc,{str(LC_SUB)} {latsnow_file}",
        options="-b F32",
        output=latmelt_file,
    )
    cdo.timmean(
        input=f"-yearmonmean -monmean -setmisstoc,0 -divc,273.15 {latmelt_file}",
        options="-b F32",
        output=aux_file,
    )
    cdo.chname(
        "prsn,smelt", input=aux_file, options="-b F32", output=meltentr_file
    )
    cdo.fldmean(
        input=meltentr_file, options="-b F32", output=meltentr_mean_file
    )
    with Dataset(meltentr_mean_file) as f_l:
        meltentr_gmean = f_l.variables["smelt"][0, 0, 0]
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: a list of files containing the height of the bondary layer top
            (htop), the masked rainfall precipitation (prrmask), the masked
            snowfall precipitation (prsnmask), the temperature of the vertical
            column between the cloud top and the ground (tcolumn);
    aux_file: the name of a dummy aux. file to be used for computations;

    Returns
    -------
    The global annual mean entropy production related to potential energy of
    the droplet, the file containing it.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    removeif(aux_file)
    htop_file = infile[0]
    prrmask_file = infile[1]
    prsnmask_file = infile[2]
    tcolumn_file = infile[3]
    poten_file = wdir + f"/{model}_potEnergy_drop.nc"
    potentr_file = wdir + f"/{model}_pot_drop_entr.nc"
    potentr_mean_file = wdir + f"/{model}_potEnergy_drop_gmean.nc"
    cdo.mulc(
        GRAV,
        input=f"-mul {htop_file} -add {prrmask_file} {prsnmask_file}",
        options="-b F32",
        output=poten_file,
    )
    flist = [poten_file, tcolumn_file, aux_file]
    potentr_gmean = entr(
        flist, "htop", "spotp", potentr_file, potentr_mean_file
    )
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: a list of file containing the masked rainfall precipitation
            (prrmask) and the temperature of the cloud (tcloud);
    aux_file: the name of a dummy aux. file to be used for computations;

    Returns
    -------
    The global annual mean entropy production related to rainfall, the
    file containing it.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    prrmask_file = infile[0]
    removeif(aux_file)
    latrain_file = wdir + f"/{model}_latentEnergy_rain.nc"
    rainentr_file = wdir + f"/{model}_rain_entr.nc"
    rainentr_mean_file = wdir + f"/{model}_rainEntropy_gmean.nc"
    cdo.mulc(
        str(L_C),
        input=f"-setmisstoc,0 {prrmask_file}",
        options="-b F32",
        output=latrain_file,
    )
    flist = [latrain_file, infile[1], aux_file]
    rainentr_gmean = entr(
        flist, "prr", "srain", rainentr_file, rainentr_mean_file
    )
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: a list of file containing hfss, the temperature at the boundary
            layer top (tabl), ts, respectively (with dimensions (time,lat,lon);
    aux_file: the name of a dummy aux. file to be used for computations;

    Returns
    -------
    The global annual mean entropy production related to sensible heat fluxes,
    the file containing it.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    difftemp_file = wdir + f"/{model}_difftemp_bl.nc"
    sensentr_file = wdir + f"/{model}_sens_entr.nc"
    sensentr_mean_file = wdir + f"/{model}_sensEntropy_gmean.nc"
    cdo.reci(
        input=f"-sub -reci {infile[1]}  -reci {infile[2]}",
        options="-b F32",
        output=difftemp_file,
    )
    flist = [infile[0], difftemp_file, aux_file]
    sensentr_gmean = entr(
        flist, "hfss", "ssens", sensentr_file, sensentr_mean_file
    )
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    infile: a list of file containing the masked snowfall precipitation
            (prsnmask) and the temperature of the cloud (tcloud);
    aux_file: the name of a dummy aux. file to be used for computations;

    Returns
    -------
    The global annual mean entropy production related to snowfall, the
    file containing it.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    prsnmask_file = infile[0]
    removeif(aux_file)
    latsnow_file = wdir + f"/{model}_latentEnergy_snow.nc"
    snowentr_file = wdir + f"/{model}_snow_entr.nc"
    snowentr_mean_file = wdir + f"/{model}_snowEntropy_gmean.nc"
    cdo.mulc(
        str(LC_SUB),
        input=f"-setmisstoc,0 {prsnmask_file}",
        options="-b F32",
        output=latsnow_file,
    )
    flist = [latsnow_file, infile[1], aux_file]
    snowentr_gmean = entr(
        flist, "prsn", "ssnow", snowentr_file, snowentr_mean_file
    )
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
    ----------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    aux_file: the name of a dummy aux. file to be used for computations;
    input_data: a dictionary of file names containing the input fields;
    auxlist: a list of auxiliary files;

    Returns
    -------
    A list containing global mean water mass and latent energy values, a list
    of files containing them.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    hfls_file = e.select_metadata(
        input_data, short_name="hfls", dataset=model
    )[0]["filename"]
    pr_file = e.select_metadata(input_data, short_name="pr", dataset=model)[0][
        "filename"
    ]
    prsn_file = e.select_metadata(
        input_data, short_name="prsn", dataset=model
    )[0]["filename"]
    wmbudg_file = wdir + f"/{model}_wmb.nc"
    wm_gmean_file = wdir + f"/{model}_wmb_gmean.nc"
    latene_file = wdir + f"/{model}_latent.nc"
    latene_gmean_file = wdir + f"/{model}_latent_gmean.nc"
    removeif(aux_file)
    cdo.sub(input=f"{auxlist[0]} {pr_file}", output=aux_file)
    wmass_gmean = write_eb("hfls", "wmb", aux_file, wmbudg_file, wm_gmean_file)
    removeif(aux_file)
    cdo.sub(
        input=f"{hfls_file} -add -mulc,{str(LC_SUB)} {prsn_file} -mulc,{str(L_C)} {auxlist[1]}",
        output=aux_file,
    )
    latent_gmean = write_eb(
        "hfls", "latent", aux_file, latene_file, latene_gmean_file
    )
    varlist = [wmass_gmean, latent_gmean]
    fileout = [wmbudg_file, latene_file]
    remove_files = [wm_gmean_file, latene_gmean_file]
    for filen in remove_files:
        os.remove(filen)
    return varlist, fileout


def write_eb(namein, nameout, aux_file, d3_file, gmean_file):
    """Change variable name in the NetCDF file and compute averages.

    Arguments:
    ----------
    namein: initial name of the variable;
    nameout: final name of the variable;
    aux_file: the name of an auxiliary file;
    d3_file: the file containing (time,lat,lon) fields;
    gmean_file: the name of a file where to put the annual and globally
                averaged fields;

    Returns
    -------
    A global annual mean value of the budget.

    @author: Valerio Lembo, Hamburg University, 2018.
    """
    cdo = Cdo()
    ch_name = f"{namein},{nameout}"
    cdo.chname(ch_name, input=aux_file, options="-b F32", output=d3_file)
    cdo.fldmean(input=f"-yearmonmean {d3_file}", output=gmean_file)
    with Dataset(gmean_file) as f_l:
        constant = f_l.variables[nameout][:]
    return constant
