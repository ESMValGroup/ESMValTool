"""AUXILIARY FIELDS RETRIEVAL.

Module for computation of the auxiliary variables needed by the tool.

It contains the following functions:
- init_mkthe_te: compute emission temperature from OLR;
- init_mkthe_wat: initialise wfluxes;
- init_mkthe_lec: compute monthly mean near-surface zonal and meridional
                  velocities when daily means are provided;
- init_mkthe_direntr: compute auxiliary files needed for material entropy
                      production retrieval with the direct method;
- input_fields: obtain input fields for mkthe_main;
- mkthe_main: obtain equivalent potential temperatures, temperatures
              representative of the sensible and latent heat exchanges in
              the lower layers of the troposphere, boundary layer height and
              lifting condensation level temperature.
- mon_from_day: obtain monthly means from daily means;
- wfluxes: obtain evaporation and precipitation from precipitation and latent
           heat fluxes;
- write_output: write auxiliary fields to NetCDF file;

@author: Valerio Lembo, valerio.lembo@uni-hamburg.de, Universitat Hamburg, 2018
"""

import os
from shutil import move

import numpy as np
from cdo import Cdo
from netCDF4 import Dataset

import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.thermodyn_diagtool import fourier_coefficients

ALV = 2.5008e6  # Latent heat of vaporization
G_0 = 9.81  # Gravity acceleration
P_0 = 100000.0  # reference pressure
RV = 461.51  # Gas constant for water vapour
T_MELT = 273.15  # freezing temp.
AKAP = 0.286  # Kappa (Poisson constant R/Cp)
GAS_CON = 287.0  # Gas constant
RA_1 = 610.78  # Parameter for Magnus-Teten-Formula
H_S = 300.0  # stable boundary layer height (m)
H_U = 1000.0  # unstable boundary layer height (m)
RIC_RS = 0.39  # Critical Richardson number for stable layer
RIC_RU = 0.28  # Critical Richardson number for unstable layer
L_C = 2501000  # latent heat of condensation
SIGMAINV = 17636684.3034  # inverse of the Stefan-Boltzmann constant


def init_mkthe_te(model, wdir, input_data):
    """Compute auxiliary fields or perform time averaging of existing fields.

    Arguments:
    ---------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    filelist: a list of file names containing the input fields;

    Returns
    -------
    A file containing annual mean emission temperature fields, the time mean
    globally averaged emission temperature, the file containing emission
    temperature fields.
    """
    cdo = Cdo()
    rlut_file = e.select_metadata(
        input_data, short_name="rlut", dataset=model
    )[0]["filename"]
    # Compute monthly mean fields from 2D surface daily fields
    # emission temperature
    te_file = wdir + f"/{model}_te.nc"
    cdo.sqrt(input=f"-sqrt -mulc,{SIGMAINV} {rlut_file}", output=te_file)
    te_ymm_file = wdir + f"/{model}_te_ymm.nc"
    cdo.yearmonmean(input=te_file, output=te_ymm_file)
    te_gmean_file = wdir + f"/{model}_te_gmean.nc"
    cdo.timmean(input=f"-fldmean {te_ymm_file}", output=te_gmean_file)
    with Dataset(te_gmean_file) as f_l:
        te_gmean_constant = f_l.variables["rlut"][0, 0, 0]
    return te_ymm_file, te_gmean_constant, te_file


def init_mkthe_wat(model, wdir, input_data, flags):
    """Compute auxiliary fields or perform time averaging of existing fields.

    Arguments:
    ---------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    filelist: a list of file names containing the input fields;
    flags: (wat: a flag for the water mass budget module (y or n),
            entr: a flag for the material entropy production (y or n);
            met: a flag for the material entropy production method
            (1: indirect, 2, direct, 3: both));

    Returns
    -------
    A list of input fields.
    """
    wat = flags[0]
    if wat == "True":
        evspsbl_file, prr_file = wfluxes(model, wdir, input_data)
        aux_files = [evspsbl_file, prr_file]
    return aux_files


def init_mkthe_lec(model, wdir, input_data):
    """Compute auxiliary fields or perform time averaging of existing fields.

    Arguments:
    ---------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    filelist: a dictionary of file names containing the input fields;
    flags: (wat: a flag for the water mass budget module (y or n),
            entr: a flag for the material entropy production (y or n);
            met: a flag for the material entropy production method
            (1: indirect, 2, direct, 3: both));

    Returns
    -------
    The file containing monthly mean near-surface velocities in the zonal and
    meridional direction.
    """
    uas_file = e.select_metadata(input_data, short_name="uas", dataset=model)[
        0
    ]["filename"]
    vas_file = e.select_metadata(input_data, short_name="vas", dataset=model)[
        0
    ]["filename"]
    uasmn_file = mon_from_day(wdir, model, "uas", uas_file)
    vasmn_file = mon_from_day(wdir, model, "vas", vas_file)
    return uasmn_file, vasmn_file


def init_mkthe_direntr(model, wdir, input_data, te_file, flags):
    """Compute the MEP with the direct method.

    Arguments:
    ---------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    input_data: a dictionary of file names containing the input fields;
    te_file: a file containing the emission temperature computed from OLR;
    flags: (wat: a flag for the water mass budget module (y or n),
            entr: a flag for the material entropy production (y or n);
            met: a flag for the material entropy production method
            (1: indirect, 2, direct, 3: both));

    Returns
    -------
    A list of files cotiaining the components of the MEP with the direct
    method.
    """
    cdo = Cdo()
    met = flags[3]
    if met in {"2", "3"}:
        evspsbl_file, prr_file = wfluxes(model, wdir, input_data)
        hfss_file = e.select_metadata(
            input_data, short_name="hfss", dataset=model
        )[0]["filename"]
        hus_file = e.select_metadata(
            input_data, short_name="hus", dataset=model
        )[0]["filename"]
        ps_file = e.select_metadata(
            input_data, short_name="ps", dataset=model
        )[0]["filename"]
        ts_file = e.select_metadata(
            input_data, short_name="ts", dataset=model
        )[0]["filename"]
        uas_file = e.select_metadata(
            input_data, short_name="uas", dataset=model
        )[0]["filename"]
        uas_tres = e.select_metadata(
            input_data, short_name="uas", dataset=model
        )[0]["mip"]
        vas_file = e.select_metadata(
            input_data, short_name="vas", dataset=model
        )[0]["filename"]
        vas_tres = e.select_metadata(
            input_data, short_name="vas", dataset=model
        )[0]["mip"]
        if uas_tres == "day":
            uasmn_file = wdir + f"/{model}_uas_mm.nc"
            uasmn_file = mon_from_day(wdir, model, "uas", uas_file)
            uas_file = uasmn_file
        if vas_tres == "day":
            vasmn_file = wdir + f"/{model}_uas_mm.nc"
            vasmn_file = mon_from_day(wdir, model, "vas", vas_file)
            vas_file = vasmn_file
        mk_list = [
            ts_file,
            hus_file,
            ps_file,
            uas_file,
            vas_file,
            hfss_file,
            te_file,
        ]
        htop_file, tabl_file, tlcl_file = mkthe_main(wdir, mk_list, model)
        # Working temperatures for the hydrological cycle
        tcloud_file = wdir + f"/{model}_tcloud.nc"
        removeif(tcloud_file)
        cdo.mulc(
            "0.5",
            input=f"-add {tlcl_file} {te_file}",
            options="-b F32",
            output=tcloud_file,
        )
        tcolumn_file = wdir + f"/{model}_t_vertav_pot.nc"
        removeif(tcolumn_file)
        cdo.mulc(
            "0.5",
            input=f"-add {ts_file} {tcloud_file}",
            options="-b F32",
            output=tcolumn_file,
        )
        # Working temperatures for the kin. en. diss. (updated)
        tasvert_file = wdir + f"/{model}_tboundlay.nc"
        removeif(tasvert_file)
        cdo.fldmean(
            input=f"-mulc,0.5 -add {ts_file} {tabl_file}",
            options="-b F32",
            output=tasvert_file,
        )
        aux_files = [
            evspsbl_file,
            htop_file,
            prr_file,
            tabl_file,
            tasvert_file,
            tcloud_file,
            tcolumn_file,
            tlcl_file,
        ]
    else:
        aux_files = []
    return aux_files


def input_fields(wdir, file_list):
    """Manipulate input fields and read datasets.

    Arguments:
    ---------
    wdir: the working directory path;
    file_list: the list of file containing ts, hus,
               ps, uas, vas, hfss, te;

    Returns
    -------
    hfss, huss, ps, te, ts and teh near-surface wind speed fields.
    """
    cdo = Cdo()
    ts_miss_file = wdir + "/ts.nc"
    removeif(ts_miss_file)
    cdo.setctomiss("0", input=file_list[0], output=ts_miss_file)
    hus_miss_file = wdir + "/hus.nc"
    removeif(hus_miss_file)
    cdo.setctomiss("0", input=file_list[1], output=hus_miss_file)
    ps_miss_file = wdir + "/ps.nc"
    removeif(ps_miss_file)
    cdo.setctomiss("0", input=file_list[2], output=ps_miss_file)
    vv_missfile = wdir + "/V.nc"
    removeif(vv_missfile)
    vv_file = wdir + "/V_miss.nc"
    removeif(vv_file)
    cdo.sqrt(
        input=f"-add -sqr {file_list[3]} -sqr {file_list[4]}",
        options="-b F32",
        output=vv_file,
    )
    cdo.setctomiss("0", input=vv_file, output=vv_missfile)
    os.remove(vv_file)
    hfss_miss_file = wdir + "/hfss.nc"
    removeif(hfss_miss_file)
    cdo.setctomiss("0", input=file_list[5], output=hfss_miss_file)
    te_miss_file = wdir + "/te.nc"
    removeif(te_miss_file)
    cdo.setctomiss("0", input=file_list[6], output=te_miss_file)
    with Dataset(ts_miss_file) as dataset:
        t_s = dataset.variables["ts"][:, :, :]
    with Dataset(hus_miss_file) as dataset:
        hus = dataset.variables["hus"][:, :, :, :]
        lev = dataset.variables["plev"][:]
    with Dataset(ps_miss_file) as dataset:
        p_s = dataset.variables["ps"][:, :, :]
    with Dataset(vv_missfile) as dataset:
        vv_hor = dataset.variables["uas"][:, :, :]
    with Dataset(hfss_miss_file) as dataset:
        hfss = dataset.variables["hfss"][:, :, :]
    with Dataset(te_miss_file) as dataset:
        t_e = dataset.variables["rlut"][:, :, :]
    huss = hus[:, 0, :, :]
    huss = np.where(lev[0] >= p_s, huss, 0.0)
    nlev = len(lev)
    for l_l in range(nlev):
        aux = hus[:, l_l, :, :]
        aux = np.where((p_s >= lev[l_l]), aux, 0.0)
        huss = huss + aux
    remove_files = [
        ts_miss_file,
        hus_miss_file,
        ps_miss_file,
        vv_missfile,
        hfss_miss_file,
        te_miss_file,
    ]
    for filen in remove_files:
        os.remove(filen)
    return hfss, huss, p_s, t_e, t_s, vv_hor


def mkthe_main(wdir, file_list, modelname):
    """Compute the auxiliary variables for the Thermodynamic diagnostic tool.

    It computes equivalent potential temperatures and temperatures
    representative of the sensible and latent heat exchanges in the lower
    layers of the troposphere. Estimates of the boundary layer height and
    lifting condensation level are also provided.

    Arguments:
    ---------
    wdir: the working directory path;
    file_list: the list of file containing ts, hus, ps, uas, vas, hfss, te;
    modelname: the name of the model from which the fields are;

    Returns
    -------
    The fields containing boundary layer top height, bondary layer mean
    temperature, temperature at the lifting condensation level (LCL).
    """
    hfss, huss, p_s, t_e, t_s, vv_hor = input_fields(wdir, file_list)
    ricr = RIC_RU
    h_bl = H_U
    ricr = np.where(hfss >= 0.75, ricr, RIC_RS)
    h_bl = np.where(hfss >= 0.75, h_bl, H_S)
    ev_p = huss * p_s / (huss + GAS_CON / RV)  # Water vapour pressure
    td_inv = (1 / T_MELT) - (RV / ALV) * np.log(ev_p / RA_1)  # Dewpoint t.
    t_d = 1 / td_inv
    hlcl = 125.0 * (t_s - t_d)  # Empirical formula for LCL height
    #  Negative heights are replaced by the height of the stable
    #  boundary layer (lower constraint to the height of the cloud layer)
    hlcl = np.where(hlcl >= 0.0, hlcl, h_bl)
    cp_d = GAS_CON / AKAP
    ztlcl = t_s - (G_0 / cp_d) * hlcl
    # Compute the pseudo-adiabatic lapse rate to obtain the height of cloud
    # top knowing emission temperature.
    gw_pa = (G_0 / cp_d) * (
        1
        + ((ALV * huss) / (GAS_CON * ztlcl))
        / (1 + ((ALV**2 * huss * 0.622) / (cp_d * GAS_CON * ztlcl**2)))
    )
    htop = -(t_e - ztlcl) / gw_pa + hlcl
    #  Use potential temperature and critical Richardson number to compute
    #  temperature and height of the boundary layer top
    ths = t_s * (P_0 / p_s) ** AKAP
    thz = ths + 0.03 * ricr * (vv_hor) ** 2 / h_bl
    p_z = p_s * np.exp((-G_0 * h_bl) / (GAS_CON * t_s))  # Barometric eq.
    t_z = thz * (P_0 / p_z) ** (-AKAP)
    outlist = [ztlcl, t_z, htop]
    htop_file, tabl_file, tlcl_file = write_output(
        wdir, modelname, file_list, outlist
    )
    return htop_file, tabl_file, tlcl_file


def mon_from_day(wdir, model, name, filein):
    """Compute monthly mean from daily mean.

    Arguments:
    ---------
    wdir: the working directory path;
    model: the model name;
    name: the name of the field to be averaged;
    filein: the input file containing the field to be averaged;

    Returns
    -------
    The name of the file containing the monthly averaged field.
    """
    cdo = Cdo()
    fileaux = wdir + "/aux.nc"
    cdo.selvar(name, input=filein, output=fileaux)
    move(fileaux, filein)
    fileout = wdir + f"/{model}_{name}_mm.nc"
    cdo.selvar(
        name, input=f"-monmean {filein}", option="-b F32", output=fileout
    )
    return fileout


def removeif(filename):
    """Remove filename if it exists."""
    try:
        os.remove(filename)
    except OSError:
        pass


def wfluxes(model, wdir, input_data):
    """Compute evaporation and rainfall mass fluxes.

    Arguments:
    ---------
    model: the model name;
    wdir: the working directory where the outputs are stored;
    input_data: a dictionary of file names containing the input fields;

    Returns
    -------
    The names of the files containing evaporation and rainfall precipitation
    fluxes.
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
    aux_file = wdir + "/aux.nc"
    evspsbl_file = wdir + f"/{model}_evspsbl.nc"
    cdo.divc(str(L_C), input=f"{hfls_file}", output=evspsbl_file)
    prr_file = wdir + f"/{model}_prr.nc"
    cdo.sub(input=f"{pr_file} {prsn_file}", output=aux_file)
    cdo.chname("pr,prr", input=aux_file, output=prr_file)
    return evspsbl_file, prr_file


def write_output(wdir, model, file_list, varlist):
    """Write auxiliary variables to new NC files, write new attributes.

    Arguments:
    ---------
    wdir: the work directory where the outputs are stored;
    model: the name of the model;
    file_list: the list containing the input fields;
    varlist: a list containing the variables to be written to NC files, i.e.
             tlcl (the temperature at the LCL), t_z (the temperature at the
             boundary layer top), htop (the height of the boundary layer top);
             their dimensions are as (time, lat, lon);

    Returns
    -------
    The names of the files containing fields of boundary layer top height,
    bondary layer mean temperature, temperature at the lifting condensation
    level (LCL).
    """
    cdo = Cdo()
    fourc = fourier_coefficients

    dataset = Dataset(file_list[0])
    ztlcl = varlist[0]
    t_z = varlist[1]
    htop = varlist[2]
    tlcl_temp = wdir + "/tlcl.nc"
    removeif(tlcl_temp)
    with Dataset(tlcl_temp, "w", format="NETCDF4") as w_nc_fid:
        w_nc_fid.description = (
            f"Monthly mean LCL temperature from {model} model. ",
            "Calculated by Thermodynamics model diagnostics ",
            "in ESMValTool. Author Valerio Lembo, ",
            "Meteorologisches Institut, Universitaet ",
            "Hamburg.",
        )
        with Dataset(file_list[0]) as dataset:
            fourc.extr_time(dataset, w_nc_fid)
            fourc.extr_lat(dataset, w_nc_fid, "lat")
            fourc.extr_lon(dataset, w_nc_fid)
        w_nc_var = w_nc_fid.createVariable(
            "tlcl", "f8", ("time", "lat", "lon")
        )
        w_nc_var.setncatts(
            {
                "long_name": "LCL Temperature",
                "units": "K",
                "level_desc": "surface",
                "var_desc": (
                    "LCL temperature from LCL ",
                    "height (Magnus formulas and dry ",
                    "adiabatic lapse ratio)",
                ),
                "statistic": "monthly mean",
            }
        )
        w_nc_fid.variables["tlcl"][:] = ztlcl
    tabl_temp = wdir + "/tabl.nc"
    removeif(tabl_temp)
    with Dataset(tabl_temp, "w", format="NETCDF4") as w_nc_fid:
        w_nc_fid.description = (
            f"Monthly mean BL top temperature for {model} model. ",
            "Calculated by Thermodynamics model diagnostics ",
            "in ESMValTool. Author Valerio ",
            "Lembo, Meteorologisches Institut, ",
            "Universitaet Hamburg.",
        )
        with Dataset(file_list[0]) as dataset_tabl:
            fourc.extr_time(dataset_tabl, w_nc_fid)
            fourc.extr_lat(dataset_tabl, w_nc_fid, "lat")
            fourc.extr_lon(dataset_tabl, w_nc_fid)
        w_nc_var = w_nc_fid.createVariable(
            "tabl", "f8", ("time", "lat", "lon")
        )
        w_nc_var.setncatts(
            {
                "long_name": "Temperature at BL top",
                "units": "K",
                "level_desc": "surface",
                "var_desc": (
                    "Temperature at the Boundary Layer ",
                    "top, from boundary layer thickness and ",
                    "barometric equation",
                ),
                "statistic": "monthly mean",
            }
        )
        w_nc_fid.variables["tabl"][:] = t_z
    htop_temp = wdir + "/htop.nc"
    removeif(htop_temp)
    with Dataset(htop_temp, "w", format="NETCDF4") as w_nc_fid:
        w_nc_fid.description = (
            f"Monthly mean height of the BL top for {model} model. ",
            "Calculated by Thermodynamics model diagnostics ",
            "in ESMValTool. Author Valerio ",
            "Lembo, Meteorologisches Institut, ",
            "Universitaet Hamburg.",
        )
        with Dataset(file_list[0]) as dataset_htop:
            fourc.extr_time(dataset_htop, w_nc_fid)
            fourc.extr_lat(dataset_htop, w_nc_fid, "lat")
            fourc.extr_lon(dataset_htop, w_nc_fid)
        w_nc_var = w_nc_fid.createVariable(
            "htop", "f8", ("time", "lat", "lon")
        )
        w_nc_var.setncatts(
            {
                "long_name": "Height at BL top",
                "units": "m",
                "level_desc": "surface",
                "var_desc": (
                    "Height at the Boundary Layer top, ",
                    "from boundary layer thickness and ",
                    "barometric equation",
                ),
                "statistic": "monthly mean",
            }
        )
        w_nc_fid.variables["htop"][:] = htop
    tlcl_file = wdir + f"/{model}_tlcl.nc"
    cdo.setrtomiss("400,1e36", input=tlcl_temp, output=tlcl_file)
    tabl_file = wdir + f"/{model}_tabl.nc"
    cdo.setrtomiss("400,1e36", input=tabl_temp, output=tabl_file)
    htop_file = wdir + f"/{model}_htop.nc"
    cdo.setrtomiss("12000,1e36", input=htop_temp, output=htop_file)
    return htop_file, tabl_file, tlcl_file
