"""Script containing cube re-naming functions for driving scripts.

Author
------
Gregory Munday (Met Office, UK)
"""


def rename_clim_variables(cube):
    """Rename variables and a coord to fit in JULES framework.

    Parameters
    ----------
    cube : cube
        input cube

    Returns
    -------
    cube : cube
        cube with renamed variables
    """
    if cube.var_name == "tas":
        cube.rename("Air Temperature")
        cube.var_name = "tl1_clim"
    if cube.var_name == "range_tl1":
        cube.rename("Diurnal Range")
        cube.var_name = "range_tl1_clim"
    if cube.var_name == "hurs":
        cube.rename("Relative Humidity")
        cube.var_name = "rh1p5m_clim"
    if cube.var_name == "huss":
        cube.rename("Specific Humidity")
        cube.var_name = "ql1_clim"
    if cube.var_name == "pr":
        cube.rename("Precipitation")
        cube.var_name = "precip_clim"
    if cube.var_name == "sfcWind":
        cube.rename("Wind Speed")
        cube.var_name = "wind_clim"
    if cube.var_name == "ps":
        cube.rename("Surface Pressure")
        cube.var_name = "pstar_clim"
    if cube.var_name == "rsds":
        cube.rename("Surface Downwelling Shortwave Radiation")
        cube.var_name = "swdown_clim"
    if cube.var_name == "rlds":
        cube.rename("Surface Downwelling Longwave Radiation")
        cube.var_name = "lwdown_clim"

    cube.coord("month_number").rename("imogen_drive")

    return cube


def rename_anom_variables(cube):
    """Rename variables and a coord to fit in JULES framework.

    Parameters
    ----------
    cube : cube
        input cube

    Returns
    -------
    cube : cube
        cube with renamed variables
    """
    if cube.var_name == "tas":
        cube.rename("Air Temperature")
        cube.var_name = "tl1_anom"
    if cube.var_name == "range_tl1":
        cube.rename("Diurnal Range")
        cube.var_name = "range_tl1_anom"
    if cube.var_name == "hurs":
        cube.rename("Relative Humidity")
        cube.var_name = "rh1p5m_anom"
    if cube.var_name == "huss":
        cube.rename("Specific Humidity")
        cube.var_name = "ql1_anom"
    if cube.var_name == "pr":
        cube.rename("Precipitation")
        cube.var_name = "precip_anom"
    if cube.var_name == "sfcWind":
        cube.rename("Wind Speed")
        cube.var_name = "wind_anom"
    if cube.var_name == "ps":
        cube.rename("Surface Pressure")
        cube.var_name = "pstar_anom"
    if cube.var_name == "rsds":
        cube.rename("Surface Downwelling Shortwave Radiation")
        cube.var_name = "swdown_anom"
    if cube.var_name == "rlds":
        cube.rename("Surface Downwelling Longwave Radiation")
        cube.var_name = "lwdown_anom"

    cube.coord("month_number").rename("imogen_drive")

    return cube


def rename_variables(cube):
    """Rename variables and a coord to fit in JULES framework.

    Parameters
    ----------
    cube : cube
        input cube

    Returns
    -------
    cube : cube
        cube with renamed variables
    """
    if cube.var_name == "tas":
        cube.rename("Air Temperature")
        cube.var_name = "tl1"
    if cube.var_name == "hurs":
        cube.rename("Relative Humidity")
        cube.var_name = "rh1p5m"
    if cube.var_name == "huss":
        cube.rename("Specific Humidity")
        cube.var_name = "ql1"
    if cube.var_name == "pr":
        cube.rename("Precipitation")
        cube.var_name = "precip"
    if cube.var_name == "sfcWind":
        cube.rename("Wind Speed")
        cube.var_name = "wind"
    if cube.var_name == "ps":
        cube.rename("Surface Pressure")
        cube.var_name = "pstar"
    if cube.var_name == "rsds":
        cube.rename("Surface Downwelling Shortwave Radiation")
        cube.var_name = "swdown"
    if cube.var_name == "rlds":
        cube.rename("Surface Downwelling Longwave Radiation")
        cube.var_name = "lwdown"

    cube.coord("month_number").rename("imogen_drive")

    return cube


def rename_regression_variables(cube):
    """Rename variables to fit in JULES framework.

    Parameters
    ----------
    cube : cube
        input cube

    Returns
    -------
    cube : cube
        cube with renamed variables
    """
    if cube.var_name == "tl1_anom":
        cube.rename("Air Temperature")
        cube.var_name = "tl1_patt"
    if cube.var_name == "range_tl1_anom":
        cube.rename("Diurnal Range")
        cube.var_name = "range_tl1_patt"
    if cube.var_name == "rh1p5m_anom":
        cube.rename("Relative Humidity")
        cube.var_name = "rh1p5m_patt"
    if cube.var_name == "ql1_anom":
        cube.rename("Specific Humidity")
        cube.var_name = "ql1_patt"
    if cube.var_name == "precip_anom":
        cube.rename("Precipitation")
        cube.var_name = "precip_patt"
    if cube.var_name == "wind_anom":
        cube.rename("Wind Speed")
        cube.var_name = "wind_patt"
    if cube.var_name == "pstar_anom":
        cube.rename("Surface Pressure")
        cube.var_name = "pstar_patt"
    if cube.var_name == "swdown_anom":
        cube.rename("Surface Downwelling Shortwave Radiation")
        cube.var_name = "swdown_patt"
    if cube.var_name == "lwdown_anom":
        cube.rename("Surface Downwelling Longwave Radiation")
        cube.var_name = "lwdown_patt"

    return cube


def rename_variables_base(cube):
    """Rename variables and a coord for imogen_mode='off'.

    Parameters
    ----------
    cube : cube
        input cube

    Returns
    -------
    cube : cube
        cube with renamed variables
    """
    if cube.var_name == "tl1_patt":
        cube.rename("Air Temperature")
        cube.var_name = "tas"
    if cube.var_name == "range_tl1_patt":
        cube.rename("Diurnal Range")
        cube.var_name = "tas_range"
    if cube.var_name == "rh1p5m_patt":
        cube.rename("Relative Humidity")
        cube.var_name = "hurs"
    if cube.var_name == "ql1_patt":
        cube.rename("Specific Humidity")
        cube.var_name = "huss"
    if cube.var_name == "precip_patt":
        cube.rename("Precipitation")
        cube.var_name = "pr"
    if cube.var_name == "wind_patt":
        cube.rename("Wind Speed")
        cube.var_name = "sfcWind"
    if cube.var_name == "pstar_patt":
        cube.rename("Surface Pressure")
        cube.var_name = "ps"
    if cube.var_name == "swdown_patt":
        cube.rename("Surface Downwelling Shortwave Radiation")
        cube.var_name = "rsds"
    if cube.var_name == "lwdown_patt":
        cube.rename("Surface Downwelling Longwave Radiation")
        cube.var_name = "rlds"

    cube.coord("imogen_drive").rename("month_number")

    return cube
