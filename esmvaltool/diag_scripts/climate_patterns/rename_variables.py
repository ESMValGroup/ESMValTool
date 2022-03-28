def rename_clim_variables(cube):
    # rename variables to fit in JULES framework
    if cube.var_name == "tas":
        cube.rename("t1p5m_clim")
        cube.var_name = "t1p5m_clim"
    if cube.var_name == "hurs":
        cube.rename("rh1p5m_clim")
        cube.var_name = "rh1p5m_clim"
    if cube.var_name == "huss":
        cube.rename("q1p5m_clim")
        cube.var_name = "q1p5m_clim"
    if cube.var_name == "pr":
        cube.rename("precip_clim")
        cube.var_name = "precip_clim"
    if cube.var_name == "sfcWind":
        cube.rename("wind_clim")
        cube.var_name = "wind_clim"
    if cube.var_name == "ps":
        cube.rename("pstar_clim")
        cube.var_name = "pstar_clim"
    if cube.var_name == "rsds":
        cube.rename("swdown_clim")
        cube.var_name = "swdown_clim"
    if cube.var_name == "rlds":
        cube.rename("lwdown_clim")
        cube.var_name = "lwdown_clim"

    return cube


def rename_anom_variables(cube):
    # rename variables to fit in JULES framework
    if cube.var_name == "t1p5m":
        cube.rename("t1p5m_anom")
        cube.var_name = "t1p5m_anom"
    if cube.var_name == "rh1p5m":
        cube.rename("rh1p5m_anom")
        cube.var_name = "rh1p5m_anom"
    if cube.var_name == "q1p5m":
        cube.rename("q1p5m_anom")
        cube.var_name = "q1p5m_anom"
    if cube.var_name == "precip":
        cube.rename("precip_anom")
        cube.var_name = "precip_anom"
    if cube.var_name == "wind":
        cube.rename("wind_anom")
        cube.var_name = "wind_anom"
    if cube.var_name == "pstar":
        cube.rename("pstar_anom")
        cube.var_name = "pstar_anom"
    if cube.var_name == "swdown":
        cube.rename("swdown_anom")
        cube.var_name = "swdown_anom"
    if cube.var_name == "lwdown":
        cube.rename("lwdown_anom")
        cube.var_name = "lwdown_anom"

    return cube


def rename_variables(cube):
    # rename variables to fit in JULES framework
    if cube.var_name == "tas":
        cube.rename("t1p5m")
        cube.var_name = "t1p5m"
    if cube.var_name == "hurs":
        cube.rename("rh1p5m")
        cube.var_name = "rh1p5m"
    if cube.var_name == "huss":
        cube.rename("q1p5m")
        cube.var_name = "q1p5m"
    if cube.var_name == "pr":
        cube.rename("precip")
        cube.var_name = "precip"
    if cube.var_name == "sfcWind":
        cube.rename("wind")
        cube.var_name = "wind"
    if cube.var_name == "ps":
        cube.rename("pstar")
        cube.var_name = "pstar"
    if cube.var_name == "rsds":
        cube.rename("swdown")
        cube.var_name = "swdown"
    if cube.var_name == "rlds":
        cube.rename("lwdown")
        cube.var_name = "lwdown"

    return cube
