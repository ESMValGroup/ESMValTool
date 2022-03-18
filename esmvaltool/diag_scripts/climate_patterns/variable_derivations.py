def diurnal_temp_range(cubelist):
    range_cube = cubelist[0] - cubelist[1]
    range_cube.rename("range_t1p5m_clim")
    range_cube.var_name = "range_t1p5m_clim"

    return range_cube
