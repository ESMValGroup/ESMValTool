import os
import logging
import calendar

import numpy as np
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
from matplotlib import colors

import iris
import iris.cube
import iris.analysis
import iris.util
from iris.analysis import SUM
from iris.coords import AuxCoord
import iris.coord_categorisation
from iris.cube import CubeList
import iris.quickplot as qplt

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))

class EnergyBudget(object):
    def __init__(self, config):
        self.cfg = config
        self.datasets = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variables = esmvaltool.diag_scripts.shared.Variables(self.cfg)

    def compute(self):
        cubes = iris.load(self.datasets)
        rsdt = cubes.extract_strict('toa_incoming_shortwave_flux') #incoming_sw_rad_toa_field
        rsut = cubes.extract_strict('toa_outgoing_shortwave_flux') #outgoing_sw_rad_toa_field
        rsds = cubes.extract_strict('surface_downwelling_shortwave_flux_in_air') #total_downward_sw_surface_field
        rsns = cubes.extract_strict('Surface Net downward Shortwave Radiation') #net_downward_sw_surface_field
        rlut = cubes.extract_strict('toa_outgoing_longwave_flux') #olr_field
        rlds = cubes.extract_strict('surface_downwelling_longwave_flux_in_air') #ilr_field
        rlns = cubes.extract_strict('Surface Net downward Longwave Radiation') #net_downward_lw_surface_field
        hfss = cubes.extract_strict('surface_upward_sensible_heat_flux') # sensible_heat_field
        hfls = cubes.extract_strict('surface_upward_latent_heat_flux') # latent_heat_field

        results = {}
        #shortwave
        up_sw_reflected_surf = rsds - rsns
        sw_refl_clouds = rsut - up_sw_reflected_surf
        sw_abs_atm = rsdt - sw_refl_clouds - rsds

        #longwave
        up_lw_emitted_surf = rlds - rlns

        #net
        net_surf_rad = rsns + rlns

        #surface fluxes
        rad_adsorbed_surface = net_surf_rad - hfss -hfls
        rad_net_toa = rsdt - rsut - rlut
        bowen_ratio = hfss / hfls
        print(bowen_ratio)

def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        EnergyBudget(config).compute()


if __name__ == '__main__':
    main()                                                                                                                                                                                                                                                                                          1,8           Top
