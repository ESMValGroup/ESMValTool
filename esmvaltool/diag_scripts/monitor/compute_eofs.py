import logging
import os

import iris
import matplotlib.pyplot as plt
from eofs.iris import Eof
from mapgenerator.plotting.plotmap import PlotMap
from monitor_base import MonitorBase

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(__name__)


class Eofs(MonitorBase):
    def __init__(self, config):
        super().__init__(config)

    def compute(self):
        for module in ['matplotlib', 'fiona']:
            module_logger = logging.getLogger(module)
            module_logger.setLevel(logging.WARNING)
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            variables = group_metadata(data[alias], 'variable_group')
            for var_name, var_info in variables.items():
                logger.info('Plotting variable %s', var_name)
                var_info = var_info[0]
                cube = iris.load_cube(var_info['filename'])
                solver = Eof(cube, weights='coslat')
                variable_options = self._get_variable_options(
                    var_info['variable_group'], '')
                plot_map = PlotMap(loglevel='INFO')
                eof = solver.eofs(neofs=1)[0, ...]
                eof.long_name = var_info.get('eof_name', eof.long_name)
                eof.standard_name = None
                plot_map.plot_cube(eof, save=False, **variable_options)
                plt.savefig(self.get_plot_path('eof', var_info, 'png'),
                            bbox_inches='tight',
                            pad_inches=.2,
                            dpi=plot_map.dpi)
                plt.close(plt.gcf())
                filename = self.get_plot_path('pc', var_info)
                pc = solver.pcs(npcs=1)[:, 0]
                pc.long_name = var_info.get('pc_name', pc.long_name)
                pc.standard_name = None
                self.plot_cube(pc, filename)


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Eofs(config).compute()


if __name__ == "__main__":
    main()
