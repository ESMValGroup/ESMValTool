"""
Diagnostic to compute and plot the first EOF of an arbitrary input.

It is part of the monitoring recipe.
"""
import logging

import iris
from eofs.iris import Eof
import matplotlib.pyplot as plt
from mapgenerator.plotting.plotmap import PlotMap

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(__name__)


class Eofs(MonitorBase):
    """
    Diagnostic to compute EOFs and plot them.

    It is also an example on how to derive the monitor class to use its
    plotting capabilities in diagnostics that can not be done using only the
    preprocessor.
    """

    def compute(self):
        """Compute the diagnostic."""
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
                pcomp = solver.pcs(npcs=1, pcscaling=1)[:, 0]
                pcomp.long_name = var_info.get('pc_name', pcomp.long_name)
                pcomp.standard_name = None
                self.plot_cube(pcomp, filename)


def main():
    """Main method."""
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Eofs(config).compute()


if __name__ == "__main__":
    main()
