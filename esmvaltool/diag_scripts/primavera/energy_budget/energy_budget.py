import os
import logging

import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox

from PIL import Image

import iris
import iris.cube
import iris.analysis
import iris.util

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


class EnergyBudget(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.template = self.cfg.get('plot_template')
        self.start = None
        self.end = None
        self.project = []
        self.exp = []
        self.dataset = []
        self.ensemble = []

        self.pos = {}
        self.pos['rsdt'] = (460, 230)
        self.pos['rsut'] = (15, 300)
        self.pos['rsds'] = (None, None)
        self.pos['rsns'] = (340, 540)
        self.pos['rlut'] = (900, 160)
        self.pos['rlds'] = (1150, 800)
        self.pos['rlns'] = (None, None)
        self.pos['hfss'] = (465, 640)
        self.pos['hfls'] = (690, 660)
        self.pos['up_sw_rfl_surf'] = (170, 830)
        self.pos['sw_rfl_clouds'] = (180, 480)
        self.pos['sw_abs_atm'] = (580, 300)
        self.pos['up_lw_emit_surf'] = (820, 910)
        self.pos['net_surf_rad'] = (None, None)
        self.pos['rad_ads_surface'] = (None, None)
        self.pos['rad_net_toa'] = (350, 62)
        self.pos['bowen_ratio'] = (None, None)

        self.textprops = dict(size=5)
        self.pad = 0
        self.sep = 1
        self.box_alignment = (0, .5)
        self.bboxprops = dict(facecolor='white',
                              boxstyle='round',
                              color='white',
                              alpha=0.9)
        self.dpi = 1000

    def compute(self):
        data = self.load()
        for index, alias in enumerate(data['rsdt']):
            # Shortwave
            data['up_sw_rfl_surf'][alias] = (data['rsds'][alias] -
                                             data['rsns'][alias])
            long_name = 'Upward Shortwave Reflected Surface'
            data['up_sw_rfl_surf'][alias].long_name = long_name

            data['sw_rfl_clouds'][alias] = (data['rsut'][alias] -
                                            data['up_sw_rfl_surf'][alias])
            long_name = 'Shortwave Reflected Clouds'
            data['sw_rfl_clouds'][alias].long_name = long_name

            data['sw_abs_atm'][alias] = (data['rsdt'][alias] -
                                         data['sw_rfl_clouds'][alias] -
                                         data['rsds'][alias])
            long_name = 'Shortwave Absorbed Atmosphere'
            data['sw_abs_atm'][alias].long_name = long_name

            # Longwave
            data['up_lw_emit_surf'][alias] = (data['rlds'][alias] -
                                              data['rlns'][alias])
            long_name = 'Upward Longwave Emitted Surface'
            data['up_lw_emit_surf'][alias].long_name = long_name

            # Net
            data['net_surf_rad'][alias] = (data['rsns'][alias] +
                                           data['rlns'][alias])
            long_name = 'Net Surface Radiation'
            data['net_surf_rad'][alias].long_name = long_name

            # Surface fluxes
            data['rad_ads_surface'][alias] = (data['net_surf_rad'][alias] -
                                              data['hfss'][alias] -
                                              data['hfls'][alias])
            long_name = 'Radiation Adsorbed Surface'
            data['rad_ads_surface'][alias].long_name = long_name

            data['rad_net_toa'][alias] = (data['rsdt'][alias] -
                                          data['rsut'][alias] -
                                          data['rlut'][alias])
            long_name = 'Radiation Net TOA'
            data['rad_net_toa'][alias].long_name = long_name

            data['bowen_ratio'][alias] = (data['hfss'][alias] /
                                          data['hfls'][alias])
            data['bowen_ratio'][alias].long_name = 'Bowen Ratio'

            cubes = iris.cube.CubeList([
                data['rsdt'][alias],
                data['rsut'][alias],
                data['rsds'][alias],
                data['rsns'][alias],
                data['rlut'][alias],
                data['rlds'][alias],
                data['rlns'][alias],
                data['hfss'][alias],
                data['hfls'][alias],
                data['up_sw_rfl_surf'][alias],
                data['sw_rfl_clouds'][alias],
                data['sw_abs_atm'][alias],
                data['up_lw_emit_surf'][alias],
                data['net_surf_rad'][alias],
                data['rad_ads_surface'][alias],
                data['rad_net_toa'][alias],
                data['bowen_ratio'][alias]
            ])
            self.save(index, cubes)
        if self.template:
            self.plot(data)

    def load(self):
        data = {}
        data['rsdt'] = {}
        data['rsut'] = {}
        data['rsds'] = {}
        data['rsns'] = {}
        data['rlut'] = {}
        data['rlds'] = {}
        data['rlns'] = {}
        data['hfss'] = {}
        data['hfls'] = {}
        data['up_sw_rfl_surf'] = {}
        data['sw_rfl_clouds'] = {}
        data['sw_abs_atm'] = {}
        data['up_lw_emit_surf'] = {}
        data['net_surf_rad'] = {}
        data['rad_ads_surface'] = {}
        data['rad_net_toa'] = {}
        data['bowen_ratio'] = {}
        for filename in self.filenames:
            alias = self.filenames.get_info('alias', filename)
            short_name = self.filenames.get_info(n.SHORT_NAME, filename)
            self.project.append(self.filenames.get_info('project', filename))
            self.exp.append(self.filenames.get_info('exp', filename))
            self.dataset.append(self.filenames.get_info('dataset', filename))
            self.ensemble.append(self.filenames.get_info('ensemble', filename))
            self.start = self.filenames.get_info('start_year', filename)
            self.end = self.filenames.get_info('end_year', filename)
            cube = iris.load_cube(filename)
            data[short_name][alias] = cube

        return data

    def save(self, index, cubes):
        file_name = ('{project}_{dataset}_{exp}_{ensemble}_{script}_variables'
                     '_{start}_{end}.nc').format(
                         project=self.project[index],
                         dataset=self.dataset[index],
                         exp=self.exp[index],
                         ensemble=self.ensemble[index],
                         script=self.cfg[n.SCRIPT],
                         start=self.start,
                         end=self.end
                         )
        file_path = os.path.join(self.cfg[n.WORK_DIR], file_name)
        iris.save(cubes, file_path)
        ancestors = []
        for filename in self.filenames:
            if self.dataset[index] in filename:
                ancestors.append(filename)

        caption = ("{script} between {start} and {end}"
                   "according to {dataset}").format(
                       script=self.cfg[n.SCRIPT].split('_'),
                       start=self.start,
                       end=self.end,
                       dataset=self.dataset[index]
                   )
        record = {
            'caption': caption,
            'domains': ['global'],
            'autors': ['vanniere_benoit'],
            'references': ['acknow_project'],
            'ancestors': ancestors
            }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(file_path, record)

    def plot(self, data):
        fig, ax = plt.subplots()
        plot_file = os.path.join(os.path.dirname(__file__), self.template)
        img = Image.open(plot_file).convert("RGBA")
        plt.imshow(img)
        for var in data:
            if None in self.pos[var]:
                continue
            text = []
            char = 'A'
            tags = []
            for index, alias in enumerate(data['rsdt']):
                text.append(TextArea(
                    ('{:}. {:.2f}'.format(char, data[var][alias].data)),
                    textprops=self.textprops))
                tags.append('{:}. {:}-{:}'.format(
                    char,
                    self.dataset[index],
                    self.ensemble[index]))
                char = chr(ord(char) + 1)
            texts_vbox = VPacker(children=text, pad=self.pad, sep=self.sep)
            ann = AnnotationBbox(texts_vbox,
                                 self.pos[var],
                                 xycoords=ax.transData,
                                 box_alignment=self.box_alignment,
                                 bboxprops=self.bboxprops)
            ax.add_artist(ann)
            plt.axis('off')
        ax.text(0, 0.95, '\n'.join(tags),
                fontsize=6,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes)
        script = self.cfg[n.SCRIPT]
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        plt_name = '{script}_Trenberth_Diagram.{out_type}'.format(
            script=script,
            out_type=out_type
        )
        fig.savefig(os.path.join(self.cfg[n.PLOT_DIR], plt_name),
                    format='{0}'.format(out_type),
                    dpi=self.dpi)


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        EnergyBudget(config).compute()


if __name__ == "__main__":
    main()
