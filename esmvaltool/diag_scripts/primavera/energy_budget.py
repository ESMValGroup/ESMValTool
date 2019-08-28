import os
import logging
import string

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.font_manager
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

import iris
import iris.cube
import iris.analysis
import iris.util

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))

class EnergyBudget(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.template = self.cfg.get('plot_template')
        self.start = None
        self.end = None
        self.project = None
        self.exp = None

        self.pos = {}
        self.pos['rsdt'] = (511, 170)
        self.pos['rsut'] = (15, 176)
        self.pos['rsds'] = (None, None)
        self.pos['rsns'] = (349, 632)
        self.pos['rlut'] = (916, 130)
        self.pos['rlds'] = (1052, 833)
        self.pos['rlns'] = (None, None)
        self.pos['hfss'] = (524, 632)
        self.pos['hfls'] = (656, 632)
        self.pos['up_sw_rfl_surf'] = (105, 632)
        self.pos['sw_rfl_clouds'] = (300, 361)
        self.pos['sw_abs_atm'] = (583, 325)
        self.pos['up_lw_emit_surf'] = (825, 833)
        self.pos['net_surf_rad'] = (None, None)
        self.pos['rad_ads_surface'] = (None, None)
        self.pos['rad_net_toa'] = (350, 62)
        self.pos['bowen_ratio'] = (None, None)

        self.textprops = dict(size=6)
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
        for dataset in data['rsdt'].keys():
            #shortwave
            data['up_sw_rfl_surf'][dataset] = (data['rsds'][dataset] -
                                               data['rsns'][dataset])
            long_name = 'Upward Shortwave Reflected Surface'
            data['up_sw_rfl_surf'][dataset].long_name = long_name

            data['sw_rfl_clouds'][dataset] = (data['rsut'][dataset] -
                                              data['up_sw_rfl_surf'][dataset])
            long_name = 'Shortwave Reflected Clouds'
            data['sw_rfl_clouds'][dataset].long_name = long_name

            data['sw_abs_atm'][dataset] = (data['rsdt'][dataset] -
                                           data['sw_rfl_clouds'][dataset] -
                                           data['rsds'][dataset])
            long_name = 'Shortwave Absorbed Atmosphere'
            data['sw_abs_atm'][dataset].long_name = long_name

            #longwave
            data['up_lw_emit_surf'][dataset] = (data['rlds'][dataset] -
                                                data['rlns'][dataset])
            long_name = 'Upward Longwave Emitted Surface'
            data['up_lw_emit_surf'][dataset].long_name = long_name

            #net
            data['net_surf_rad'][dataset] = (data['rsns'][dataset] +
                                             data['rlns'][dataset])
            long_name = 'Net Surface Radiation'
            data['net_surf_rad'][dataset].long_name = long_name

            #surface fluxes
            data['rad_ads_surface'][dataset] = (data['net_surf_rad'][dataset] -
                                                data['hfss'][dataset] -
                                                data['hfls'][dataset])
            long_name = 'Radiation Adsorbed Surface'
            data['rad_ads_surface'][dataset].long_name = long_name

            data['rad_net_toa'][dataset] = (data['rsdt'][dataset] -
                                            data['rsut'][dataset] -
                                            data['rlut'][dataset])
            long_name = 'Radiation Net TOA'
            data['rad_net_toa'][dataset].long_name = long_name

            data['bowen_ratio'][dataset] = (data['hfss'][dataset] /
                                            data['hfls'][dataset])
            data['bowen_ratio'][dataset].long_name = 'Bowen Ratio'

            cubes = iris.cube.CubeList([
                data['rsdt'][dataset],
                data['rsut'][dataset],
                data['rsds'][dataset],
                data['rsns'][dataset],
                data['rlut'][dataset],
                data['rlds'][dataset],
                data['rlns'][dataset],
                data['hfss'][dataset],
                data['hfls'][dataset],
                data['up_sw_rfl_surf'][dataset],
                data['sw_rfl_clouds'][dataset],
                data['sw_abs_atm'][dataset],
                data['up_lw_emit_surf'][dataset],
                data['net_surf_rad'][dataset],
                data['rad_ads_surface'][dataset],
                data['rad_net_toa'][dataset],
                data['bowen_ratio'][dataset]
            ])
            self.save(dataset, cubes)
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
            dataset = self.filenames.get_info(n.DATASET, filename)
            short_name = self.filenames.get_info(n.SHORT_NAME, filename)
            self.project = self.filenames.get_info('project', filename)
            self.exp = self.filenames.get_info('exp', filename)
            self.start = self.filenames.get_info('start_year', filename)
            self.end = self.filenames.get_info('end_year', filename)
            cube = iris.load_cube(filename)
            data[short_name][dataset] = cube

        return data

    def save(self, dataset, cubes):
        file_name = ('{project}_{dataset}_{exp}_{script}_variables'
                    '_{start}_{end}.nc').format(
                        project=self.project,
                        dataset=dataset,
                        exp=self.exp,
                        script=self.cfg[n.SCRIPT],
                        start=self.start,
                        end=self.end
                    )
        iris.save(cubes, os.path.join(self.cfg[n.WORK_DIR], file_name))

    def plot(self, data):
        fig, ax = plt.subplots()
        plot_file = os.path.join(os.path.dirname(__file__), self.template)
        img = Image.open(plot_file).convert("RGBA")
        plt.imshow(img)
        for var in data.keys():
            if None in self.pos[var]:
                continue
            text = []
            char = 'A'
            tags = []
            for dataset in data['rsdt'].keys():
                text.append(TextArea((
                    '{:}. {:.2f}'.format(char,
                                         data[var][dataset].data)),
                                         textprops=self.textprops))
                tags.append('{:}. {:}'.format(char, dataset))
                char = chr(ord(char) + 1)
            texts_vbox = VPacker(children=text, pad=self.pad, sep=self.sep)
            ann = AnnotationBbox(texts_vbox,
                                 self.pos[var],
                                 xycoords=ax.transData,
                                 box_alignment=self.box_alignment,
                                 bboxprops=self.bboxprops)
            ax.add_artist(ann)
            plt.axis('off')
        ax.text(1, 0, '\n'.join(tags), fontsize=10)
        script = self.cfg[n.SCRIPT]
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        models = '_'.join(data['rsdt'].keys())
        plt_name = '{script}_{models}_Trenberth_Diagram.{out_type}'.format(
            script=script,
            models=models,
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
