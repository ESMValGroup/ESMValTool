import os

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as names
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic


class DecadalExample
    def __init__(self, config):
        self.cfg = config

    def get_provenance_record(title, ancestor_files):
        """Create a provenance record describing the diagnostic data and plot."""
        caption = ("Comparison of {title} between a DCPP experiment"
                   "and an observational dataset.").format(title)

        record = {
            'caption': caption,
            'statistics': ['mean'],
            'domains': ['global'],
            'plot_types': ['timeseries'],
            'authors': [
                'loosveldt-tomas_saskia',
            ],
            'references': [
                'acknow_project',
            ],
            'ancestors': ancestor_files,
        }
        return record


    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'project')
        time = []
        ancestors = []
        for dataset in data['OBS6']:
            cube = iris.load_cube(dataset['filename'])
            name = dataset['dataset']
            iris.coord_categorisation.add_year(cube, 'time')
            cube.coord('time').bounds = None
            plt.plot(cube.coord('time').points, cube.data, label=f'{name}')
            ancestors.append(dataset['filename']

        for dataset in data['CMIP6']:
            cube = iris.load_cube(dataset['filename'])
            iris.coord_categorisation.add_year(cube, 'time')
            name = dataset['dataset']
            sub_exp = dataset['sub_experiment']
            cube.coord('time').bounds = None
            plt.plot(cube.coord('time').points, cube.data, label=f'{name}-{sub_exp}')
            ancestors.append(dataset['filename'])
        
        plt.rcParams["figure.figsize"] = (40,6)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 5.5})
        plt.xlabel("time (days since 01-01-1850)")
        plt.ylabel("Temperature (K)")
        title = 'Global mean of Near-Surface Air Temperature (tas)'
        plt.title(title)
        plt.tight_layout()
    
        plt.grid(True)


        extension = self.cfg['output_file_type']
        plotname = 'decadal_test' + f'.{extension}'
        plt.savefig(os.path.join(self.cfg[names.PLOT_DIR], plotname), dpi=3000)

        provenance_record = get_provenance_record(title, ancestors)
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plotname, provenance_record)
        plt.close()

        


def main():
    with run_diagnostic() as config:
        DecadalExample(config).compute()

if __name__ == "__main__":
    main()


