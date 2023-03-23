"""Generate SLURM run scripts to run recipes."""
import os
import subprocess
from pathlib import Path

import esmvaltool

# Name of the conda environment in which esmvaltool is installed
env = 'coretool26rc4'
# Mail notifications when a submitted job fails or finishes
mail = False
# Name of the conda environment in which esmvaltool is installed
submit = False
# Name of the DKRZ account in which the job will be billed
account = ''  # Select a compute project to be billed
# Name of the directory in which the job outputs files) are saved.
# The outputs will be saved in /home/user/<outputs>
outputs = 'output_rc4'
# Default Levante computing partition used
partition = 'compute'
# Default amount of memory used
memory = '64G'
# Default walltime
time = '04:00:00'
# Full path to the mambaforge/etc/profile.d/conda.sh executable
# Set the path to conda
conda_path = 'PATH_TO/mambaforge/etc/profile.d/conda.sh'
# Full path to config_file
# If none, ~/.esmvaltool/config-user.yml is used
config_file = ''

# List of recipes that require non-default SLURM options set above
SPECIAL_RECIPES = {
    'recipe_anav13jclim': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
    },
    'recipe_bock20jgr_fig_6-7': {
        'partition': '#SBATCH --partition=shared \n',
        'time': '#SBATCH --time=48:00:00 \n',
        'memory': '#SBATCH --mem=50000 \n',
    },
    'recipe_bock20jgr_fig_8-10': {
        'partition': '#SBATCH --partition=shared \n',
        'time': '#SBATCH --time=48:00:00 \n',
        'memory': '#SBATCH --mem=50000 \n',
    },
    'recipe_check_obs': {
        'partition': '#SBATCH --partition=compute \n',
        'memory': '#SBATCH --constraint=512G \n',
    },
    'recipe_climate_change_hotspot': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_collins13ipcc': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
        'memory': '#SBATCH --constraint=512G \n',
    },
    'recipe_eady_growth_rate': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_ecs': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_ecs_constraints': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_extreme_index': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_eyring06jgr': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_eyring13jgr_12': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_flato13ipcc_figures_938_941_cmip6': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_ipccwg1ar6ch3_atmosphere': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_ipccwg1ar6ch3_fig_3_9': {
        'partition': '#SBATCH --partition=shared \n',
        'time': '#SBATCH --time=15:00:00 \n',
        'memory': '#SBATCH --mem=150000 \n',
    },
    'recipe_ipccwg1ar6ch3_fig_3_19': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_ipccwg1ar6ch3_fig_3_42_a': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
        'memory': '#SBATCH --constraint=512G \n',
    },
    'recipe_ipccwg1ar6ch3_fig_3_42_b': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_ipccwg1ar6ch3_fig_3_43': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_lauer22jclim_fig3-4_zonal': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_lauer22jclim_fig5_lifrac': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_meehl20sciadv': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_mpqb_xch4': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_perfmetrics_CMIP5': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_perfmetrics_CMIP5_4cds': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_perfmetrics_land_CMIP5': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_russell18jgr': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_schlund20esd': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
        'memory': '#SBATCH --constraint=512G \n',
    },
    'recipe_schlund20jgr_gpp_abs_rcp85': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_schlund20jgr_gpp_change_1pct': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_schlund20jgr_gpp_change_rcp85': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_smpi': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_smpi_4cds': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_tebaldi21esd': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
    },
    'recipe_wenzel16jclim': {
        'partition': '#SBATCH --partition=compute \n',
    },
    'recipe_wenzel16nat': {
        'partition': '#SBATCH --partition=compute \n',
    },
}

# These recipes either use CMIP3 input data
# (see https://github.com/ESMValGroup/ESMValCore/issues/430)
# and recipes where tasks require the full compute node memory.
ONE_TASK_RECIPES = {
    'recipe_bock20jgr_fig_1-4',
    'recipe_bock20jgr_fig_6-7',
    'recipe_bock20jgr_fig_8-10',
    'recipe_flato13ipcc_figure_96',
    'recipe_flato13ipcc_figures_938_941_cmip3',
    'recipe_ipccwg1ar6ch3_fig_3_9',
    'recipe_ipccwg1ar6ch3_fig_3_42_a',
    'recipe_ipccwg1ar6ch3_fig_3_43',
    'recipe_check_obs'
    'recipe_collins13ipcc',
    'recipe_lauer22jclim_fig3-4_zonal',
    'recipe_lauer22jclim_fig5_lifrac',
    'recipe_smpi',
    'recipe_smpi_4cds',
    'recipe_wenzel14jgr',
}

def generate_submit():
    """Generate and submit scripts."""
    home = os.path.expanduser('~')
    # Fill the list with the names of the recipes to be excluded
    exclude = ['recipe_schlund20jgr_gpp_abs_rcp85',
               'recipe_schlund20jgr_gpp_change_1pct',
               'recipe_schlund20jgr_gpp_change_rcp85']
    dir_recipes = Path('/'.join((esmvaltool.__path__[0], 'recipes')))

    for recipe in Path(dir_recipes).rglob('*.yml'):
        filename = f'launch_{recipe.stem}.sh'
        if recipe.stem in exclude:
            continue
        with open(f'{filename}', 'w', encoding='utf-8') as file:
            file.write('#!/bin/bash -l \n')
            file.write('\n')
            file.write(f'#SBATCH --job-name={recipe.stem}.%J\n')
            file.write(
                f'#SBATCH --output={home}/{outputs}/{recipe.stem}.%J.out\n'
            )
            file.write(
                f'#SBATCH --error={home}/{outputs}/{recipe.stem}.%J.err\n'
            )
            file.write(f'#SBATCH --account={account}\n')
            if not SPECIAL_RECIPES.get(recipe.stem, None):
                # continue
                file.write(f'#SBATCH --partition={partition}\n')
                file.write('#SBATCH --time=08:00:00\n')
            else:
                file.write(SPECIAL_RECIPES[recipe.stem]['partition'])
                file.write(SPECIAL_RECIPES[recipe.stem]['time'])
                file.write(SPECIAL_RECIPES[recipe.stem]['memory'])
            if mail:
                file.write('#SBATCH --mail-type=FAIL,END \n')
            file.write('\n')
            file.write('set -eo pipefail \n')
            file.write('unset PYTHONPATH \n')
            file.write('\n')
            file.write(f'. {conda_path}\n')
            file.write(f'conda activate {env}\n')
            file.write('\n')
            file.write(f'esmvaltool run {str(recipe)}')

        if submit:
            subprocess.check_call(['sbatch', filename])


if __name__ == '__main__':
    generate_submit()
