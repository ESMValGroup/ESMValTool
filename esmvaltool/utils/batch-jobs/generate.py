"""Generate SLURM run scripts to run recipes."""
import os
from pathlib import Path

import esmvaltool

env = 'coretool26rc4'
mail = False
submit = False
account = ''  # Select a compute project to be billed
outputs = 'output_rc4'
partition = 'compute'

SPECIAL_RECIPES = {
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
    'recipe_schlund20esd': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
        'memory': '#SBATCH --constraint=512G \n',
    },
    'recipe_collins13ipcc': {
        'partition': '#SBATCH --partition=compute \n',
        'time': '#SBATCH --time=08:00:00 \n',
        'memory': '#SBATCH --constraint=512G \n',
    },
    'recipe_smpi': {
        'partition': '#SBATCH --partition=interactive \n',
        'time': '#SBATCH --time=04:00:00 \n',
        'memory': '#SBATCH --ntasks=8 \n',
    },
    'recipe_smpi_4cds': {
        'partition': '#SBATCH --partition=interactive \n',
        'time': '#SBATCH --time=04:00:00 \n',
        'memory': '#SBATCH --ntasks=8 \n',
    },
}

home = os.path.expanduser('~')

exclude = []  # Fill the list with the names of the recipes to be excluded

dir_recipes = Path('/'.join((esmvaltool.__path__[0], 'recipes')))

# Set the path to conda
conda_path = 'PATH_TO/mambaforge/etc/profile.d/conda.sh'

for recipe in Path(dir_recipes).rglob('*.yml'):
    filename = f'launch_{recipe.stem}.sh'
    if recipe.stem in exclude:
        continue
    with open(f'{filename}', 'w', encoding='utf-8') as file:
        file.write('#!/bin/bash -l \n')
        file.write('\n')
        file.write(f'#SBATCH --job-name={recipe.stem}.%J\n')
        file.write(f'#SBATCH --output={home}/{outputs}/{recipe.stem}.%J.out\n')
        file.write(f'#SBATCH --error={home}/{outputs}/{recipe.stem}.%J.err\n')
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
        os.system(f'sbatch {filename}')
