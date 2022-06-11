import esmvaltool
import esmvalcore
import os
import glob

from pathlib import Path

env = 'coretool2205'
mail = 'saskia.loosveldt@bsc.es'
account = 'bk1088'

SPECIAL_RECIPES = {
    'recipe_bock20jgr_fig_6-7': {
        'partition': '#SBATCH --partition=interactive \n',
        'time': '#SBATCH --time=12:00:00 \n',
        'memory': '#SBATCH --memory=50000 \n',
        },
    'recipe_bock20jgr_fig_8-10': {
        'partition': '#SBATCH --partition=interactive \n',
        'time': '#SBATCH --time=12:00:00 \n',
        'memory': '#SBATCH --memory=50000 \n',
    },
    'recipe_schlund20esd': {
        'partition': '#SBATCH --partition=interactive \n',
        'time': '#SBATCH --time=12:00:00 \n',
        'memory': '#SBATCH --memory=50000 \n',        
    },
    'recipe_collins13ipcc': {
        'partition': '#SBATCH --partition=interactive \n',
        'time': '#SBATCH --time=12:00:00 \n',
        'memory': '#SBATCH --memory=50000 \n',
    }
}

dir_recipes = '/'.join((esmvaltool.__path__[0], 'recipes'))

recipes = glob.glob(dir_recipes)

home = os.path.expanduser('~')

for recipe in Path(dir_recipes).rglob('*.yml'):
    filename = 'launch_{recipe.stem}.sh'
    file = open(f'{filename}', 'w')
    file.write('#!/bin/bash -l \n')
    file.write('\n')
    file.write(f'#SBATCH --job-name={recipe.stem}.%J\n')
    file.write(f'#SBATCH --output={home}/outputs/{recipe.stem}.%J.out\n')
    file.write(f'#SBATCH --error={home}/outputs/{recipe.stem}.%J.err\n')
    file.write(f'#SBATCH --account={account}\n')
    if not SPECIAL_RECIPES.get(recipe.stem, None):
       file.write('#SBATCH --partition=compute\n')
       file.write('#SBATCH --time=08:00:00\n')
    else:
        file.write(SPECIAL_RECIPES[recipe.stem]['partition'])
        file.write(SPECIAL_RECIPES[recipe.stem]['time'])
        file.write(SPECIAL_RECIPES[recipe.stem]['memory'])
    if mail:
        file.write(f'#SBATCH --mail-user={mail} \n')
        file.write(f'#SBATCH --mail-type=ALL \n')
    file.write('\n')
    file.write('set -eo pipefail \n')
    file.write('unset PYTHONPATH \n')
    file.write('\n')
    file.write(f'. {home}/mambaforge/etc/profile.d/conda.sh\n')
    file.write(f'conda activate {env}\n')
    file.write(f'\n')
    file.write(f'esmvaltool run {recipe.name}')
    file.close()

    #if not SPECIAL_RECIPES.get(recipe.stem, None):
    os.system(f'sbatch {filename}')
