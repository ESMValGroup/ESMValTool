#!/bin/bash -l 

#SBATCH --job-name=recipe_bock20jgr_fig_1-4.%J
#SBATCH --output=/home/b/b380971//work/bd1083/b380971/output/daily_ozone//recipe_bock20jgr_fig_1-4.%J.out
#SBATCH --error=/home/b/b380971//work/bd1083/b380971/output/daily_ozone//recipe_bock20jgr_fig_1-4.%J.err
#SBATCH --account=bd1083
#SBATCH --partition=interactive
#SBATCH --time=04:00:00
#SBATCH --mem=64G

set -eo pipefail 
unset PYTHONPATH 

. /work/bd0854/b380971/python/mambaforge/etc/profile.d/conda.sh
conda activate esmvaltool_sweden

esmvaltool run bock20jgr/recipe_bock20jgr_fig_1-4.yml --max_parallel_tasks=1
