#!/bin/bash -e
###############################################################################
### BATCH SCRIPT TO RUN THE ESMVALTOOL AT DKRZ MISTRAL
### Author: Mattia Righi (DLR)
###############################################################################
#SBATCH --partition=compute
#SBATCH --time=5:00:00
#SBATCH --mail-type=FAIL
#SBATCH --account=bd1083
#SBATCH --output=job_%j.out.log
#SBATCH --error=job_%j.err.log
###############################################################################

# Submit job with: sbatch job_DKRZ-MISTRAL.sh

# Input arguments
RECIPE=/work/bd0854/b380971/Arctic-midlat/2023/github_review_april/ESMValTool/esmvaltool/recipes/recipe_galytska23jgr_test.yml
CONFIG=/home/b/b380971/.esmvaltool/config-user_EG_levante.yml

# Set environment
CONDAPATH=/work/bd0854/b380971/python/mambaforge

# Changes below this line should not be required
source $CONDAPATH/etc/profile.d/conda.sh
conda activate esmvaltool_apr
conda info --envs

esmvaltool run --config-file $CONFIG $RECIPE
