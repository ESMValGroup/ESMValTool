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
RECIPE=/work/bd0854/b380971/v2_copied/ESMValTool-private/esmvaltool/recipes/galytska22/recipe_Arctic_telecon_tests.yml
CONFIG=/work/bd0854/b380971/v2_copied/ESMValTool-private/config-user_EG_levante.yml

# Set environment
CONDAPATH=/home/b/b380971/mambaforge

# Changes below this line should not be required
source $CONDAPATH/etc/profile.d/conda.sh
conda activate esmval_aug
conda info --envs

esmvaltool run --config-file $CONFIG $RECIPE
