#!/bin/bash -e
###############################################################################
### BATCH SCRIPT TO RUN THE ESMVALTOOL AT DKRZ MISTRAL
### Author: Mattia Righi (DLR)
###############################################################################
#SBATCH --partition=prepost
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --account=XXXXXX
#SBATCH --output=job_%j.out.log
#SBATCH --error=job_%j.err.log
###############################################################################

# Submit job with: sbatch job_DKRZ-MISTRAL.sh

# Input arguments
RECIPE= # e.g. recipe_perfmetrics_CMIP5.yml
CONFIG= # e.g. config-user.yml

# Set environment
CONDAPATH=  # e.g. /home/soft/miniconda3/

# Changes below this line should not be required
source $CONDAPATH/etc/profile.d/conda.sh
conda activate esmvaltool
conda info --envs

esmvaltool run --config-file $CONFIG $RECIPE
