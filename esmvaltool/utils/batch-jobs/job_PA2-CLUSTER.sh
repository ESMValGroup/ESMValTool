#!/bin/bash -e
###############################################################################
### BATCH SCRIPT TO RUN THE ESMVALTOOL AT THE DLR PA2-CLUSTER
### Author: Mattia Righi (DLR)
###############################################################################
################# shell to use
#PBS -S /bin/sh
################# export all  environment  variables to job-script
#PBS -V
################# name of the log file
#PBS -o ./$PBS_JOBNAME.$PBS_JOBID.log
################# join standard and error stream (oe, eo) ?
#PBS -j oe
################# do not rerun job if system failure occurs
#PBS -r n
################# send e-mail when [(a)borting|(b)eginning|(e)nding] job
#PBS -m ae
#PBS -M your.email@here
################# ressources (nodes, optional: number of cores; max. runtime)
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
###############################################################################
 
# Submit job with: qsub -q <queue> job_PA2-CLUSTER.sh

# Input arguments
RECIPE=recipe_perfmetrics_CMIP5.yml
CONFIG=config-user.yml

# Set environment
CONDAPATH=  # e.g. /home/soft/miniconda3/
CONDAENV=   # e.g. $CONDAPATH/envs/esmvaltool/bin
ESMVALPATH= # e.g. /home/ESMValTool/esmvaltool

# Changes below this line should not be required
export PATH=$PATH:$CONDAPATH/bin/
conda info --envs
module load ncl
$CONDAENV/esmvaltool $ESMVALPATH/recipes/$RECIPE -c $ESMVALPATH/$CONFIG
