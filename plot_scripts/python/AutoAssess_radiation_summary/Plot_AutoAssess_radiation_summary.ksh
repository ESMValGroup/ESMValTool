#!/usr/bin/ksh
# This is a plotting script for AutoAssess radiation summary plot.
# Preparation for input file: 
# summary files are in directories work/'variable-name' (e.g. work/rlut) under ESMValTool.
# To create a single summary file, do following under 'work/AutoAssess_radiation_rms_summary' directory.
# %cat */summary_global_'model'.csv > summary_global_'model'_all.csv (e.g. 'model'=MPI-ESM-LR)
# Input file should contain error values from more than one variable. Otherwise this program fails.

area=radiation
region=summary_global
cntl=HadGEM2-A
expt=MPI-ESM-LR
obs=obs2
period=ANN
metrics=~/esmvaltool/work/AutoAssess_radiation_rms_summary/$region
ometrics=~/esmvaltool/plot_scripts/python/AutoAssess_radiation_summary/$area
outdir=~/esmvaltool/work/AutoAssess_radiation_rms_summary

plot_norm_ac=~/esmvaltool/plot_scripts/python/AutoAssess_radiation_summary/plot_norm_ac_ideal.py

python2.7 $plot_norm_ac --plot=${outdir}/AutoAssess_radiation_summary.png                              \
                        --exp=${expt} --ctl=${cntl}                  \
                        --file_obs=${ometrics}_${obs}_cmipname_${period}.csv            \
                        --file_exp=${metrics}_${expt}_all.csv   \
                        --file_ctl=${metrics}_${cntl}_all.csv

exit
