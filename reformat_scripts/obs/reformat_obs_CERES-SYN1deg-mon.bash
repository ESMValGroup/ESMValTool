#!/usr/bin/env bash -eu
###############################################################################
## REFORMAT SCRIPT FOR CERES SYNOPTIC OBSERVATIONAL DATA (monthly means)
###############################################################################
##
## Tier
##    Tier 2: other freely available datasets (only downloads larger than 2GB
##            require registration).
##
## Source
##    https://eosweb.larc.nasa.gov/project/ceres/syn1deg-month_ed3a_table
##
## Last access
##    2017-06-08
##
## Download and processing instructions
##    On above link, 
##       *) select "computed surface fluxes"
##       *) download data
##       *) select "computed TOA fluxes"
##       *) download data
##       *) ppdate INFILES below
##       *) run script
##
## Caveats
##    File names of the input/output files have to be adjusted manually
##    including specification of the time range contained in the files ($YRS).
##
## Modification history
##    20170609-A_laue_ax: written, based on reformat_obs_CERES-SYN1deg-SFC.bash
##
###############################################################################

inpath="${ESMValTool_RAWOBSPATH}/Tier2/CERES-SYN1deg"
outpath="${ESMValTool_OBSPATH}/Tier2/CERES-SYN1deg"

# **************************************************************************
# ******** change input file names / base of output file names here ********
# **************************************************************************
YRS=200003-201609
SFC_SRCFILE=CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Computed-SFC_${YRS}.nc
TOA_SRCFILE=CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed3A_Computed-TOA_${YRS}.nc
OUTBASE=CERES-SYN1deg_comp_Ed3A
# **************************************************************************

sfc_orig=(sfc_comp_sw-up_all_mon
          sfc_comp_sw-down_all_mon
          sfc_comp_lw-up_all_mon
          sfc_comp_lw-down_all_mon
          sfc_comp_lw-down_clr_mon)

sfc_vars=(rsus
          rsds
          rlus
          rlds
          rldscs)

toa_orig=(toa_comp_sw-up_all_mon
          toa_comp_lw-up_all_mon
          toa_comp_sw-up_clr_mon
          toa_comp_lw-up_clr_mon
          toa_comp_sw-down_all_mon)

toa_vars=(rsut
          rlut
          rsutcs
          rlutcs
          rsdt)

if [[ ! -e $outpath ]]; then
    mkdir -p $outpath
fi
rm -f ${outpath}/tmp1.nc

# process surface variables

INFILE=$inpath/${SFC_SRCFILE}
for idx in $(seq 0 $((${#sfc_orig[@]}-1)));do
    echo ncks -v ${sfc_orig[$idx]} $INFILE ${outpath}/OBS_${OUTBASE}_T2Ms_${sfc_vars[$idx]}_${YRS}.nc
    ncks -v ${sfc_orig[$idx]} $INFILE ${outpath}/tmp1.nc
    ncrename -v ${sfc_orig[$idx]},${sfc_vars[$idx]} ${outpath}/tmp1.nc
    ncks --mk_rec_dmn time ${outpath}/tmp1.nc ${outpath}/OBS_${OUTBASE}_T2Ms_${sfc_vars[$idx]}_${YRS}.nc
    rm -f ${outpath}/tmp1.nc
done

# process TOA variables

INFILE=$inpath/${TOA_SRCFILE}
for idx in $(seq 0 $((${#toa_orig[@]}-1)));do
    echo ncks -v ${toa_orig[$idx]} $INFILE ${outpath}/OBS_${OUTBASE}_T2Ms_${toa_vars[$idx]}_${YRS}.nc
    ncks -v ${toa_orig[$idx]} $INFILE ${outpath}/tmp1.nc
    ncrename -v ${toa_orig[$idx]},${toa_vars[$idx]} ${outpath}/tmp1.nc
    ncks --mk_rec_dmn time ${outpath}/tmp1.nc ${outpath}/OBS_${OUTBASE}_T2Ms_${toa_vars[$idx]}_${YRS}.nc
    rm -f ${outpath}/tmp1.nc
done

