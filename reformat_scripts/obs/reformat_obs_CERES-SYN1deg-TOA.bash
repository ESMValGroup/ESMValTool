#!/usr/bin/env bash -eu
###############################################################################
## REFORMAT SCRIPT FOR CERES SYNOPTIC OBSERVATIONAL DATA (TOA)
###############################################################################
##
## Tier
##    Tier 3: restricted dataset (downloads larger than 2GB require
##            registration).
##
## Source
##    http://ceres-tool.larc.nasa.gov/ord-tool/
##
## Last access
##
## Download and processing instructions
##    On above link, 
##       *) choose SYN1deg
##       *) expand "Computed TOA Fluxes" and check
##          ') Shortwave Flux (All sky + Clear sky)
##          ') Longwave Flux (All sky + Clear sky)
##          ') Shortwave Flux Down (All sky + Clear sky)
##       *) Download
##       *) Update INFILES below
##       *) Run script as "./reformat_obs_ceres-synoptic-toa.bash"
##
## Caveats
##    Script only handles one year at the time
##
## Modification history
##    20151002-A_eval_ma: written.
##
###############################################################################

inpath="${ESMValTool_RAWOBSPATH}/Tier3/CERES-SYN1deg/toa"
outpath="${ESMValTool_OBSPATH}/Tier3/CERES-SYN1deg"

#INFILES=(CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040101-20040816.nc
#         CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040817-20041231.nc)
INFILES=(CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040101-20040701.nc
         CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040702-20041231.nc)

orig=(toa_comp_sw-up_clr_3h
      toa_comp_sw-up_all_3h
      toa_comp_lw-up_clr_3h
      toa_comp_lw-up_all_3h)

vars=(rsutcs
      rsut
      rlutcs
      rlut)

if [[ ! -e $outpath ]]; then
    mkdir -p $outpath
fi
rm -f ${outpath}/tmp1.nc

for ifile_idx in $(seq 0 $((${#INFILES[@]}-1)));do
    INFILE=$inpath/${INFILES[ifile_idx]}
    yrs=$(echo $INFILE | sed 's/.*_\([0-9][0-9][0-9][0-9][0-9][0-9]\).*-\([0-9][0-9][0-9][0-9][0-9][0-9]\).*.nc/\1-\2/')
    for idx in $(seq 0 $((${#orig[@]}-1)));do
        echo ncks -v ${orig[$idx]} $INFILE ${outpath}/${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc
        ncks -v ${orig[$idx]} $INFILE ${outpath}/tmp1.nc
        ncrename -v ${orig[$idx]},${vars[$idx]} ${outpath}/tmp1.nc
        ncks --mk_rec_dmn time ${outpath}/tmp1.nc ${outpath}/${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc
        rm -f ${outpath}/tmp1.nc
    done
done

yrs1=$(echo ${INFILES[0]} | sed 's/.*_\([0-9][0-9][0-9][0-9][0-9][0-9]\).*-\([0-9][0-9][0-9][0-9][0-9][0-9]\).*.nc/\1-\2/')
yrs_start=$(echo $yrs1 | sed 's/\([0-9][0-9][0-9][0-9][0-9][0-9]\)-\([0-9][0-9][0-9][0-9][0-9][0-9]\)/\1/')
yrs2=$(echo ${INFILES[${#INFILES[@]}-1]} | sed 's/.*_\([0-9][0-9][0-9][0-9][0-9][0-9]\).*-\([0-9][0-9][0-9][0-9][0-9][0-9]\).*.nc/\1-\2/')
yrs_end=$(echo $yrs2 | sed 's/\([0-9][0-9][0-9][0-9][0-9][0-9]\)-\([0-9][0-9][0-9][0-9][0-9][0-9]\)/\2/')
yrs=${yrs_start}-${yrs_end}

for idx in $(seq 0 $((${#orig[@]}-1)));do
    ncrcat -h ${outpath}/${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_* ${outpath}/${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc
done

for idx in $(seq 0 $((${#vars[@]}-1)));do
    mv -f ${outpath}/${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc $outpath/OBS_CERES-SYN1deg_sat_3hr_T2Is_${vars[$idx]}_${yrs}.nc
    rm -f ${outpath}/${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_*
done
