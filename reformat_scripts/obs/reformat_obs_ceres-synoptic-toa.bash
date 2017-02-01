#! /bin/bash -eu
 
#;;#############################################################################
#;; REFORMAT SCRIPT FOR CERES SYNOPTIC OBSERVATIONAL DATA (TOA)
#;;#############################################################################
#;;
#;; Tier
#;;    Tier 3: restricted dataset (downloads larger than 2GB require
#;;            registration).
#;;
#;; Source
#;;    http://ceres-tool.larc.nasa.gov/ord-tool/
#;;
#;; Download and processing instructions
#;;    On above link, 
#;;       *) choose SYN1deg
#;;       *) expand "Computed TOA Fluxes" and check
#;;          ') Shortwave Flux (All sky + Clear sky)
#;;          ') Longwave Flux (All sky + Clear sky)
#;;          ') Shortwave Flux Down (All sky + Clear sky)
#;;       *) Download
#;;       *) Update INFILES below
#;;       *) Run script as "./reformat_obs_ceres-synoptic-toa.bash"
#;;
#;; Caveats
#;;    Script only handles one year at the time
#;;
#;; Modification history
#;;    20151002-A_eval_ma: written.
#;;
#;;#############################################################################

#INFILES=(CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040101-20040816.nc CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040817-20041231.nc)
INFILES=(CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040101-20040701.nc CERES_SYN1deg-3H_Terra-Aqua-MODIS_Ed3A_Subset_20040702-20041231.nc)

orig=(toa_comp_sw-up_clr_3h   toa_comp_sw-up_all_3h    toa_comp_lw-up_clr_3h    toa_comp_lw-up_all_3h)
vars=(rsutcs                  rsut                     rlutcs                   rlut)

rm -f tmp1.nc

inpath=/export/pa_data01/ESMVal/obs/RAW/Tier3/CERES-SYN1deg/toa
outpath=/export/pa_data02/ESMVal/obs/Tier3/CERES-SYN1deg

for ifile_idx in $(seq 0 $((${#INFILES[@]}-1)));do
    INFILE=$inpath/${INFILES[ifile_idx]}
    yrs=$(echo $INFILE | sed 's/.*_\([0-9][0-9][0-9][0-9][0-9][0-9]\).*-\([0-9][0-9][0-9][0-9][0-9][0-9]\).*.nc/\1-\2/')
    for idx in $(seq 0 $((${#orig[@]}-1)));do
        echo ncks -v ${orig[$idx]} $INFILE ${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc
        ncks -v ${orig[$idx]} $INFILE tmp1.nc
        ncrename -v ${orig[$idx]},${vars[$idx]} tmp1.nc
        ncks --mk_rec_dmn time tmp1.nc ${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc
        rm -f tmp1.nc
    done
done

yrs1=$(echo ${INFILES[0]} | sed 's/.*_\([0-9][0-9][0-9][0-9][0-9][0-9]\).*-\([0-9][0-9][0-9][0-9][0-9][0-9]\).*.nc/\1-\2/')
yrs_start=$(echo $yrs1 | sed 's/\([0-9][0-9][0-9][0-9][0-9][0-9]\)-\([0-9][0-9][0-9][0-9][0-9][0-9]\)/\1/')
yrs2=$(echo ${INFILES[${#INFILES[@]}-1]} | sed 's/.*_\([0-9][0-9][0-9][0-9][0-9][0-9]\).*-\([0-9][0-9][0-9][0-9][0-9][0-9]\).*.nc/\1-\2/')
yrs_end=$(echo $yrs2 | sed 's/\([0-9][0-9][0-9][0-9][0-9][0-9]\)-\([0-9][0-9][0-9][0-9][0-9][0-9]\)/\2/')
yrs=${yrs_start}-${yrs_end}

echo $yrs1
echo $yrs_start
echo $yrs2
echo $yrs_end
echo $yrs
exit

for idx in $(seq 0 $((${#orig[@]}-1)));do
    ncrcat -h ${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_* ${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc
done

for idx in $(seq 0 $((${#vars[@]}-1)));do
    mv -i ${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc $outpath/${vars[$idx]}_3hr_CERES-SYN1deg_observation_r1i1p1_${yrs}.nc
#    mv -i ${vars[$idx]}_3hr_CERES-SYN1deg-3H-computed_r1i1p1_${yrs}.nc $outpath/OBS_CERES-SYN1deg-3H-computed_obs_1_T2Is_${vars[$idx]}.nc
done
