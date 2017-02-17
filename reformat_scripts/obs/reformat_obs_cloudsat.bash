#!/usr/bin/env bash -eu
###############################################################################
## REFORMAT SCRIPT FOR THE CLOUDSAT OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 1: obs4mips
##
## Source
##    Download at ESGF-nodes, documenation at
##    https://www.earthsystemcog.org/projects/obs4mips/satellite_products
##
## Last access
##
## Download and processing instructions
##    Download from ESGF, expected file name is along the lines, 
##    cltcloudsat_obs4MIPs_CloudSat_L3_V2.0_20060601_20060630.nc
##    Run this script
##
## Caveats
##
## Modification history
##    20151009-A_eval_ma: written.
##
###############################################################################

# USER SETTINGS
declare -A cloudsat2cmip5  # map variable names
cloudsat2cmip5['cltcloudsat']='clt'
# USER SETTINGS ENDS

rm -f tmp1.nc
rm -f *-time-is-record-dimension.nc

function usage {
 echo "Usage: ./reformat_obs_cloudsat.bash [-h] -s <SRC_DIR> -d <DST_DIR> -v <CLOUDSAT_VARIABLE>
" 1>&2 
}

function log {
 echo "[ $(date -u '+%Y-%m-%d  %H:%M') ]: " $*
}

info()
{
    log "*II* $*"
}

warning()
{
    log "*WW* $*"
}

error()
{
    log "*EE* $*" 1>&2
    exit 1
}

function usage_and_exit {
  exit_code=${1:-0}
  usage
  exit $exit_code
}

if [ $# -eq 0 ]; then
  usage_and_exit 0
fi

while getopts "hs:d:v:" opt; do
  case $opt in
    h)
      usage_and_exit 0
      ;;
    v)
      var=${OPTARG}
      ;;
    s)
      src_dir=${OPTARG}
      ;;
    d)
      dst_dir=$OPTARG
      ;;
  esac
done

YYYYMMDD1=$(ls -1 ${src_dir}/${var}_* | sed 's/.*\([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\)_\([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\1/' | sort -n | head -1)
YYYYMMDD2=$(ls -1 ${src_dir}/${var}_* | sed 's/.*\([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\)_\([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\).*/\2/' | sort -n | tail -1)

# Round first year YYYYMMDD1 upwards if needed
month=$(date +%m -d $YYYYMMDD1)
if [[ $month -ne '01' ]]
then
    year=$(date +%Y -d $YYYYMMDD1)
    year=$((year + 1))
    month=01
    day=01
    YYYYMMDD1=${year}${month}${day}
fi

# Round last year YYYYMMDD2 downwards if needed
month=$(date +%m -d $YYYYMMDD2)
if [[ $month -ne '12' ]]
then
    year=$(date +%Y -d $YYYYMMDD2)
    year=$((year - 1))
    month=12
    day=31
    YYYYMMDD2=${year}${month}${day}
fi


for fil in $(ls ${src_dir}/${var}_*)
do
    ncks --mk_rec_dmn time $fil ${fil%*.nc}-time-is-record-dimension.nc
done
ncrcat -h *time-is-record-dimension.nc tmp1.nc
rm -f ${dst_dir}/${cloudsat2cmip5[$var]}_Amon_Cloudsat-L3_obs4mips_v2.0_${YYYYMMDD1}-${YYYYMMDD2}.nc
ncrename -v ${var},${cloudsat2cmip5[$var]} tmp1.nc ${dst_dir}/${cloudsat2cmip5[$var]}_Amon_Cloudsat-L3_obs4mips_v2.0_${YYYYMMDD1}-${YYYYMMDD2}.nc
rm -f tmp1.nc
rm -f *-time-is-record-dimension.nc


# Rename to OBS-class
mv ${dst_dir}/${cloudsat2cmip5[$var]}_Amon_Cloudsat-L3_obs4mips_v2.0_${YYYYMMDD1}-${YYYYMMDD2}.nc ${dst_dir}/OBS_Cloudsat_sat_2_T2Ms_${cloudsat2cmip5[$var]}.nc
