#!/bin/bash

#loop for Z500 files preparation
#interpolation on regolar 2.5x2.5 grid, NH selection, daily averages.
cdo=/usr/bin/cdo
cdonc="$cdo -f nc"
cdo4="$cdo -f nc4 -z zip"

#define experiment and years
exp=$1
year1=$2
year2=$3
infile=$4
z500filename=$5

DATADIR=$(dirname $z500filename)
TEMPDIR=$DATADIR/tempdir_${exp}_$RANDOM
mkdir -p $TEMPDIR

if [ ! -f $z500filename ] ; then

	echo "Z500 data are missing... full extraction is performed"

	#create a single huge file: not efficient but universal
#	$cdonc cat $INDIR/*.nc $TEMPDIR/fullfile.nc
	#$cdonc sellonlatbox,0,360,0,90 -remapcon2,r144x73 -setlevel,50000 -setname,zg -selyear,$year1/$year2 $TEMPDIR/fullfile.nc $TEMPDIR/smallfile.nc
	$cdonc sellonlatbox,0,360,0,90 -remapcon2,r144x73 -setlevel,50000 -setname,zg $infile $TEMPDIR/smallfile.nc

	#in order to avoid issues, all data are forced to be geopotential height in case geopotential is identified (i.e. values too large for a Z500
	sanityvalue=$($cdonc outputint -fldmean -seltimestep,1 $TEMPDIR/smallfile.nc)
	echo $sanityvalue

	#sanity check
	if [[ $sanityvalue -gt 10000 ]] ; then 
		echo "Values too high: considering as geopotential, dividing by g"
		$cdonc divc,9.80665 $TEMPDIR/smallfile.nc $TEMPDIR/smallfile2.nc
		mv $TEMPDIR/smallfile2.nc $TEMPDIR/smallfile.nc
	else
		echo "Geopotential height identified."
	fi	

	#pippo=($(cdo showyear $TEMPDIR/smallfile.nc) )
	#echo ${pippo[0]}
	#echo ${pippo[-1]}
	$cdo4 -a copy $TEMPDIR/smallfile.nc $z500filename

else
	echo "Z500 NetCDF data seems there, avoid z500_prepare.sh"
fi

#check cleaning
rm -f $TEMPDIR/*.nc
rmdir $TEMPDIR


