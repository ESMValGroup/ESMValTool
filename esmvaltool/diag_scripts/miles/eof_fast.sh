#!/bin/bash

# -----------------------------------------------#
# NB: due to change in CDO code after version 1.6.4
# eigenvectors need to be normalized and area weighting must be specified for eofcoeffs
# if you use a previous version please be aware that inconsistencies may arise
# Beware: PCs are not standardized (this is done in the plotting tool)
# -----------------------------------------------#
cdo=/usr/bin/cdo
cdonc="$cdo -f nc"
cdo4="$cdo -f nc4 -z zip"

exp=$1
year1=$2
year2=$3
seasons=$4
teles=$5
z500filename=$6
FILESDIR=$7

DATADIR=$(dirname $z500filename)
TEMPDIR=$DATADIR/tempdir_${exp}_$RANDOM
mkdir -p $TEMPDIR

#number of EOFs
neofs=4

#preparing unique netcdf file
$cdonc monmean -selyear,$year1/$year2 $z500filename $TEMPDIR/monthly_file.nc

for tele in $teles ; do
for season in $seasons ; do
	echo $season

	#fix folders and file names
	EOFDIR=$FILESDIR/$exp/EOFs/${tele}/${year1}_${year2}/${season}
	mkdir -p $EOFDIR
        suffix=${exp}_${year1}_${year2}_${season}

	#select seasons, compute monthly anomalies
        $cdonc selseas,$season $TEMPDIR/monthly_file.nc $TEMPDIR/season_monthly.nc
        $cdonc -r ymonmean $TEMPDIR/season_monthly.nc $TEMPDIR/season_mean.nc
        $cdo4 -r -b 64 ymonsub $TEMPDIR/season_monthly.nc $TEMPDIR/season_mean.nc $EOFDIR/Z500_monthly_anomalies_${suffix}.nc

	#fix borders for EOFs
        if [ "${tele}" = NAO ]; then
                box=-90,40,20,85
        fi
        if [ "${tele}" = AO ]; then
                box=0,360,20,85
        fi

	#select box for anomalies
	$cdonc sellonlatbox,${box}  $EOFDIR/Z500_monthly_anomalies_${suffix}.nc $TEMPDIR/box_anomalies_monthly.nc	

	# compute EOFs
        $cdonc -r eof,$neofs $TEMPDIR/box_anomalies_monthly.nc $EOFDIR/${tele}_Z500_eigenvalues_${suffix}.nc $TEMPDIR/pattern.nc

	# normalize eigenvectors
	$cdonc enlarge,$TEMPDIR/pattern.nc -sqrt -fldsum -sqr $TEMPDIR/pattern.nc $TEMPDIR/factor.nc
	$cdo4 div $TEMPDIR/pattern.nc $TEMPDIR/factor.nc $EOFDIR/${tele}_Z500_pattern_${suffix}.nc

	#------------------------------------------------#
	# Beta version to fix the signs of the main EOFS
	#------------------------------------------------#
	
	# define regions for sign control: boxes where values should be positive
	if [ "${tele}" = NAO ]; then
		eof1=-30,30,40,50 	#NAO
		eof2=-60,0,40,60 	#East Atlantic Pattern 
		eof3=-30,30,50,70	#Scandinavian Blocking
	fi

	if [ "${tele}" = AO ]; then
		eof1=0,360,20,50 	#Arctic Oscillation
		eof2=-120,-60,40,60	#PNA
	fi

	cp $EOFDIR/${tele}_Z500_pattern_${suffix}.nc $TEMPDIR/temp_pattern.nc

	for n in $(seq 1 $neofs) ; do
		
		name=eof$n
		
		# if sign control is not fixed, skip
		if [[ -z ${!name} ]] ; then $cdonc seltimestep,$n $TEMPDIR/temp_pattern.nc $TEMPDIR/new_$n.nc ; continue ; fi
		
		#evaluate the value in the box
		value=$($cdonc output -fldmean -sellonlatbox,${!name} -seltimestep,$n $TEMPDIR/temp_pattern.nc)

		#if negative, flip sign, otherwise ignore
		if (( $(echo "$value < 0" |bc -l) )) ; then
			$cdonc mulc,-1 -seltimestep,$n $TEMPDIR/temp_pattern.nc $TEMPDIR/new_$n.nc
		else
			$cdonc seltimestep,$n $TEMPDIR/temp_pattern.nc $TEMPDIR/new_$n.nc
		fi
	done

	#create new file
	rm -f $TEMPDIR/corrected.nc
	$cdo4 cat $TEMPDIR/new_*.nc $TEMPDIR/corrected.nc	
	mv $TEMPDIR/corrected.nc $EOFDIR/${tele}_Z500_pattern_${suffix}.nc 
	rm $TEMPDIR/new_*.nc

	# evalute grid area weights (implicitely used in eofs) and compute principal components
	$cdonc gridweights $TEMPDIR/box_anomalies_monthly.nc $TEMPDIR/ww.nc
        $cdo4 -r eofcoeff $EOFDIR/${tele}_Z500_pattern_${suffix}.nc -mul $TEMPDIR/ww.nc $TEMPDIR/box_anomalies_monthly.nc $EOFDIR/${tele}_monthly_timeseries_${suffix}_

done
done

rm -f $TEMPDIR/*.nc
rmdir $TEMPDIR
