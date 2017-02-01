#!/bin/csh
###############################################################################
# SCRIPT FOR SUBMITTING main.py AS A SERIAL JOB TO A QUEUING SYSTEM
###############################################################################
#
# Description
#    Run this script from the tool main directory as
#      csh util/jobs/main.csh nml/<namelist>
#
# History
#    20140611-A_RiMa: written.
#
###############################################################################

if ($1 == "") then
    echo "Usage csh util/jobs/main.csh nml/<namelist>"
    exit 1
endif

set outlog = `basename $1`.out
set errlog = `basename $1`.err

switch ($HOST)
    # DLR PA1 CLUSTER
    case 'pa1':
        set currpath = `pwd`
        cat <<EOF | qsub -
#!/bin/sh
#PBS -V
#PBS -S /bin/sh
#PBS -N ${1}
#PBS -o ${currpath}/$outlog
#PBS -e ${currpath}/$errlog
#PBS -l nodes=1:ppn=1,walltime=24:00:00
module load ncl/6.2.0-nodap
cd ${currpath}
python main.py ${1}
EOF
        breaksw

    # DEFAULT
    default:
        echo "Unknown host $HOST"
        exit 1
        breaksw
endsw

exit 0
