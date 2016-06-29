#!/bin/bash 
#PBS -l nodes=2:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N run.sh
#PBS -j oe
#PBS -k oe
THIS_DIR=$HOME/plas/scripts/pbs
EXEC=$HOME/plas/CRs.diff.x #new code
INPUT=$HOME/plas/inputs
cd $THIS_DIR
INP_TURB="$INPUT/turb.in"
#INP_ORI="../inputs/ori_242.in"
INP_ORI="$INPUT/ori_122.in"
INP_GRAL="$INPUT/plas.in"
INPUTS="${INP_TURB} ${INP_ORI} ${INP_GRAL}"
DIR_OUTPUT="../../output/r.0.50_NmS.128_Nm2d.256"
MON1=$DIR_OUTPUT/info/mon1.log
MON2=$DIR_OUTPUT/info/mon2.log
NPROCS="-np 48"

# backupeamos los inputs
cp -p ${INP_TURB} "${DIR_OUTPUT}/info/turb.in"
cp -p ${INP_ORI} "${DIR_OUTPUT}/info/orientations.in"
cp -p ${INP_GRAL} "${DIR_OUTPUT}/info/plas.in"

# corrida
mpirun $NPROCS $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#EOF
