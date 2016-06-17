#!/bin/bash 
#PBS -l nodes=3:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N run.sh
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INP_TURB="../inputs/INPUT_TURB.Nm128.inp"
#INP_ORI="../inputs/orientations_isotropic_Nth16_Nph8.in"
INP_ORI="../inputs/orientations_isotropic.inp"
INP_GRAL="../inputs/INPUT_GRAL_Ek.1e8_atol.1e-6.inp"
INPUTS="${INP_TURB} ${INP_ORI} ${INP_GRAL}"
DIR_OUTPUT="../output/output_Ek.1e8eV_atol.1e-6_nB.50"
MON1=$DIR_OUTPUT/info/mon1.log
MON2=$DIR_OUTPUT/info/mon2.log
NPROCS="-np 72"

# backupeamos los inputs
cp -p ${INP_TURB} "${DIR_OUTPUT}/info/turb.in"
cp -p ${INP_ORI} "${DIR_OUTPUT}/info/orientations.in"
cp -p ${INP_GRAL} "${DIR_OUTPUT}/info/plas.in"

# corrida
mpirun $NPROCS $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#EOF
