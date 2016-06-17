#!/bin/bash 
#PBS -l nodes=2:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N run.sh
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/INPUT_TURB.Nm128.inp ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/INPUT_GRAL_Ek.1e6_ii.inp"
DIR_OUTPUT="../output/output_Ek.1e6eV_rtol.1e-5"
MON1=$DIR_OUTPUT/info/mon1.log
MON2=$DIR_OUTPUT/info/mon2.log
NPROCS="-np 48"

mpirun $NPROCS $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
