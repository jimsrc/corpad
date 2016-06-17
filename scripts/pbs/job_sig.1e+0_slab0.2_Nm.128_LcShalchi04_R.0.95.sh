#!/bin/bash 
#PBS -l nodes=3:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N R_0.95
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+0_Nm128_slab.0.2_LcShalchi04.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_R.0.95.in"
NAMES_NODES="-machinefile ./nodes_R.0.95"
DIR_OUTPUT="../output/R_0.95"
MON1=$DIR_OUTPUT/info/mon1.log
MON2=$DIR_OUTPUT/info/mon2.log
NPROCS="-np 72"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
