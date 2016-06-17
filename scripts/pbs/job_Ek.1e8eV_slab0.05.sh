#!/bin/bash 
#PBS -l nodes=2:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N 1e8ev_slab0.05
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+0_Nm128_slab.0.05.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.1e8eV_alone_nB50.in"
NAMES_NODES="-machinefile ./names_1e8eV_slab0.05"
DIR_OUTPUT="../output/output_Ek.1.0e+08eV/slab0.05"
MON1=$DIR_OUTPUT/info/mon1_i.log
MON2=$DIR_OUTPUT/info/mon2_i.log
NPROCS="-np 48"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#
# NOTES: parametros normales de giacalone99
