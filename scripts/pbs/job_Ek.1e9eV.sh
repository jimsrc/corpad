#!/bin/bash 
#PBS -l nodes=4:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N 1e9eV
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+0_Nm128_slab.0.2.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.1e9eV_alone_nB50.in"
NAMES_NODES="-machinefile ./names_1e9eV"
DIR_OUTPUT="../output/output_Ek.1.0e+09eV/"
MON1=$THIS_DIR/mon1_1e9eV.log
MON2=$THIS_DIR/mon2_1e9eV.log
NPROCS="-np 96"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#
# NOTES: para metros normales de giacalone99
