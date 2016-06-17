#!/bin/bash 
#PBS -l nodes=5:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N 1e9ev
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1.0_L2d1e-3.Lslab1e-2_Nm128.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.1e9eV_alone_nB50.in"
NAMES_NODES="-machinefile ./nodes_1e9eV"
DIR_OUTPUT="../output/output_Ek.1.0e+09eV/Lc2d1e-3.Lslab1e-2"
MON1=$THIS_DIR/mon1_1e9eV.log
MON2=$THIS_DIR/mon2_1e9eV.log
NPROCS="-np 120"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
