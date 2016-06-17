#!/bin/bash 
#PBS -l nodes=2:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N 1e6eV/s.1e+2/Nm256
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+2_Nm256.in ../inputs/orientations_isotropic.inp ../inputs/GRAL_Ek.1e6eV_alone.in"
NAMES_NODES="-machinefile ./nodes_sig.1e+2_1e6eV_Nm256"
DIR_OUTPUT="../output/output_Ek.1e6eV/sig.1e+2/Nm256"
MON1=$THIS_DIR/mon1_sig.1e+2_1e6eV_Nm256.log
MON2=$THIS_DIR/mon2_sig.1e+2_1e6eV_Nm256.log
NPROCS="-np 48"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
