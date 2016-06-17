#!/bin/bash 
#PBS -l nodes=1:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N 1e6eV.Nm2
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/main_mpi.x
cd $THIS_DIR
INPUTS="../inputs/INPUT_TURB.Nm2.inp ../inputs/orientations_isotropic.inp ../inputs/INPUT_GRAL_Ek.1e6.inp"
NAMES_NODES="-machinefile ./names.list_Nm2"
DIR_OUTPUT="../output/output_Ek.1e6eV/Nm2"
MON1=$THIS_DIR/mon1_1e6eV_Nm2.log
MON2=$THIS_DIR/mon2_1e6eV_Nm2.log
NPROCS="-np 24"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
