#!/bin/bash 
#PBS -l nodes=1:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N slab0.01_Nm256
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+0_Nm256_slab.0.01_LcShalchi04.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.0.6e6eV_alone.in"
NAMES_NODES="-machinefile ./nodes_slab0.01_shalchi04"
DIR_OUTPUT="../output/output_Ek.6.0e+05eV/slab0.01/shalchi04/Nm256"
MON1=$DIR_OUTPUT/info/mon1_slab0.01_shalchi04.log
MON2=$DIR_OUTPUT/info/mon2_slab0.01_shalchi04.log
NPROCS="-np 24"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
