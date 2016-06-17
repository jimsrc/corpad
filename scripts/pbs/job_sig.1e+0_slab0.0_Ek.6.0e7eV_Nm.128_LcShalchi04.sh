#!/bin/bash 
#PBS -l nodes=2:ppn=24
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N slab0.0
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+0_Nm128_slab.0.0_LcShalchi04.in ../inputs/orientations_isotropic.inp ../inputs/GRAL_Ek.6e7eV_alone.in"
NAMES_NODES="-machinefile ./nodes_slab0.0_shalchi06"
DIR_OUTPUT="../output/output_Ek.6.0e+07eV/slab0.0/shalchi06"
MON1=$THIS_DIR/mon1_slab0.0_shalchi06.log
MON2=$THIS_DIR/mon2_slab0.0_shalchi06.log
NPROCS="-np 48"

mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
