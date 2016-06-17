#!/bin/bash 
#PBS -l nodes=1:ppn=24
#PBS -q larga
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N run.sh
#PBS -j oe
#PBS -k oe

THIS_DIR=/home/jimmy.meza/code/pure/scripts
EXEC=/home/jimmy.meza/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+1_Nm128_slab.0.20.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.1e8eV_alone_nB50.in"
#NAMES_NODES="-machinefile ./names_1e8eV_slab0.20"
DIR_OUTPUT="../output/output_Ek.1.0e+08eV/slab0.20/sig.1.0e+01"
MON1=$DIR_OUTPUT/info/mon1_i.log
MON2=$DIR_OUTPUT/info/mon2_i.log
NPROCS="-np 24"

mpirun $NPROCS $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
