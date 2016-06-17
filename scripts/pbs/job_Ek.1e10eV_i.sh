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
INPUTS="../inputs/TURB_sig1e+0_Nm128_slab.0.2.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.1e10eV_alone_nB50.in"
#NAMES_NODES="-machinefile ./names_1e10eV"
DIR_OUTPUT="../output/output_Ek.1.0e+10eV/"
MON1=$DIR_OUTPUT/info/mon1_iv.log
MON2=$DIR_OUTPUT/info/mon2_iv.log
NPROCS="-np 48"

#mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
mpirun $NPROCS $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#
# NOTES: para metros normales de giacalone99
