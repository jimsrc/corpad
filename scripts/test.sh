#-----------------------------------------------
EXEC=../CRs.diff.x 
MPIRUN=/usr/local/bin/mpirun
#EXEC=../beta.x
INPUTS="../inputs/INPUT_TURB.ALONE.inp ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/INPUT_GRAL.xx1.inp"
OUTPUT_DIR="../output/xx"
NPROCS="-np 4"
${MPIRUN} $NPROCS $EXEC $INPUTS $OUTPUT_DIR #1> mon1.log 2> mon2.log 
