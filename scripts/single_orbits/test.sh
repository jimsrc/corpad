#-----------------------------------------------
EXEC=../CRs.diff.x
INPUTS="../inputs/INPUT_TURB.ALONE.inp ../inputs/orientations_isotropic.inp ../inputs/INPUT_GRAL.xx1.inp"
OUTPUT_DIR="../output/xx"
NPROCS="-np 3"
mpirun $NPROCS $EXEC $INPUTS $OUTPUT_DIR #1> mon1.log 2> mon2.log 
