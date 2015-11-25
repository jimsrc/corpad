#-----------------------------------------------
EXEC=../main_mpi.x
INPUTS="../inputs/TURB_sig1.0.in ../inputs/orientations_isotropic.inp ../inputs/GRAL_Ek.1e6eV_alone.in"
OUTPUT_DIR="../output/output_Ek.1e6eV"
NPROCS="-np 2"
MON1="mon1_1e6eV.log"
MON2="mon2_1e6eV.log"
mpirun $NPROCS $EXEC $INPUTS $OUTPUT_DIR 1> $MON1 2> $MON2
