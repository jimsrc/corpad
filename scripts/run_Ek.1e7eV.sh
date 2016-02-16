EXEC=../main_mpi.x
NP=4
INPUTS="../inputs/INPUT_TURB.giacalone.inp ../inputs/init_orientations_15por15.inp ../inputs/INPUT_PLA.Ek_1e7eV.inp"
OUTPUT_DIR="../output_mpi/output_Ek.1e7eV"
mpirun -np $NP $EXEC $INPUTS $OUTPUT_DIR 1> mon1_Ek.1e7eV.log 2> mon2_Ek.1e7eV.log 
