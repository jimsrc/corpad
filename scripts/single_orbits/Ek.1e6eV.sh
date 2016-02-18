#-----------------------------------------------
EXEC=../../CRs.diff.x
INPUTS="../../inputs/single_orbits/TURB.sig1.0e+00_slab0.2_Bo.05.in ../../inputs/orientations_isotropic_Nth16_Nph8.in ../../inputs/single_orbits/GRAL_Ek.1e6eV_alone.in"
OUTPUT_DIR="../../output/single_orbits/Ek.1.0e+06eV/Nm128/slab0.20/sig.1.0e+00/Lc2D.1.0e-02_LcSlab.1.0e-02"
MON1=$DIR_OUTPUT/info/mon1.log
MON2=$DIR_OUTPUT/info/mon2.log
NPROCS="-np 2"
mpirun $NPROCS $EXEC $INPUTS $OUTPUT_DIR #1> $MON1 2> $MON2
