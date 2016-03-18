#!/bin/bash
EXE=../CRs.diff.x #../CRs.diff_wgprof.x  #../CRs.diff.x 
MPIRUN=/usr/local/bin/mpirun
INP_TURB="../inputs/INPUT_TURB.ALONE.inp"
#INP_ORI="../inputs/orientations_isotropic_Nth16_Nph8.in"
INP_ORI="../inputs/orientations_small.in" # 2nd half of ...Nth16_Nph8.in
INP_GRAL="../inputs/INPUT_GRAL.xx1.inp"
INPUTS="$INP_TURB $INP_ORI $INP_GRAL"
OUTPUT_DIR="../output/xx"
NPROCS="-np 4"
${MPIRUN} $NPROCS $EXE $INPUTS $OUTPUT_DIR #1> mon1.log 2> mon2.log 
#EOF
