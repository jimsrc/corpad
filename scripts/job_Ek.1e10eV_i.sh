#PBS -l nodes=1:ppn=8
#PBS -l walltime=38:00:00
#PBS -M jimmy.ilws@gmail.com
#PBS -m abe 
#PBS -N 1e10eV
#PBS -j eo

THIS_DIR=$HOME/code/pure/scripts
EXEC=$HOME/code/pure/CRs.diff.x
cd $THIS_DIR
INPUTS="../inputs/TURB_sig1e+0_Nm128_slab.0.2.in ../inputs/orientations_isotropic_Nth16_Nph8.in ../inputs/GRAL_Ek.1e10eV_alone_nB50.in"
#NAMES_NODES="-machinefile ./names_1e10eV"
DIR_OUTPUT="../out/output_Ek.1.0e+10eV/"
MON1=$DIR_OUTPUT/info/mon1_iv.log
MON2=$DIR_OUTPUT/info/mon2_iv.log
#NPROCS="-np 24"

#+++++++++++++++++++++++++++++++++++++++++++++
NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')
NODES=$(sort $PBS_NODEFILE | uniq | wc -l)

echo "nodes ($NP cpu total):"
sort $PBS_NODEFILE | uniq

cd $PBS_O_WORKDIR
#echo $NP > file_mon
#echo $NODES >> file_mon
#+++++++++++++++++++++++++++++++++++++++++++++
# Ejecuta el trabajo:
mpiexec --comm=pmi -np $NP $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#mpiexec -np $NP $EXE $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2

#mpirun $NPROCS $NAMES_NODES $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
#mpirun $NPROCS $EXEC $INPUTS $DIR_OUTPUT 1> $MON1 2> $MON2
exit 0

