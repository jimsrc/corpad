# el directorio oficial esta en "master"
DIR_SRC=./src
DIR_BETA=./src/beta
EXE_WMPI=./CRs.diff.x 
EXE_WOMPI=./CRs.diff_woMPI.x 
EXE_WPROF=./CRs.diff_wgprof.x 	# profiling / debugging
MPICXX=/usr/local/bin/mpic++
OPTIM=-O3 #-O3
OPT_GP=-O0 #-O3 #-O0 	# for use w/ gprof

default:
	make w_mpi

beta:
	${MPICXX} ${DIR_BETA}/main_mpi.cc -o ./beta.x


w_mpi:
	#${MPICXX} ${DIR_SRC}/main_mpi.cc -O2 -o ${EXE_WMPI}
	${MPICXX} ${DIR_SRC}/main_mpi.cc ${OPTIM} -o ${EXE_WMPI}


mpi-w-gprof:
	${MPICXX} -pg ${OPT_GP} ${DIR_SRC}/main_mpi.cc -o ${EXE_WPROF}


wo_mpi:
	g++ ${DIR_SRC}/main.cc -o ${EXE_WOMPI}


clean:
	rm *.x
