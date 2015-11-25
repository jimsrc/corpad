# el directorio oficial esta en "master"
DIR_SRC=./src/master
DIR_BETA=./src/beta
EXE_WMPI=./CRs.diff.x 
EXE_WOMPI=./CRs.diff_woMPI.x 
EXE_WPROF=./CRs.diff_wgprof.x 	# profiling / debugging
MPICXX=/usr/local/bin/mpic++

default:
	make w_mpi

beta:
	${MPICXX} ${DIR_BETA}/main_mpi.cc -o ./beta.x


w_mpi:
	${MPICXX} ${DIR_SRC}/main_mpi.cc -o ${EXE_WMPI}


mpi-w-gprof:
	${MPICXX} -pg ${DIR_SRC}/main_mpi.cc -o ${EXE_WPROF}


wo_mpi:
	g++ ${DIR_SRC}/main.cc -o ${EXE_WOMPI}


clean:
	rm *.x
