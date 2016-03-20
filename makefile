# el directorio oficial esta en "master"
DIR_SRC=./src
DIR_BETA=./src/beta
EXE_WMPI=./CRs.diff.x 
EXE_WOMPI=./CRs.diff_woMPI.x 
EXE_WPROF=./CRs.diff_wgprof.x 	# profiling / debugging
MPICXX=/usr/local/bin/mpic++
OPTIM=-O3 #-O3
OPT_GP=-O0 #-O3 #-O0 	# for use w/ gprof

.PHONY: default help object executable all clean
SOURCE_C = $(wildcard *.cc)
OBJECTS_C = $(patsubst %.cc, %.o, $(SOURCE_C))


default:
	make w_mpi

%.o: %.cc
	g++ -c $^ -o $@



beta:
	${MPICXX} ${DIR_BETA}/main_mpi.cc -o ./beta.x


w_mpi:
	${MPICXX} ${DIR_SRC}/main_mpi.cc ${OPTIM} -o ${EXE_WMPI} && ls -lhtr ${EXE_WMPI}


mpi-w-gprof:
	${MPICXX} -pg ${OPT_GP} ${DIR_SRC}/main_mpi.cc -o ${EXE_WPROF}


wo_mpi:
	g++ ${DIR_SRC}/main.cc -o ${EXE_WOMPI}


clean:
	rm *.x


#++++++++++++++++++++++++++++++++++++++++

