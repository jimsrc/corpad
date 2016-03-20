SHELL=/bin/bash
.PHONY: default help object executable all clean

#SOURCE_C = $(wildcard *.cc)
#SOURCE_C = tt.cc nr3_.cc #dummy.cc #nrr3_.cc
SOURCE_C = main_mpi.cc general.cc defs_turb.cc funcs.cc odeintt.cc stepperbs.cc
OBJECTS_C = $(patsubst %.cc, %.o, $(SOURCE_C))
CFLAGS = -O3
CXX = /usr/local/bin/mpic++
EXE = ../CRs.diff.x

${EXE}: $(OBJECTS_C)
	${CXX} -o $@  $^  ${CFLAGS}

%.o: %.cc
	${CXX} ${CFLAGS} -c $^ 

clean:
	rm *.[xo]
