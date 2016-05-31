SHELL=/bin/bash
.PHONY: default help object executable all clean

#SOURCE_C = $(wildcard *.cc)
SOURCE_C = main_mpi.cc general.cc defs_turb.cc funcs.cc odeintt.cc stepperbs.cc
OBJECTS_C = $(patsubst %.cc, %.o, $(SOURCE_C))
OBJECTS_C_GPROF = $(patsubst %.cc, %_gp.o, $(SOURCE_C))
CFLAGS = -O3 -Wlogical-op -fprofile-use #-fprofile-use #-fprofile-generate
CXX = /usr/local/bin/mpic++
#EXE = ../CRs.diff.x
EXE = ./main.x
# gprof stuff
EXE_WPROF=./main_wgprof.x 	# profiling / debugging
OPT_GP=-O0 #-O3 #-O0 	# for use w/ gprof

all: $(OBJECTS_C)
	${CXX} -o $(EXE)  $^  ${CFLAGS}
	ls -lhtr $(EXE)

%.o: %.cc
	${CXX} ${CFLAGS} -c $^ 

%_gp.o: %.cc
	$(CXX) $(OPT_GP) -c $^ -o $@

gprof_exe: $(OBJECTS_C_GPROF)
	${CXX} -pg ${OPT_GP} $^ -o $(EXE_WPROF)

clean:
	rm *.[xo] #*.gcda

# INFO ON MORE CFLAGS:
#-Wlogical-op
#   Warn about suspicious uses of logical
#   operators in expressions.  This
#   includes using logical operators in
#   contexts where a bit-wise operator is
#   likely to be expected.
#
#-fprofile-generate
#-fprofile-generate=path
#   Enable options usually used for
#   instrumenting application to produce
#   profile useful for later recompilation
#   with profile feedback based
#   optimization.  You must use
#   -fprofile-generate both when compiling
#   and when linking your program.
#
#   The following options are enabled:
#   "-fprofile-arcs", "-fprofile-values",
#   "-fvpt".
#
#   If path is specified, GCC will look at
#   the path to find the profile feedback
#   data files. See -fprofile-dir.
#
#-fprofile-use
#-fprofile-use=path
#   Enable profile feedback directed
#   optimizations, and optimizations
#   generally profitable only with profile
#   feedback available.
#
#   The following options are enabled:
#   "-fbranch-probabilities", "-fvpt",
#   "-funroll-loops", "-fpeel-loops",
#   "-ftracer"
#
#   By default, GCC emits an error message
#   if the feedback profiles do not match
#   the source code.  This error can be
#   turned into a warning by using
#   -Wcoverage-mismatch.  Note this may
#   result in poorly optimized code.
#
#   If path is specified, GCC will look at
#   the path to find the profile feedback
#   data files. See -fprofile-dir.
