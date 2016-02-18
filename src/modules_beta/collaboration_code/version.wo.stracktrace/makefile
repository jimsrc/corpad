# Use GNU C++ compiler
CC=g++
STACKTRACE_DIR=../../shared_lib
EXE=./test.x 
SRC=./com3.cpp
# Compiler settings (warning on all)
CFLAGS=-c -Wall -I${STACKTRACE_DIR}


this:
	mpic++ -Wall -I${STACKTRACE_DIR} ${SRC} -ldl -rdynamic -o ${EXE}

all: call_stack_gcc.o test_code.o
	$(CC) ./call_stack_gcc.o ./test_code.o -ldl -rdynamic -o ${EXE}

./call_stack_gcc.o:
	$(CC) $(CFLAGS) ../stacktrace/call_stack_gcc.cpp -o ./call_stack_gcc.o

./test_code.o:
	mpic++ -Wall -I${STACKTRACE_DIR} ${SRC} -o ./test_mpi.x
	mpic++ $(CFLAGS)  ${SRC} -o ./test_code.o

clean:
	rm -rf *.o *.x
