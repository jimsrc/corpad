//#include "funcs1.cpp"
/*#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include "mpi.h"

#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
//#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include "mpi.h"
//#include "my_typedefs.h"
//#include "nr3.h"

using namespace std;

int main(int argc, char* argv[]){
	MPI::Status  status;
    int rank, np;
    int tag=44;
    int flag;
    int *snd, *rcv;
    snd = new int;

	/* Inicio */
	MPI::Init(argc, argv);
	MPI::Intracomm * const comm = &(MPI::COMM_WORLD);

	rank 	= comm->Get_rank();
	np	    = comm->Get_size();


    int src=0, dst=1;
    printf("hhhhhhhhhhhhhhhhhhh \n");

    if( rank==dst ){
        MPI::Request req1;
        flag = comm->Iprobe(src, tag);
        req1 = comm->Irecv(rcv, 1, MPI::INT, src, tag);
        printf(" [r:%d] recibi '%d' de r:%d\n", rank, *rcv, src);
    }
	comm->Barrier(); usleep(10);

    if( rank==src ){
        MPI::Request req0;
        *snd = 22;

        printf(" ----> %d\n", *snd);
        req0 = comm->Isend(snd, 1, MPI::INT, dst, tag);
        req0.Wait();
        printf(" [r:%d] envie '%d' a r:%d\n", rank, *snd, dst);
    }
    comm->Barrier(); usleep(10);

    if( rank==dst ){   
        printf(" [r:%d] tengo '%d' de r:%d\n", rank, *rcv, src);
    }
    comm->Barrier();

	usleep(500);
    printf(" ----> ok.\n");
	//--------------------------
	MPI::Finalize();
	return EXIT_SUCCESS;
}

