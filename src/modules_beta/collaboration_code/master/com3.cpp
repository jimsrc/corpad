#include "funcs1.cpp"

//using namespace std;

int main(int argc, char* argv[]){
	int         rank;       // Rank del proceso
	int         p;             // Numero de procesos
	int         source, src;        // Rank del que envia
	int         dst;          // Rank del que recibe
	int         tag_stat=0, tag_coord=1;       // Tag del mensaje
	int 	msg_out, msg_rcv, nx, ny, njobs_per_proc;
	MPI::Status  status;
	MPI::Request req, req0, req1;
	//char        hostname[10];
	//gethostname(hostname, 10);


	/* Inicio */
	MPI::Init(argc, argv);
	//	MPI::Intracomm comm;
	MPI::Intracomm * const comm = &(MPI::COMM_WORLD);
	//comm = &(MPI::COMM_WORLD);

	rank 	= comm->Get_rank();
	p	    = comm->Get_size();
	nx	    = 10;  // nro total de plas
	ny	    = 8;  // nro total de Bs (uso este para dsps equipartir entre los procesadores)
	njobs_per_proc = nx*ny;
	//jobs.build(nx, ny, p, rank);
	//jobs.set_job_flags();
	printf(" ******\n");

    //printf(" rankkkkk: %d \n", rank);
    job jobs;
	msg_handler msgr(nx, ny, &jobs, comm);
    //printf("#@#rankkkkk: %d \n", rank);
	
	int ok, q, flag, aux, rr;
	//double tjob[] = {4, 6, 9, 18}; // [usec]
	double tjob[] = {18, 8, 7, 5}; // [usec]
    for(int i=0; i<4; i++) 
        tjob[i]*=1.0e6; //[usec]
	int *coords;
	comm->Barrier();
	usleep(500);
    try{   
        while(msgr.next()){
            // informo cuanto jobs voy procesando
            //msgr.bcast_my_status();	// envio mi nro total de jobs procesados y coords
            //printf(" [r:%d] *&&nnow: %d\n", rank, jobs.n_now[rank]);
            //const int nnow_ = jobs.n_now[rank];
            //-------- (... work ...)
            usleep(tjob[rank]); //printf(" sleeping... \n");
            //-------- (... work ...)
            msgr.help_routine();
        }
    }
    catch(NRerror s) { 
        printf(" [r:%d] @main() ---> LINE: %d \n", rank, __LINE__);
        NRcatch(s); 
        exit(1);
    }
	//printf(" [r:%d] >>> njobs: %d\n", rank, jobs.n_now[rank]);
	comm->Barrier();
	//msgr.print_nh();
	usleep(500);
    for(int i=0; i<p; i++){
        comm->Barrier();
        if(i==rank){
            printf("\n [r:%d] >>> njobs: %d", rank, jobs.n_now[rank]);
            jobs.print_tsk();
            usleep(500);
        }
    }
	//--------------------------
	MPI::Finalize();
	return EXIT_SUCCESS;
}

