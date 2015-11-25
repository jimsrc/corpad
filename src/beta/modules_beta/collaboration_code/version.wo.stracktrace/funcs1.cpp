#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include "mpi.h"
//#include "my_typedefs.h"
#include <iostream>
#include <stacktrace/call_stack_gcc.cpp>
#include <stacktrace/stack_exception.hpp>
#include "nr3.h"
//#include "backward.hpp"

//typedef int dtp;    // datatype (DTp)


// se usa: LiberaMat<double>(mat, 4) por ejemplo.
template <typename dtp>
void LiberaMat(dtp **M, int iFilas){
    int j;
    for(j=0; j<iFilas; j++){
        //if((Mat[i] = (double *) malloc(nCol* sizeof(double)))==NULL) {^M
        free(M[j]);//&M+i?
    }   
free(M);
}


// se usa a=Allocmat<double>(3,4) por ejemplo.
template <typename dtp>
dtp **Allocmat(int nFilas, int nCol){
    dtp **Mat;
    int i;
    if((Mat = (dtp**) malloc(nFilas * sizeof(dtp *))) == NULL)
        return NULL;//dice que, desde ya, no se pudo alojar memoria pa la matriz
    for(i=0; i<nFilas; i++)
        if((Mat[i] = (dtp *) malloc(nCol* sizeof(dtp)))==NULL) {
            LiberaMat(Mat,i);//libera la matriz (aborta mision).
            return NULL;
        }   
        return Mat;
}


class job{
	public:
		job(void){};
		void build(int, int, int, int);
		int nx, ny;
		MatInt tsk; //int **tsk;	// array de tasks para c/procesador
		int p, my_rank;
		void set_job_flags();
		int next();
		void print_tsk();
		int *i_now, *j_now, *n_now;
		void check_balance(void);
		int rank_helped;
		bool should_i_help();
		int who_needs_help();
		int *coords_helped;
	private:
		int nmin, nmax;
};

void job::build(int nxx, int nyy, int pp, int myr){
	nx = nxx;
	ny = nyy;
	p  = pp;
	my_rank = myr;
	//tsk  = (int**) calloc(nx, sizeof(int*));
    tsk = MatInt(nx,ny);
	for(int i=0; i<nx; i++){
		//tsk[i] = (int*) calloc(ny, sizeof(int));
		for(int j=0; j<ny; j++)
			tsk[i][j] = -1;
	}
	n_now = new int[p];	// nro de tareas en progress
	i_now = new int[p];
	j_now = new int[p];
	coords_helped = new int[2];	// coords a las q ayudo
	for(int i=0; i<pp; i++) 
		n_now[i] = 0;	// tareas que marco como "pendientes"
}

bool job::should_i_help(){
	// hallo los nmin/nmax globales:
	// NOTA: mas de uno puede tener el mismo avance
	// que el minimo global 'nmin'. Idem para 'nmax'
	nmin = 9999999;
	nmax = 0;
	for(int r=0; r<p; r++){
		nmin = std::min(nmin, n_now[r]);
		nmax = std::max(nmax, n_now[r]);
	}
	// determino si alguien necesita ayuda:
	if(nmax-nmin>=2){			// alguien esta avanzando mas q otro?
		printf("\033[31m [r:%d] ===UNBALANCED===\033[0m \n", my_rank);

		if(nmax==n_now[my_rank]){	// me corresponde a mi ayudar?? (*)
			rank_helped 	 = who_needs_help(); // a quien necesito ayudar??
            //posiciones candidato a ayudar 
            if(i_now[rank_helped]+1 < nx-1){
                coords_helped[0] = i_now[rank_helped]+1; // (**)
                coords_helped[1] = j_now[rank_helped];
            }
            else if(j_now[rank_helped]+1 < ny-1){
                coords_helped[0] = i_now[rank_helped]; //
                coords_helped[1] = j_now[rank_helped]+1; // (**)
            }
            else{ // no me imagino en q caso podria entrar aqui
                printf("\n -----> NO SABEMOS EN Q POSICION AYUDAR!!\n");
                exit(1);
            }
			return true;
		}
		else				// entonces le corresponde a otro
			return false;
	}
	else					// todos los jobs estan avanzando "igual"
		return false;
	// (*) mas de uno puede tener el mismo avance ('n_now[my_rank]') que
	// el maximo global ('nmax')
    // (**) notar q 'i_now[rank_helped] + 1' nunca puede dar una posicion ilegal, ya que la condicion "nmax-nmin>=2" (recordar q 'rank_helped' es igual a 'nmin' por definicion) garantiza q el procesador de rango 'rank_helped' le falta al menos 2 posiciones para terminar todo!
}

int job::who_needs_help(){
	// busco el rango 'r' q necesita ayuda
	for(int r=0; r<p; r++){	
		if( n_now[r] == nmin )
			return r;	// el 1ero q tenga el mismo avance que
					// el minimo global (*)
	}
	// (*): notar q 'my_rank' nunca puede ser 'r', xq en
	// "should_i_help()" ya nos aseguramos q nmax!=nmin y 
	// que n_now[my_rank]==nmax.
}

void job::check_balance(){
	nmin = 9999999;
	nmax = 0;
	for(int i=0; i<p; i++){
		nmin = std::min(nmin, n_now[i]);
		nmax = std::max(nmax, n_now[i]);
	}
	if(nmax - nmin>=2)
		printf(" [r:%d] ===UNBALANCED===\n", my_rank);
}

void job::set_job_flags(){
	int np = ny/p;
	int nyy = np*(my_rank+1);
	int jj  = np*(my_rank);
	for(int i=0; i<nx; i++){
		for(int j=jj; j<nyy; j++){
			tsk[i][j] = 0;
		}
	}
	print_tsk();	// imprime flags en alas coordenadas de los jobs
}

void job::print_tsk(){
	int np = ny/p;
	int nyy = np*(my_rank+1);
	int jj  = np*(my_rank);

	printf("\n tsk(r:%d):\n", my_rank);
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
            //if(j==my_rank)
            if(j>=jj && j<nyy)
                printf("\033[33m %d\033[0m", tsk[i][j]);       
            else
                printf(" %d", tsk[i][j]);       
        }
		printf("\n");
	}
}

// 0: si ya termine todos mis jobs
// 1: si aun faltan jobs
int job::next(){
	//calc_progress();
    for(int j=0; j<ny; j++) // busco si hay mas B-realizations pendientes
        for(int i=0; i<nx; i++) // primero busca si hay mas plas pendientes
			if(tsk[i][j]==0){
				//printf(" -----> buscando..\n");
				i_now[my_rank] = i;
				j_now[my_rank] = j;
				tsk[i][j] = 1;
				n_now[my_rank]++;//nro de jobs hechos (incluyendo el q esta en progreso)
			    printf(" [r:%d] nnow++: %d\n", my_rank, n_now[my_rank]);
				return 1;
			}
	printf("\033[33m [r:%d] termine todo! :) \033[0m\n", my_rank);
	return 0;
}

/*********************************************************************************
**********************************************************************************/

class msg_handler{
	public:
		msg_handler(int, int, job *, MPI::Intracomm *);
		void bcast_my_status(void);
		job* jobs;
		int p;
        int rank;
		int tag_stat;
        int tag_coord;
		int *qsnd;
        MatInt qrcv;
		void recv_status(void);
		void bcast_coords_for_helped_guy(void);
		int *coords;
		int dst;
		MPI::Intracomm* comm;
		void gather_coords_for_helped_guy(void);
		void update_my_progress_after_help(int*, int);
        bool should_i_help(void);
        int next(void);
};


msg_handler::msg_handler(int nx, int ny, job * const jobss, MPI::Intracomm * const commm){
	comm = commm;
	p    = comm->Get_size();//pp;
    rank = comm->Get_rank();
    //printf(" ---> @msg_handler: rank:%d \n", rank);
	//rank = rankk;

	jobs = jobss;
    //jobs = (job*) malloc(sizeof(job*));
    jobs->build(nx, ny, p, rank);
    jobs->set_job_flags();

    // banderas de comunicacion MPI
	tag_stat  = 0;
	tag_coord = 1;

	qsnd = new int[3];

    qrcv = MatInt(p,3);//Allocmat<int>(p,3);//aloco matriz de enteros
    for(int i=0; i<p; i++)
        for(int j=0; j<3; j++)
            qrcv[i][j] = -1;
    //printf(" --->$@msg_handler: rank:%d \n", rank);
}


bool msg_handler::should_i_help(){
    return jobs->should_i_help();
}


int msg_handler::next(){
    return jobs->next();
}

/*
* Esto es un msje dirigido al que queremos ayudar (cuyo 
* rango es 'jobs->rank_helped'). Le mandamos las coordenadas en 
* que queremos ayudarle para que el pueda actualizar su "tabla"
* de tareas.
*/
void msg_handler::bcast_coords_for_helped_guy(void){//MPI::Intracomm* comm){
	MPI::Request req;
    int *crds;// = new int[2];
	dst  		= jobs->rank_helped;	// [int] rank del q necesita ayuda
	crds 		= jobs->coords_helped;	// [pointer] coords del "job-helped"
	req = comm->Isend(crds, 2, MPI::INT, dst, tag_coord);
	req.Wait();
    jobs->tsk[crds[0]][crds[1]] = -2; // flag de "yo ayude aqui"
}


/* 
 * Esto es una rutina para q el "ayudado" reciba las coordenadas de 
 * la posicion en donde recibe ayuda.
 * Si no hay nada q recibir, o nadie esta ayudando realmente, esta 
 * rutina no hace absolutamente NADA
 */
void msg_handler::gather_coords_for_helped_guy(void){
	int flag;
	//int * const crds = new int[2]; // puntero con direcc constante
	int * crds = new int[2]; // puntero 
    //int crds[2];
	MPI::Request req;
    //const int nnow_ = jobs->n_now[rank];
    //if(nnow_ > jobs->n_now[rank]) throw(" ---> MMMMMMMMM!\n");
	for(int src_=0; src_<p; src_++) if(src_!=rank){
        const int src = src_;
		crds[0] = -1;	// para luego chekar si recibi algo
		flag = comm->Iprobe(src, tag_coord);
		req  = comm->Irecv(crds, 2, MPI::INT, src, tag_coord); usleep(1e6);
		//printf(" ---> crds: %d\n", crds[0]);
		if(crds[0]!=-1){// si recibo algo, alguien me quiere ayudar, por 
				        // lo tanto crds[0] debe ser !=-1
            printf("\033[1;34m [r:%d] esta intentando ayudarme r:%d @(%d,%d)  \033[0m\n", rank, src, crds[0], crds[1]);
			update_my_progress_after_help(crds, src);   // sumo 1 para contar la ayuda q recibo,
					//marco flag "1" en tsk[][], y hago un bcast() a todos 
					//para q actualizarles de mi avance; asi evito q alguien 
					//mas me ayude en las mismas coordenadas! :)
		}
	}
    //if(nnow_ > jobs->n_now[rank]){
    //    printf(" [r:%d] bef:%d, aft:%d\n", rank, nnow_, jobs->n_now[rank]);
    //    throw("aaaaaaaa");
    //}
}


// esto se ejecuta si recibi ayuda: solo me importa actualizar
// mi estado de avance y avisar al resto de ello. No importa *quien*
// me ayuda.
void msg_handler::update_my_progress_after_help(int *crds, int src){
    //+++++++ tal vez estoy recibiendo una ayuda antigua
	if(jobs->tsk[crds[0]][crds[1]]==0){ // si realmente me esta ayudando (*)
		printf("\033[1;34m [r:%d] ME AYUDA r:%d @(%d,%d) \033[0m\n", rank, src, crds[0], crds[1]);
		jobs->n_now[rank]++;
		jobs->tsk[crds[0]][crds[1]] = 2; // flag de "me ayudaron aqui"
		bcast_my_status();	// comunico a todos de mi nuevo avance
	}
    //+++++++ para ver q no este pasando NADA raro
	else if(jobs->tsk[crds[0]][crds[1]]==-1){ // chekeo de consistencia (**)
		printf(" [r:%d] ---> ERROR: [%d] me esta tratando de ayudar donde no debe!", rank, src);
        char *msg; msg = new char[50];
        sprintf(msg, "NRmatrix subscript out of bounds, i/nn: %d/%d", i, nn);
        throw(msg);
        //call_stack st;
        //cout << st.to_string();
	}
    //+++++++ si es q realmente fue una ayuda antigua 
    else{
        printf("\033[1;34m [r:%d] THANKU, BUT NO THANKU!! r:%d @(%d,%d) \033[0m\n", rank, src, crds[0], crds[1]);
    }
	// (*) las viejas ayudas del mismo procesaor 'src' no
	// pasaran de este if()
    // (**) solo chequeo q no me ayude en "-1". Si es otro flag, 
    // no importa. Usualmente siginifica q trata de ayudar de nuevo
    // en las mismas coordenadas (ayudas viejas). 
}


void msg_handler::bcast_my_status(){
	MPI::Request req;
    const int nnow_ = jobs->n_now[rank];
	qsnd[0] = jobs->n_now[rank];
	qsnd[1] = jobs->i_now[rank];
	qsnd[2] = jobs->j_now[rank];
    printf(" [r:%d] ----- begin bcast() \n", rank);
	for(int dst=0; dst<p; dst++) if(dst!=rank){
		printf(" [r:%d] ---> [r:%d]\n", rank, dst);
		req = comm->Isend(qsnd, 3, MPI::INT, dst, tag_stat);
		req.Wait();
	}
    printf(" [r:%d] ----- finish bcast() \n", rank);
    if(nnow_ > jobs->n_now[rank]){
        printf(" [r:%d] bef:%d, aft:%d\n", rank, nnow_, jobs->n_now[rank]);
        throw("aaaaaaaa");
    }
}


void msg_handler::recv_status(){
	int flag, cc; 
	MPI::Request req;
	for(int src_=0; src_<p; src_++) if(src_!=rank){
        const int src = src_; // para evitar raros cambios automaticos (locos! :$ WTF???!!!!)
        //if(src>=p || src<0 || src==rank) printf(" ~~~src/p: %d/%d\n", src, p);
		for(int i=0; i<10; i++){ //para asegurarme de recibir el ultimo send
            //printf(" [r:%d] aqui-1a, src=%d ", rank, src);
			flag = comm->Iprobe(src, tag_stat); //printf(" ok1, ");
			req  = comm->Irecv(qrcv[src], 3, MPI::INT, src, tag_stat); //printf(" ok2.\n");
		}
        //cc = src!=rank?1:0;
        //cc *= src<p?1:0;
        //printf(" ok:%d, ", cc);
		jobs->n_now[src] = qrcv[src][0];
		jobs->i_now[src] = qrcv[src][1];
		jobs->j_now[src] = qrcv[src][2];
		printf(" [r:%d] qrcv(r=%d): %d   @(%d,%d); \n", rank, src,
			qrcv[src][0], qrcv[src][1], qrcv[src][2]);
	}
}

//
