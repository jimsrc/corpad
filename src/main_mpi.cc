#include <mpi.h>    // (**1) 

//#include <iostream>
//#include <math.h>
//#include <cstdlib>
//#include "general.cc"
//#include "nr3.h"

//#include "stepper.h"
//#include "defs_turb.cc"
//#include "funcs.cc"   //--
#include "control.h"
#include "general.h"
#include "funcs.h"
#include "odeintt.h"    // "odeint.h"
#include "stepperbs.h"
//#include <mpi.h>  // (**2)
//************************ COMENTARIOS *****************************
// (**1) y (**2): si uso las librerias de anaconda, "mpi.h" debe ir abajo. Si uso las librerias 
// standard de linux, debo ponerlo arriba
//******************************************************************

#ifdef KILL_HANDLER
// I use the static "self-pointer" to access the
// members of the instance 'bsode' inside the main()  =D =D
template<class Stepper> 
Odeint<Stepper>* Odeint<Stepper>::_thisptr = NULL;


// to be called when ctrl-c (SIGINT) signal is sent to process
void signal_handler(int signum){
    //printf(" Caught signal %d\n",signum);
    // Cleanup and close up stuff here
    Odeint<StepperBS<rhs> >::_thisptr->abort_mission(signum);

    // Terminate program
    exit(signum);
}
#endif //KILL_HANDLER


int main(int argc, char* argv[]){
    Int ord, NPOINTS, n_ori, Nplas_rank, w_rank, w_size, imin, imax, n_Brealiz, nHistTau, nThColl;
    const Int nvar=6;       // nmbr of y-variables
    Doub FRAC_GYROPERIOD, NMAX_GYROPERIODS, ATOL, RTOL, tmaxHistTau;
    Doub **array_ori;
    char *fname_turb, *fname_orientations, *fname_gral, *dir_out;
    Doub atol, rtol;        // absolute and relative tolerance
    VecDoub ystart(nvar);       // allocate initial x, y[0, 1] values
    bool quest_all, quest_last, cond, exist_file;
    string str_timescale;
    //--- inputs args
    fname_turb         = argv[1]; // parametros de turbulencia SW: solar wind
    fname_orientations = argv[2]; // orientaciones de las plas
    fname_gral         = argv[3]; // input de propiedades de pla y otros
    dir_out            = argv[4]; // directorio donde guardo las salidas

    //--- leo parametros grales
    read_params(fname_gral, FRAC_GYROPERIOD, NMAX_GYROPERIODS, 
        NPOINTS, ATOL, RTOL, n_Brealiz, str_timescale, tmaxHistTau,
        nHistTau, nThColl);
    printf(" -------------------->>>>\n");
    //--- construyo escalas fisicas
    cout << " frac_gy: "<< FRAC_GYROPERIOD << endl;
    cout << " npoints: " << NPOINTS << endl;
    cout << " ATOL: "<< ATOL << endl;
    cout << " RTOL: " << RTOL << endl;

    //--------------- inicializa modelos de turbulencia!
    PARAMS par(fname_turb);            // inputs para c/ region

    Doub h1=FRAC_GYROPERIOD, hmin=0.0; // init stepsize and minimum-value
    Doub x1=0.0, x2=NMAX_GYROPERIODS;  // init and final x-points

    // settings para las condic iniciales de las plas
    array_ori = read_orientations(fname_orientations, n_ori);   // lista de orientaciones de la veloc inicial

    //------ output objects w/ 'NPOINTS' points in output
    Output<StepperBS<rhs> > outbs;
    outbs.set_Bmodel(&par); 
    rhs d;          // object representing system-of-equations
    //---------------------------------------------------
    atol    = ATOL;
    if(RTOL==-1){
        rtol = atol;                // as suggested in p.914
    }
    else if(RTOL>0.0){  // scrictly grater than zero!
        rtol = RTOL;
    }
    cout << " ABSOLUTE TOLERANCE = " << atol << endl;
    cout << " RELATIVE TOLERANCE = " << rtol << endl;
    // ------------------------ initiate Open MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &w_size);

    Nplas_rank      = n_ori / w_size;
    printf(" Nplas_rank (plas para c/proc): %d/%d", Nplas_rank, n_ori);
    
    bool helper=false; // i don't help by default
    int i=0, j=0;
    
    while(j<n_Brealiz){
        #ifdef KILL_HANDLER
        signal(SIGTERM, signal_handler);
        signal(SIGINT, signal_handler);
        #endif //KILL_HANDLER

        // new B-field realization
        par.fix_B_realization(j);
        printf("\n [rank:%d] ****** NEXT Bfield realization j:%d/%d ******", w_rank, j, n_Brealiz-1);
        i = 0;  // restart ple counting
        while(i<n_ori){
            outbs.build(str_timescale, NPOINTS, i, j, dir_out); // reset "Output" for each ple
            exist_file  = outbs.file_exist(); // i ask 'outbs' because he's who handles output-filenames.
            imin        = w_rank*Nplas_rank;     // i-minimo para c/procesador
            imax        = (w_rank+1)*Nplas_rank; // i-maximo para c/procesador
            quest_all   = i>=imin && i<imax;
            quest_last  = i == (w_rank + w_size*Nplas_rank);
            cond        = (quest_all || quest_last || helper) && (!exist_file);
            if(cond){
                outbs.claim_own();  // tell everyone I'm working w/ this ple
                //----------------------------------- Bulirsch-Stoer
                printf(" [rank:%d] i:%d\n", w_rank, i);
                // nueva pla
                printf(" [rank:%d] corriendo: %s\n", w_rank, outbs.fname_out);
                init_orientation(i, array_ori, ystart); // values for initial y-values

                Odeint<StepperBS<rhs> > bsode(ystart,x1,x2,atol,rtol,
                    h1,hmin,outbs,d,par,nHistTau,nThColl,w_rank); // initialize the solver for each ple

                #ifdef KILL_HANDLER
                bsode._thisptr = &bsode; // 'bsode' address 
                #endif //KILL_HANDLER

                outbs.tic();                // medimos el tiempo de simulac de c/pla
                bsode.integrate();
                outbs.toc();

                // output --> file
                printf(" [r:%d] ---> writing: %s\n",w_rank,outbs.fname_out);
                outbs.save2file();                              // guardo en archivo
                printf(" [r:%d] nok/nbad: %d / %d\n",w_rank,bsode.nok,bsode.nbad);
            }
            // solo entra si no soy ayudador (helper=false)
            if((i==imax-1) & (j==n_Brealiz-1) & (!helper)){ // si fue mi ultimo trabajo asignado
                helper = true;  // i become a helper
                i = -1;      // reset particle counter
                j = 0;       // reset B-field counter
                par.fix_B_realization(j);
            }
            i++;
        }
        j++;
    }
    
    //------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

/*
 * TODO:
 * [x] for() ---> while()
 * [x] generar files "owned.B13.pla003_r.007_h.0"
 * [x] flag "helper=0,1"
 * [x] Ahora, file_exist() debe chekear el archivo "owned.."
 * [ ] rann0() ----> ZIGGURATZ
 * [ ] se pueden evitar las declaraciones en void MODEL_TURB::calc_dB_SLAB(..)??
 * [ ] hacer q las variables tipo VeDoub (e.g. ysave, ym, yn, etc) sean 
 *     miembros privados de la clase, para asi evitar q se 
 *     inicialicen c/vez q use StepperBS::step(), StepperBS::dy(), 
       Odeint::Odeint(.., Output out,..), etc.
 * [ ] en las mismas rutinas del punto anterior, convertir los argumentos
 *     tipo &ptr --> *ptr. Es decir q el argumento NO sea tipo *referencia*,
 *     sino tipo *puntero*!! (las referencias generan copias!)
 *
 * +++++ despues de esta version oficial-temporal +++++
 * - eliminar los srand(), rand()
 *
 */
