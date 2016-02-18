#include <mpi.h>	// (**1) 

#include <iostream>
#include <math.h>
#include <cstdlib>
#include "general.cc"
#include "nr3.h"
#include "stepper.h"
#include "defs_turb.cc"
#include "funcs.cc"	//--
#include "odeintt.h"	// "odeint.h"
#include "stepperbs.h"
//#include <mpi.h>	// (**2)
//************************ COMENTARIOS *****************************
// (**1) y (**2): si uso las librerias de anaconda, "mpi.h" debe ir abajo. Si uso las librerias 
// standard de linux, debo ponerlo arriba
//******************************************************************

// main program largely copied from Numerical Recipes Ch 17
int main(int argc, char* argv[]){
	int ord, NPOINTS, n_ori, Nplas_rank, w_rank, w_size, imin, imax, n_Brealiz, nHistTau;
	const Int nvar=6;		// nmbr of y-variables
	double RIGIDITY, FRAC_GYROPERIOD, NMAX_GYROPERIODS, ATOL, RTOL, tmaxHistTau;
	double **array_ori;
	char *fname_turb, *fname_orientations, *fname_gral, *dir_out;
	//char fname_out[100];
	Doub atol, rtol;		// absolute and relative tolerance
	VecDoub ystart(nvar);		// allocate initial x, y[0, 1] values
	bool quest_all, quest_last, cond, exist_file;
	string str_timescale;
	//--------- declaro clases

	fname_turb        	= argv[1];	// parametros de turbulencia SW: solar wind
	fname_orientations	= argv[2];	// orientaciones de las plas
	fname_gral		= argv[3];	// input de propiedades de pla y otros
	dir_out			= argv[4];	// directorio donde guardo las salidas

	//----------------------- leo parametros grales
	read_params(fname_gral, RIGIDITY, FRAC_GYROPERIOD, NMAX_GYROPERIODS, 
			NPOINTS, ATOL, RTOL, n_Brealiz, str_timescale, tmaxHistTau, nHistTau);
	printf(" -------------------->>>>\n");
	// construyo escalas fisicas
	scl.build(RIGIDITY);				// ya declare la variable 'scl' en "./general.cc"
	cout << " scl.beta: "<<scl.beta << endl;
	cout << " rigid: "<< RIGIDITY << endl;
	cout << " frac_gy: "<< FRAC_GYROPERIOD << endl;
	cout << " npoints: " << NPOINTS << endl;
	cout << " ATOL: "<< ATOL << endl;
	cout << " RTOL: " << RTOL << endl;

	//---------------------- inicializa modelos de turbulencia!
	PARAMS par(fname_turb); 			// inputs para c/ region

	Doub h1=FRAC_GYROPERIOD, hmin=0.0;		// init stepsize and minimum-value
	Doub x1=0.0, x2=NMAX_GYROPERIODS;		// init and final x-points

	// settings para las condic iniciales de las plas
	array_ori = read_orientations(fname_orientations, n_ori);	// lista de orientaciones de la veloc inicial
	//xini	  = -2.*par.psw.p_turb.lambda_min / scl.rl;		// [1] posic inic x para todas las plas

	//------------------ output objects w/ 'NPOINTS' points in output
	Output<StepperBS<rhs> > outbs;
	outbs.set_Bmodel(par);
	rhs d;						// object representing system-of-equations
	//---------------------------------------------------------------
	atol	= ATOL;
	if(RTOL==-1){
		rtol = atol;				// as suggested in p.914
	}
	else if(RTOL>=.0){
		rtol = RTOL;
	}
	cout << " ABSOLUTE TOLERANCE = " << atol << endl;
	cout << " RELATIVE TOLERANCE = " << rtol << endl;

	// -------------------------------------- iniciamos Open MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &w_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &w_size);

	Nplas_rank      = n_ori / w_size;
	printf(" Nplas_rank (plas para c/proc): %d/%d", Nplas_rank, n_ori);

    bool helper=false; // no ayudo por defecto
    int i, j; i = j = 0;
    while(j<n_Brealiz){
		// nueva realizacion de Bfield
		par.fix_B_realization(j);
		printf("\n [rank:%d] **************** NEXT Bfield realization j:%d/%d ****************", w_rank, j, n_Brealiz);
        i = 0;  // reinicio el conteo de plas
		//for(int i=0; i<n_ori; i++){
        while(i<n_ori){
			outbs.build(str_timescale, NPOINTS, tmaxHistTau, nHistTau, i, j, dir_out);	// reseteo objeto "Output" para c/pla

			exist_file	= outbs.file_exist();		// le pregunto con "outbs" xq el es quien maneja los nombres finales de los archivos de salida
			imin		= w_rank*Nplas_rank;		// i-minimo para c/procesador
			imax		= (w_rank+1)*Nplas_rank;	// i-maximo para c/procesador
			quest_all	= i>=imin && i<imax;
			quest_last	= i == (w_rank + w_size*Nplas_rank);
			cond		= (quest_all || quest_last || helper) && (!exist_file);
			if(cond){
                outbs.claim_own();  // aviso q YO estyo trabajando con esta pla
				//----------------------------------------- Bulirsch-Stoer
				printf(" [rank:%d] i:%d\n", w_rank, i);
				// nueva pla
				printf(" [rank:%d] corriendo: %s\n", w_rank, outbs.fname_out);
				init_orientation(i, array_ori, ystart);		// values for initial y-values

				Odeint<StepperBS<rhs> > bsode(ystart,x1,x2,atol,rtol,h1,hmin,outbs,d,par,w_rank);	// inicializo el integrador para c/pla
				outbs.tic();				// medimos el tiempo de simulac de c/pla
				bsode.integrate();
				outbs.toc();

				// output --> file
				printf(" [rank:%d] ----> escribiendo: %s\n", w_rank, outbs.fname_out);
				outbs.save2file();								// guardo en archivo
				printf(" [rank:%d] nok/nbad: %d / %d\n", w_rank, bsode.nok, bsode.nbad);
			}
    
            // solo entra si no soy ayudador (helper=false)
            if((i==imax-1) & (j==n_Brealiz-1) & (!helper)){ // si fue mi ultimo trabajo asignado
                helper = true;  // me convierto en un ayudador
                i = -1;      // reseteo conteo de simulaciones
                j = 0;
                par.fix_B_realization(j);
            }
            i++;
		}
        j++;
	}
	//------------------------------------------------
	printf(" [rank:%d] scl.rl [AU]: %g\n", w_rank, scl.rl/AU_in_cm);
	printf(" [rank:%d] scl.wc [s-1]: %g\n", w_rank, scl.wc);
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return EXIT_SUCCESS;
}

/*
 * TODO:
 * - for() ---> while()
 * - generar files "owned.B13.pla003_r.007_h.0"
 * - flag "helper=0,1"
 * - Ahora, file_exist() debe chekear el archivo "owned.."
 *
 * +++++ despues de esta version oficial-temporal +++++
 * - eliminar los srand(), rand()
 *
 */
