#include <iostream>
#include <math.h>
#include <cstdlib>
#include "general.cc"
#include "nr3.h"
#include "stepper.h"
#include "defs_turb.cc"
#include "funcs.cc"	//--
#include "odeint.h"
#include "stepperbs.h"

// main program largely copied from Numerical Recipes Ch 17
int main(int argc, char* argv[]){
	int ord, NPOINTS, n_ori, n_Brealiz;
	const Int nvar=6;		// nmbr of y-variables
	double RIGIDITY, FRAC_GYROPERIOD, NMAX_GYROPERIODS, ATOL, RTOL;
	double **array_ori;
	char *fname_input_turb, *fname_orientations, *fname_input_pla;
	char dir_out[100], fname_out[100];
	Doub atol, rtol;		// absolute and relative tolerance
	VecDoub ystart(nvar);		// allocate initial x, y[0, 1] values
	//--------- declaro clases

	fname_input_turb        = argv[1];
	fname_orientations	= argv[2];
	fname_input_pla		= argv[3];

	read_params(fname_input_pla, RIGIDITY, FRAC_GYROPERIOD, NMAX_GYROPERIODS, NPOINTS, ATOL, RTOL, n_Brealiz);
	scl.build(RIGIDITY);				// ya declare la variable 'scl' en "./general.cc"
	cout << " scl.beta: "<<scl.beta << endl;
	cout << " rigid: "<< RIGIDITY << endl;
	cout << " frac_gy: "<< FRAC_GYROPERIOD << endl;
	cout << " npoints: " << NPOINTS << endl;
	cout << " ATOL: "<< ATOL << endl;
	cout << " RTOL: " << RTOL << endl;

	PARAMS par(fname_input_turb);			// inicializa el modelo de turbulencia!

	Doub h1=FRAC_GYROPERIOD, hmin=0.0;		// init stepsize and minimum-value
	Doub x1=0.0, x2=NMAX_GYROPERIODS;		// init and final x-points

	array_ori = read_orientations(fname_orientations, n_ori);	// lista de orientaciones de la veloc inicial
	//----- output objects w/10 points in output
	Output outbs;
	rhs d;						// object representing system-of-equations
	//------------------------------------------------
	atol	= ATOL;
	if(RTOL==-1){
		rtol = atol;				// as suggested in p.914
	}
	else if(RTOL>=.0){
		rtol = RTOL;
	}
	cout << " ABSOLUTE TOLERANCE = " << atol << endl;
	cout << " RELATIVE TOLERANCE = " << rtol << endl;

	ofstream file_out;
	sprintf(dir_out, "./output");

	PARAMS_TURB *pt;
	pt = (PARAMS_TURB *) malloc(sizeof(PARAMS_TURB));
	pt = &par.p_turb;
	//---------------------------------------------------------------- Bulirsch-Stoer
	double y[6];
	for(int j=0; j<n_Brealiz; j++){
		y[0]=0.; y[2]=0.; y[4]=0.;
		par.calc_dB(y);
		cout <<"dB0: " << par.dB[0] << endl;

		par.next_B_realization();
	/*	srand((*pt).sem.two[0]);
		(*pt).sem.slab[0] = rand();
		(*pt).sem.slab[1] = rand();
		(*pt).sem.slab[2] = rand();
		(*pt).sem.two[0]  = rand();
		(*pt).sem.two[1]  = rand();
		
		(*pt).fases.construir_fases_random((*pt).sem);
*/

		/*for(int i=0; i<n_ori; i++){
			// nueva pla
			cout << " i: "<< i << endl;
			sprintf(fname_out, "%s/pla_%03d.out", dir_out, i);
			cout << " f: "<< fname_out << endl;
			file_out.open(fname_out);
			outbs.build(NPOINTS);
			init_orientation(i, array_ori, ystart);		// values for initial y-values

			Odeint<StepperBS<rhs> > bsode(ystart,x1,x2,atol,rtol,h1,hmin,outbs,d,par);
			bsode.integrate();
			cout << "@@@@@@@@@@@@@@ BS @@@@@@@@@@@@ " << endl;
			output(outbs, file_out);		// imprime en ---> "cerr"
			
			cout << ystart[0] << "\n" << ystart[1] <<endl;
			cout << bsode.nok << endl;
			cout << bsode.nbad << endl;
			cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ " << endl;

			file_out.close();
		}*/
	}
	//------------------------------------------------
	cout <<" scl.rl [AU]: "<<scl.rl/AU_in_cm <<endl;
	cout <<" scl.wc [s-1]: "<<scl.wc<<endl;
	return 0;
}

