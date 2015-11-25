# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include "stdlib.h"
#include <string.h>
using namespace std;
//-----------------CTES UNIVERSALES
double c = 3.0*1e10, mo = 1.6726*1e-24;		//PROTON [cm/s], [gr]
double E_reposo = 938272013.0, q=4.8032*1e-10;		//PROTON [eV], [statC]
double pi = 3.1415926535, AU_in_cm=1.5e13;		// [1], [cm]
//-----------------DEFs para "rann0"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float rann0(long &idum){        // lo nombre diferente a "ran0" para evitar 
  long k;                       // posibles conflictos con las librerias
  float ans;
  idum ^= MASK;                   // XORing with MASK allows use of 0 and
  k = idum/IQ;                    //     other simple bit patterns for idum.
  idum = IA * (idum-k*IQ) - IR*k; // Compute idum = (IA*idum) % IM without
  if (idum < 0) idum += IM;       //     overflows by Schrage's method.
  ans = AM * idum;                // Convert idum to a floating result.
  idum ^= MASK;                   // Unmask before return.
  return ans;
}
void leer_input_turb(string INPUT_TURBULENCE, int& n_modos, double& lambda_max, double& lambda_min, double& Lc, double& Bunif, double& sigma_Bo_ratio, double& percent_slab, double& percent_2D, long *SEM_SLAB, long *SEM_2D){

	ifstream filein(INPUT_TURBULENCE.c_str());
	
	if (!filein.good()) {
                cout << "problema al abrir " << INPUT_TURBULENCE << endl;
                exit(1);
        }

        filein >> n_modos;			// nro de modos
        filein >> lambda_max;			// escala minima de fluctuaciones [AU]
        filein >> lambda_min;			// escala maxima de fluctuaciones [AU]
	filein >> Lc;				// longitud de correlacion Lc [AU]
	filein >> Bunif;			// [G]; Recordar 1nT = 1e-5G
        filein >> sigma_Bo_ratio;		// sigma^2 = (sigma_Bo_ratio) * Bo^2
        filein >> percent_slab;			// sigma_SLAB^2 = percent_slab * sigma^2
        filein >> percent_2D;			// sigma_2D^2 = percent_2D * sigma^2
	filein >> SEM_SLAB[0];			// semilla #0 para campo SLAB
	filein >> SEM_SLAB[1];			
	filein >> SEM_SLAB[2];
	filein >> SEM_2D[0];			// semilla #0 para campo 2D
	filein >> SEM_2D[1];

	// conversion de unidades
	lambda_max = AU_in_cm * lambda_max;     // [cm]
	lambda_min = AU_in_cm * lambda_min;	// [cm]
	Lc = AU_in_cm * Lc;			// [cm]

	filein.close();
}

void leer_input_pla(string INPUT_PLAS, double& RIGIDITY, double& FRACTION_GYROPERIOD, double& N_GYROPERIODS, int& N_FINAL_POSITIONS){

	ifstream filein(INPUT_PLAS.c_str());

	if (!filein.good()) {
		cout << "problema al abrir " << INPUT_PLAS << endl;
		exit(1);
	}

	filein >> RIGIDITY;			// [V]
	filein >> FRACTION_GYROPERIOD;		// DT = FRACTION_GYROPERIOD / omega;
	filein >> N_GYROPERIODS;		// T  = N_GYROPERIODS / omega
	filein >> N_FINAL_POSITIONS;

	filein.close();
}

void imprime_en_pantalla(double RL, double omega, double gamma, double beta){
        cout << " gamma = " << gamma << endl;
        cout << " beta = " << beta << endl;
	cout << " radio larmor = " << RL << " [cm]" << endl;
        cout << " omega = " << omega << " [s^-1]" << endl << endl;
}

void calc_relativ(double RIGIDITY, double &beta, double &gamma){
        gamma = pow(pow(RIGIDITY/E_reposo,2) + 1 , 0.5);
        beta = pow(1. - 1./(gamma*gamma) , 0.5);
}

void derivatives(double Bunif, double omega, double *dB, double *y, double *dydt){
        //ec's diferenciales
        dydt[0] = y[1];
        dydt[1] = omega * (y[3]*(1.+dB[2]/Bunif) - y[5] * dB[1]/Bunif);
        dydt[2] = y[3];
        dydt[3] = omega * (-y[1]*(1.+dB[2]/Bunif) + y[5] * dB[0]/Bunif);
        dydt[4] = y[5];
        dydt[5] = omega * (-y[3] * dB[0]/Bunif + y[1] * dB[1]/Bunif);
}

void runge_kutta_4(double Bunif, double omega, double *dB, double *y, double *dydx, int& n, double DT, double t){
        int i;
        double h, th, hh, h6;
	double dym[n], dyt[n], yt[n], yout[n];
        //double *dym, *dyt, *yt, *yout;

        /*dym 	= new double[6];
        dyt 	= new double[n];
        yt 	= new double[n];
	yout	= new double[n];*/

        h = DT;
        hh = h*0.5;
        h6 = h/6;
        th = t + hh;

        for (i=0; i<n; i++){
                yt[i] = y[i] + hh*dydx[i];
        }
        derivatives(Bunif, omega, dB, yt, dyt);               //computa K2
        for (i=0; i<n; i++){
                yt[i] = y[i] + hh*dyt[i];
        }
        derivatives(Bunif, omega, dB, yt, dym);               //computa K3
        for (i=0; i<n; i++){
                yt[i] = y[i] + h*dym[i];
                dym[i] += dyt[i];
        }
        derivatives(Bunif, omega, dB, yt, dyt);      //computa K4
        for (i=0; i<n; i++){
                yout[i] = y[i] + h6*(dydx[i] + dyt[i] + 2.0*dym[i]);
//              printf("i=%d> %f        %f\n", i,y[i],yout[i]); 
        }
        for(i=0; i<n; i++){
	        y[i] = yout[i];					// actualiza valores de 'y'
        }
        //delete []dym; delete []dyt; delete []yt; delete []yout;
}

void encabezado(long *SEM_SLAB, long *SEM_2D, double FRACTION_GYROPERIOD, double RIGIDITY, double PITCH, double RL, double omega, double beta, double gamma, ofstream& ofile){
	ofile << "# gamma = " << gamma << endl;
	ofile << "# beta = " << beta << endl;
	ofile << "# omega = " << omega << " [s^-1]" << endl;
	ofile << "# radio lamor = beta*c/omega = " << RL << " [cm]" << endl;
	ofile << "# ------------------particula" << endl;
	ofile << "# FRACTION_GYROPERIOD = " << FRACTION_GYROPERIOD << " " <<  endl;
	ofile << "# RIGIDITY = " << RIGIDITY << " [V]" << endl;
	ofile << "# PITCH = " << PITCH << endl;
	ofile << "# -----------------turbulencia" << endl;
	ofile << "# SEM_SLAB_0 = " << SEM_SLAB[0] << endl;
	ofile << "# SEM_SLAB_1 = " << SEM_SLAB[1] << endl;
	ofile << "# SEM_SLAB_2 = " << SEM_SLAB[2] << endl;
	ofile << "# SEM_2D_0   = " << SEM_2D[0] << endl;
	ofile << "# SEM_2D_0   = " << SEM_2D[1] << endl;
	ofile << "# ---------------------------" << endl;
	ofile << "# col1   : t [seg]"<< endl;
	ofile << "# col2-4 : x, y, z [km]" << endl;
	ofile << "# col5   : pitch_cosine" << endl;
	ofile << "# col6   : error relativo de velocidad" << endl;
	ofile << "# col7-9 : Bx, By, Bz [G]" << endl;
	ofile << "# col10  : |B| [G]" << endl;
}

void calc_pitch(double *y, double& pitch){
	double v, v_parall;
	v        = pow(y[1]*y[1] + y[3]*y[3] + y[5]*y[5], 0.5);
	v_parall = y[5];
	pitch = v_parall / v;
}

void output(double t, double *y, double vo, double Bunif, double *dB, ofstream& ofile){
	// calcula pitch-cosine
	double pitch;
	calc_pitch(y, pitch);
	// escribe archivo
        ofile << setiosflags(ios :: showpoint | ios :: uppercase);
        ofile << setw(15) << setprecision(8) << t << " ";                                               // time [s]
        ofile << setw(15) << setprecision(8) << y[0]/1e5 << " ";                                        // x [km]
        ofile << setw(15) << setprecision(8) << y[2]/1e5 << " ";                                        // y [km]
        ofile << setw(15) << setprecision(8) << y[4]/1e5 << " ";                                        // z [km]
	ofile << setw(15) << setprecision(8) << pitch << " ";						// pitch cosine
        ofile << setw(15) << setprecision(8) << pow(y[1]*y[1] + y[3]*y[3] + y[5]*y[5], 0.5)/vo - 1. << " ";     // error de veloc relat
        ofile << setw(15) << setprecision(8) << dB[0] << " ";                                           // Bx [G]
        ofile << setw(15) << setprecision(8) << dB[1] << " ";                                           // By [G]
        ofile << setw(15) << setprecision(8) << Bunif + dB[2];                                          // Bz [G]
        ofile << setw(15) << setprecision(8) << pow(dB[0]*dB[0] + dB[1]*dB[1] + (Bunif+dB[2])*(Bunif+dB[2]), 0.5) << endl;
}

void calc_freq_ciclotron(double Bunif, double RIGIDITY, double& omega)
{
        double beta, gamma;
	calc_relativ(RIGIDITY, beta, gamma);
        omega = q * Bunif / (gamma * mo * c);                   //hereda el signo de 'Bunif' inclusive
}

void calc_k_and_dk(double lambda_min, double lambda_max, int n_modos, double *k, double *dk){
	double kmin = 2. * pi / lambda_max;
	double kmax = 2. * pi / lambda_min;
	for(int i=0; i<n_modos; i++){
		k[i]  = kmin * pow(kmax/kmin, 1.*i/(n_modos-1.));
		dk[i] = k[i] * (pow(kmax/kmin, 1./(n_modos-1)) - 1.);
		//cout << "calc_k_and_sk / inside: " << k[i] << endl;
	}
}

// calcula las amplitudes de los modos en la turbulencia
void calc_Bk_SLAB(double gS, double sigma_S, double Lc, int n_modos, double *k, double *dk, double *Bk_SLAB){
	int i;
	double DENOMINADOR=0.0, FACTOR;
	for(i=0; i<n_modos; i++){
		DENOMINADOR += dk[i] / (1. + pow(k[i]*Lc, gS));
		//cout << "calc_Bk_SLAB / other " << dk[i] << endl;
	}
	for(i=0; i<n_modos; i++){
		FACTOR = dk[i] / (1. + pow(k[i]*Lc, gS)) / DENOMINADOR;
		Bk_SLAB[i] = sigma_S * pow(FACTOR, 0.5);
		//cout << "calc_Bk_SLAB / inside: " << Bk_SLAB[i] << endl;
	}
}

void calc_Bk_2D(double g2D, double sigma_2D, double Lc, int n_modos, double *k, double *dk, double *Bk_2D){
	int i;
	double DENOMINADOR=0.0, FACTOR, dV;
	for(i=0; i<n_modos; i++){
		dV = 2*pi*k[i]*dk[i];
		DENOMINADOR += dV / (1.0 + pow(k[i]*Lc, g2D));
	}
	for(i=0; i<n_modos; i++){
		dV = 2*pi*k[i]*dk[i];
		FACTOR = dV / (1.0 + pow(k[i]*Lc, g2D)) / DENOMINADOR;
		Bk_2D[i] = sigma_2D * pow(FACTOR, 0.5);
	}
}

void calc_dB_SLAB(double *pos, double *phi_s, double *a_s, double *b_s, int n_modos, double *k, double *Bk_SLAB, double *dB_SLAB){
	double COSCOS, SINSIN, FACTOR_X, FACTOR_Y;
	dB_SLAB[0] = 0.; dB_SLAB[1]=0.;				// reseteo el vector en ceros

	for(int i=0; i<n_modos; i++){
		COSCOS = cos(k[i]*pos[2] + b_s[i]) * cos(a_s[i]);
		SINSIN = sin(k[i]*pos[2] + b_s[i]) * sin(a_s[i]);
		FACTOR_X = COSCOS * cos(phi_s[i]) + SINSIN * sin(phi_s[i]);
		FACTOR_Y = COSCOS * sin(phi_s[i]) - SINSIN * cos(phi_s[i]);
		dB_SLAB[0] += Bk_SLAB[i] * FACTOR_X;
		//cout << "slab / inside: " << Bk_SLAB[i] << endl;
		dB_SLAB[1] += Bk_SLAB[i] * FACTOR_Y;
	}
}

void calc_dB_2D(double *pos, double *phi_2d, double *b_2d, int n_modos, double *k, double *Bk_2D, double *dB_2D){
	double FACTOR;
	dB_2D[0]=0.; dB_2D[1]=0.;				// reseteo el vector en ceros

	for(int i=0; i<n_modos; i++){
		FACTOR = k[i]*(pos[0] * cos(phi_2d[i]) + pos[1] * sin(phi_2d[i])) + b_2d[i];
		dB_2D[0] += Bk_2D[i] * sin(phi_2d[i]) * sin(FACTOR);
		dB_2D[1] +=-Bk_2D[i] * cos(phi_2d[i]) * sin(FACTOR);
	}
}

void calc_dB(double *y, double *phi_s, double *a_s, double *b_s, double *phi_2d, double *b_2d, int n_modos, double *k, double *Bk_SLAB, double *Bk_2D, double *dB_SLAB, double *dB_2D, double *dB){
	double pos[3] = {y[0], y[2], y[4]};
	//double *pos;
	//pos = new double[3];
	//pos[0] = y[0];
	//pos[1] = y[2];
	//pos[2] = y[4];

	calc_dB_SLAB(pos, phi_s, a_s, b_s, n_modos, k, Bk_SLAB, dB_SLAB);	// componente SLAB
	calc_dB_2D(pos, phi_2d, b_2d, n_modos, k, Bk_2D, dB_2D);		// componente 2D

	for(int i=0; i<2; i++){				// solo las componentes "x, y" son turbulentas
		dB[i] = dB_SLAB[i] + dB_2D[i];
	}
}

void guarda_posic_final(double *y, double THETA, double PHI, ofstream& OFILE){
        OFILE << setiosflags(ios :: showpoint | ios :: uppercase);
        OFILE << setw(15) << setprecision(8) << y[0]/1e5 << " ";                // [km]
        OFILE << setw(15) << setprecision(8) << y[2]/1e5 << " ";                // [km]
        OFILE << setw(15) << setprecision(8) << y[4]/1e5 << " ";                // [km]
        OFILE << setw(15) << setprecision(8) << THETA << " ";                   // [rad]
        OFILE << setw(15) << setprecision(8) << PHI << " " << endl;             // [rad]
}

void nueva_pla(int N_FINAL_POSITIONS, int& pla, int *flag, double beta, double PHI, double THETA, double& t, double *y){
	pla++;
	cout << " -----------------------pla #" << pla << endl;
	for(int i=0; i<N_FINAL_POSITIONS; i++){
		flag[i]=0;
	}
        //-------------reseteamos tiempo
        t = 0.0;
        //-------------posic en origen
        y[0] = 0;
        y[2] = 0;
        y[4] = 0;
        //-------------nueva orientacion
        double V = beta * c;
        y[1] = V * sin(THETA) * cos(PHI);
        y[3] = V * sin(THETA) * sin(PHI);
        y[5] = V * cos(THETA);
}

void guardo_en_tmax_partial(int m, int *flag, double *y, double THETA, double PHI, ofstream& OFILE){
	if(m>0 && flag[m-1]==0){
        	guarda_posic_final(y, THETA, PHI, OFILE);
		flag[m-1] = 1;
		cout << "   tmax_partial #" << m << endl;
	}
}

void construyo_angulos_random(int n_modos, long *SEM_SLAB, long *SEM_2D, double *phi_s, double *a_s, double *b_s, double *phi_2d, double *b_2d){
	for(int j=0; j<n_modos; j++){
		// va's para SLAB
		phi_s[j] 	= 2.0*pi * rann0(SEM_SLAB[0]);
		a_s[j]		= 2.0*pi * rann0(SEM_SLAB[1]);
		b_s[j]		= 2.0*pi * rann0(SEM_SLAB[2]);
		// va's para 2D
		phi_2d[j] 	= 2.0*pi * rann0(SEM_2D[0]);
		b_2d[j]		= 2.0*pi * rann0(SEM_2D[1]);
        }
}
/*----------------------------------------------------------------------------------
--------------------------------- M A I N ------------------------------------------
-----------------------------------------------------------------------------------*/
int main(int argc, char* argv[]){
        int i, j, m, n, pla, step=0, EVERY, n_modos, N_FINAL_POSITIONS, *flag;
        double *y, *dydt, *yout, vo, beta, gamma, t, DT, T_MAX, T_MAX_partial;
        double *phi_s, *a_s, *b_s, *phi_2d, *b_2d;	// angulos aleatorios
        double lambda_min, lambda_max, Lc, kmin, kmax, gS=5.0/3.0, g2D=8.0/3.0;
	double percent_slab, percent_2D, sigma_Bo_ratio, sigma, sigma_S, sigma_2D;
	double *Bk_SLAB, *Bk_2D, *Bk, *dB_SLAB, *dB_2D, *dB, *k, *dk, THETA, PHI;
	double Bunif, RIGIDITY, PITCH, RL, omega, FRACTION_GYROPERIOD, N_GYROPERIODS, pitch;
        char *INPUT_PLAS, *INPUT_TURBULENCE, *INPUT_ORIENTATIONS_PLA;
	long *SEM_SLAB, *SEM_2D;

        //Read in output file, abort if there are too few command-line arguments
	if(argc<=2){
                cout << "Bad Usage: " << argv[0] << " read also TWO input and ONE output filename on same line!" << endl;
                exit(1);
        }
        else{
		INPUT_TURBULENCE	= argv[1];
		INPUT_PLAS		= argv[2];
		INPUT_ORIENTATIONS_PLA	= argv[3];
        }
	//ofstream *ofile;
        //ofile.open(OUT_FILENAME);

	//este es el numero de ec's diferenciales
        n=6;

	SEM_SLAB= new long[3];                          // 3 semillas para generar campo 'dB_S'
	SEM_2D	= new long[2];				// 2 semillas para generar campo 'DB_2D'
        dydt	= new double[n];
        y       = new double[n];
        yout    = new double[n];

	//-----------------------------------leer input del archivo 'INP_FILENAME'
	leer_input_turb(INPUT_TURBULENCE, n_modos, lambda_max, lambda_min, Lc, Bunif, sigma_Bo_ratio, percent_slab, percent_2D, SEM_SLAB, SEM_2D);
        leer_input_pla(INPUT_PLAS, RIGIDITY, FRACTION_GYROPERIOD, N_GYROPERIODS, N_FINAL_POSITIONS);
	sigma		= sqrt(sigma_Bo_ratio) * Bunif;
	sigma_S		= sqrt(percent_slab) * sigma;
	sigma_2D	= sqrt(percent_2D) * sigma;

        //----------------------------------------alocamos espacio para vectores
	// va's
	phi_s 	= new double[n_modos];
        a_s 	= new double[n_modos];
        b_s	= new double[n_modos];
	phi_2d	= new double[n_modos];
	b_2d	= new double[n_modos];
	// modos fourier Bk's
        Bk_SLAB	= new double[n_modos];
	Bk_2D	= new double[n_modos];
	k	= new double[n_modos];
	dk	= new double[n_modos];
	// componentes turbulentas de B
	dB_SLAB = new double[3];
	dB_2D	= new double[3];
        dB 	= new double[3];
	flag	= new int[N_FINAL_POSITIONS];


	//--------------------nombres de las salidas
	char ofile_name[N_FINAL_POSITIONS][60];
	ofstream *ofile;
	ofile 	= new ofstream[N_FINAL_POSITIONS];

	for(i=0; i<N_FINAL_POSITIONS; i++){
		sprintf(ofile_name[i], "plas_%3.4fGV_tmax_%03d.out", RIGIDITY/1e9, i);
		ofile[i].open(ofile_name[i]);
	}

	//-----------------------------------constantes del sistema
        calc_freq_ciclotron(Bunif, RIGIDITY, omega);		// calcula freq ciclotron
	calc_relativ(RIGIDITY, beta, gamma);			// calcula beta y gamma telativista
	vo = beta*c;						// calcula el modeulo de la velocidad
	RL = vo / omega;					// calcula radio larmor

	//-----------------------------------en pantalla
	imprime_en_pantalla(RL, omega, gamma, beta);

	//-----------------------------------valores iniciales y constantes del sistema
        t  		= 0.0;
	DT 		= FRACTION_GYROPERIOD / omega;
	T_MAX 		= N_GYROPERIODS / omega;
	T_MAX_partial 	= T_MAX / (1.*N_FINAL_POSITIONS);	

        //-----------------------------------contruimos los 'n_modos' angulos aleatorios
	construyo_angulos_random(n_modos, SEM_SLAB, SEM_2D, phi_s, a_s, b_s, phi_2d, b_2d);

        //--------------------------------construyo las amplitudes 'Bk' de c/modo
	calc_k_and_dk(lambda_min, lambda_max, n_modos, k, dk);
	calc_Bk_SLAB(gS, sigma_S, Lc, n_modos, k, dk, Bk_SLAB);
	calc_Bk_2D(g2D, sigma_2D, Lc, n_modos, k, dk, Bk_2D);

	//--------------------------------leemos archivo donde estan las orientaciones
        //                                iniciales de las particulas
        FILE *ORIENTATIONS_FILE;
        string strfile = INPUT_ORIENTATIONS_PLA;
        ORIENTATIONS_FILE = fopen(strfile.c_str(), "r");
	pla = 0;

	while(feof(ORIENTATIONS_FILE)==0){
                fscanf(ORIENTATIONS_FILE, "%lf %lf", &THETA, &PHI);	// nueva orientacion inicial
                nueva_pla(N_FINAL_POSITIONS, pla, flag, beta, PHI, THETA, t, y);		// cos(THETA) = "PITCH" (\mu)
 	        //--------------------------------integro paso a paso
        	while(t <= T_MAX){
			calc_dB(y, phi_s, a_s, b_s, phi_2d, b_2d, n_modos, k, Bk_SLAB, Bk_2D, dB_SLAB, dB_2D, dB);
			derivatives(Bunif, omega, dB, y, dydt);		     // halla derivadas
			runge_kutta_4(Bunif, omega, dB, y, dydt, n, DT, t);  // avanza un paso en la soluc.
                	t = t + DT;

			m = int(t / T_MAX_partial);
			guardo_en_tmax_partial(m, flag, y, THETA, PHI, ofile[m-1]);
		}
	}

	// cierra archivos y borra arrays de memoria
	for(i=0; i<N_FINAL_POSITIONS; i++){
		ofile[i].close();
	}
	delete []y; delete []dydt; delete []yout; 
	delete []phi_s; delete []a_s; delete []b_s;
	delete []phi_2d; delete []b_2d;
	delete []k; delete []dk;
	delete []Bk_SLAB; delete []Bk_2D;
	delete []dB_SLAB; delete []dB_2D; delete []dB;
	delete []SEM_SLAB; delete []SEM_2D;
}

