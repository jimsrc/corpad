# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include "stdlib.h"
#include <string.h>
using namespace std;
//-----------------CTES UNIVERSALES
#define c 3.0*1e10			// [cm/s]
#define mo 1.6726*1e-24         	// PROTON [gr]
#define E_reposo 938272013.0		// PROTON [eV]
#define q 4.8032*1e-10	          	// PROTON [statC]
#define pi 3.1415926535			// [1]
#define AU_in_cm 1.5e13			// [cm]

float rann0(long&);
void leer_input_turb(string, int&, double&, double&, double&, double&, double&, double&, double&, long *, long *);
void leer_input_pla(string, double&, double&, double&, int&);
void imprime_en_pantalla(double, double, double, double);
void calc_relativ(double, double&, double&);
void derivatives(double, double, double *, double *, double *);
void runge_kutta_4(double, double, double *, double *, double *, int&, double, double);
void encabezado(long *, long *, double, double, double, double, double, double, double, ofstream&);
void calc_pitch(double *, double&);
void output(double, double *, double, double, double *, ofstream&);
void calc_freq_ciclotron(double, double, double&);
void calc_k_and_dk(double, double, int, double *, double *);
void calc_Bk_SLAB(double, double, double, int, double *, double *, double *);
void calc_Bk_2D(double, double, double, int, double *, double *, double *);
void calc_dB_SLAB(double *, double *, double *, double *, int, double *, double *, double *);
void calc_dB_2D(double *, double *, double *, int, double *, double *, double *);
void calc_dB(double *, double *, double *, double *, double *, double *, int, double *, double *, double *, double *, double *, double *);
void guarda_posic_final(double *, double, double, ofstream&);
void nueva_pla(int, int&, int *, double, double, double, double&, double *);
void guardo_en_tmax_partial(int, int *, double *, double, double, ofstream&);
void construyo_angulos_random(int, long *, long *, double *, double *, double *, double *, double *);

/*----------------------------------------------------------------------------------
--------------------------------- M A I N ------------------------------------------
-----------------------------------------------------------------------------------*/
int main(int argc, char* argv[]){
        int i, j, m, n, pla, step=0, EVERY, n_modos, N_FINAL_POSITIONS, *flag;
        double *y, *dydt, *yout, vo, beta, gamma, t, DT, T_MAX, T_MAX_partial;
        double *phi_s, *a_s, *b_s, *phi_2d, *b_2d;      // angulos aleatorios
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
                INPUT_TURBULENCE        = argv[1];
                INPUT_PLAS              = argv[2];
                INPUT_ORIENTATIONS_PLA  = argv[3];
        }
        //ofstream *ofile;
        //ofile.open(OUT_FILENAME);

        //este es el numero de ec's diferenciales
        n=6;

        SEM_SLAB= new long[3];                          // 3 semillas para generar campo 'dB_S'
        SEM_2D  = new long[2];                          // 2 semillas para generar campo 'DB_2D'
        dydt    = new double[n];
        y       = new double[n];
        yout    = new double[n];

        //-----------------------------------leer input del archivo 'INP_FILENAME'
        leer_input_turb(INPUT_TURBULENCE, n_modos, lambda_max, lambda_min, Lc, Bunif, sigma_Bo_ratio, percent_slab, percent_2D, SEM_SLAB, SEM_2D);
        leer_input_pla(INPUT_PLAS, RIGIDITY, FRACTION_GYROPERIOD, N_GYROPERIODS, N_FINAL_POSITIONS);
        sigma           = sqrt(sigma_Bo_ratio) * Bunif;
        sigma_S         = sqrt(percent_slab) * sigma;
        sigma_2D        = sqrt(percent_2D) * sigma;

        //----------------------------------------alocamos espacio para vectores
        // va's
        phi_s   = new double[n_modos];
        a_s     = new double[n_modos];
        b_s     = new double[n_modos];
        phi_2d  = new double[n_modos];
        b_2d    = new double[n_modos];
        // modos fourier Bk's
        Bk_SLAB = new double[n_modos];
        Bk_2D   = new double[n_modos];
        k       = new double[n_modos];
        dk      = new double[n_modos];
        // componentes turbulentas de B
        dB_SLAB = new double[3];
        dB_2D   = new double[3];
        dB      = new double[3];
        flag    = new int[N_FINAL_POSITIONS];


        //--------------------nombres de las salidas
        char ofile_name[N_FINAL_POSITIONS][60];
        ofstream *ofile;
        ofile   = new ofstream[N_FINAL_POSITIONS];

        for(i=0; i<N_FINAL_POSITIONS; i++){
                sprintf(ofile_name[i], "plas_%3.4fGV_tmax_%03d.out", RIGIDITY/1e9, i);
                ofile[i].open(ofile_name[i]);
        }

        //-----------------------------------constantes del sistema
        calc_freq_ciclotron(Bunif, RIGIDITY, omega);            // calcula freq ciclotron
        calc_relativ(RIGIDITY, beta, gamma);                    // calcula beta y gamma telativista
        vo = beta*c;                                            // calcula el modeulo de la velocidad
        RL = vo / omega;                                        // calcula radio larmor

        //-----------------------------------en pantalla
        imprime_en_pantalla(RL, omega, gamma, beta);

        //-----------------------------------valores iniciales y constantes del sistema
        t               = 0.0;
        DT              = FRACTION_GYROPERIOD / omega;
        T_MAX           = N_GYROPERIODS / omega;
        T_MAX_partial   = T_MAX / (1.*N_FINAL_POSITIONS);

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
                fscanf(ORIENTATIONS_FILE, "%lf %lf", &THETA, &PHI);     // nueva orientacion inicial
                nueva_pla(N_FINAL_POSITIONS, pla, flag, beta, PHI, THETA, t, y);                // cos(THETA) = "PITCH" (\mu)
                //--------------------------------integro paso a paso
                while(t <= T_MAX){
                        calc_dB(y, phi_s, a_s, b_s, phi_2d, b_2d, n_modos, k, Bk_SLAB, Bk_2D, dB_SLAB, dB_2D, dB);
                        derivatives(Bunif, omega, dB, y, dydt);              // halla derivadas
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

