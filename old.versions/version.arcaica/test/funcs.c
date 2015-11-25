# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include "stdlib.h"
#include <string.h>
using namespace std;
//-----------------CTES UNIVERSALES
#define c 3.0*1e10                      // [cm/s]         
#define mo 1.6726*1e-24                 // PROTON [gr]    
#define E_reposo 938272013.0            // PROTON [eV]    
#define q 4.8032*1e-10                  // PROTON [statC] 
#define pi 3.1415926535                 // [1]            
#define AU_in_cm 1.5e13                 // [cm]           

void leer_input_pla(string INPUT_PLAS, double& RIGIDITY, double& FRACTION_GYROPERIOD, double& N_GYROPERIODS, int& N_FINAL_POSITIONS){

        ifstream filein(INPUT_PLAS.c_str());

        if (!filein.good()) {
                cout << "problema al abrir " << INPUT_PLAS << endl;
                exit(1);
        }

        filein >> RIGIDITY;                     // [V]
        filein >> FRACTION_GYROPERIOD;          // DT = FRACTION_GYROPERIOD / omega;
        filein >> N_GYROPERIODS;                // T  = N_GYROPERIODS / omega
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

        /*dym   = new double[6];
        dyt     = new double[n];
        yt      = new double[n];
        yout    = new double[n];*/

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
                y[i] = yout[i];                                 // actualiza valores de 'y'
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
        ofile << setw(15) << setprecision(8) << pitch << " ";                                           // pitch cosine
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

