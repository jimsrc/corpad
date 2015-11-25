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
                                                                                                          
        filein >> n_modos;                      // nro de modos                                           
        filein >> lambda_max;                   // escala minima de fluctuaciones [AU]                    
        filein >> lambda_min;                   // escala maxima de fluctuaciones [AU]                    
        filein >> Lc;                           // longitud de correlacion Lc [AU]                        
        filein >> Bunif;                        // [G]; Recordar 1nT = 1e-5G                              
        filein >> sigma_Bo_ratio;               // sigma^2 = (sigma_Bo_ratio) * Bo^2                      
        filein >> percent_slab;                 // sigma_SLAB^2 = percent_slab * sigma^2                  
        filein >> percent_2D;                   // sigma_2D^2 = percent_2D * sigma^2                      
        filein >> SEM_SLAB[0];                  // semilla #0 para campo SLAB                             
        filein >> SEM_SLAB[1];                                                                            
        filein >> SEM_SLAB[2];                                                                            
        filein >> SEM_2D[0];                    // semilla #0 para campo 2D                               
        filein >> SEM_2D[1];                                                                              
                                                                                                          
        // conversion de unidades                                                                         
        lambda_max = AU_in_cm * lambda_max;     // [cm]                                                   
        lambda_min = AU_in_cm * lambda_min;     // [cm]                                                   
        Lc = AU_in_cm * Lc;                     // [cm]                                                   
                                                                                                          
        filein.close();                                                                                   
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
        dB_SLAB[0] = 0.; dB_SLAB[1]=0.;                         // reseteo el vector en ceros                
                                                                                                             
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
        dB_2D[0]=0.; dB_2D[1]=0.;                               // reseteo el vector en ceros                
                                                                                                             
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
                                                                                                             
        calc_dB_SLAB(pos, phi_s, a_s, b_s, n_modos, k, Bk_SLAB, dB_SLAB);       // componente SLAB           
        calc_dB_2D(pos, phi_2d, b_2d, n_modos, k, Bk_2D, dB_2D);                // componente 2D             
                                                                                                             
        for(int i=0; i<2; i++){                         // solo las componentes "x, y" son turbulentas       
                dB[i] = dB_SLAB[i] + dB_2D[i];                                                               
        }                                                                                                    
}                                                                                                            

void construyo_angulos_random(int n_modos, long *SEM_SLAB, long *SEM_2D, double *phi_s, double *a_s, double *b_s, double *phi_2d, double *b_2d){                                                                         
        for(int j=0; j<n_modos; j++){                                                                       
                // va's para SLAB                                                                           
                phi_s[j]        = 2.0*pi * rann0(SEM_SLAB[0]);                                              
                a_s[j]          = 2.0*pi * rann0(SEM_SLAB[1]);                                              
                b_s[j]          = 2.0*pi * rann0(SEM_SLAB[2]);                                              
                // va's para 2D                                                                             
                phi_2d[j]       = 2.0*pi * rann0(SEM_2D[0]);                                                
                b_2d[j]         = 2.0*pi * rann0(SEM_2D[1]);                                                
        }                                                                                                   
}                                                                                                           

