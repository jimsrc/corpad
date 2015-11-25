#ifndef GENERAL_H
#define GENERAL_H
// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

#include <cstdlib>
//--------------------------------------CTES UNIVERSALES
#define clight          (3.0*1e10)              // [cm/s]
#define AU_in_cm        1.5e13                  // [cm]
#define nT_in_G         (1.0*1e-5)              // [1G=1e5nT]

using namespace std;

class ESCALAS{
        public:
                ESCALAS(void){};
                void build(double);
                double Bo;              // [] campo B *constante*
                double rl;              // [cm] radio larmor
                double wc;              // [s^-1] freq ciclotron
                double vel;             // [cm/s] velocidad
		double beta;
		double gamma;
        private:
                double mo;// = (1.6726*1e-24);        // [gr] masa PROTON
                double Ereposo;// = 938272013.0;      // [eV] energia de reposo PROTON
                double Z;// = +1.;                    // [e] carga (positiva) del proton en unidades de carga electronica
                double q;// = (4.8032*1e-10);         // [statC] carga PROTON
                double B;// = 5e-5;                   // [G] 5nT en Gauss
};
#endif //SCALS
