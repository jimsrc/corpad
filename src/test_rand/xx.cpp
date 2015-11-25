#include <iostream>
#include <math.h>
#include <cstdlib>
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

/*---------------------DEFs para "rann0" ---------------*/
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
/*-------------------------------------------------------*/

int main(int argc, char* argv[]){
    long SEMILLA = atoi(argv[1]);
    int N       = atoi(argv[2]);
    int random;
    double random_f;
    //srand(SEMILLA);

    printf(" SEMILLA: %ld \n", SEMILLA);
    for(int i=0;i<N;i++){
        //printf(" %d %d\n", i, rand());
        random = int(100.*rann0(SEMILLA));
        random_f = 100.*rann0(SEMILLA);
        printf(" %d %f %ld\n", i, random_f, SEMILLA);
    }

    printf(" SEMILLA: %ld \n", SEMILLA);
    return 0;
}
