# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include "stdlib.h"
#include <string.h>
using namespace std;
//output file as global variable
ofstream ofile;

int main(int argc, char* argv[]){
	int i, j, Nth, Nphi;
	double dphi, dth, pi=3.1415926535;
	char *outfilename;

	Nth  = 15;			// nro de valores theta
	Nphi = 15;			// nro de valores phi
	dth  = pi / Nth;		// delta theta
	dphi = 2. * pi / Nphi;		// delta phi
	
	if(argc<=1){
                cout << "Bad Usage: " << argv[0] << "read also output file on same line!" << endl;
                exit(1);
        }
        else{
                outfilename = argv[1];
        }

	ofile.open(outfilename);
	
	for(i=0; i<Nth; i++){
		for(j=0; j<Nphi; j++){
			ofile << setiosflags(ios :: showpoint | ios :: uppercase);
			ofile << setw(15) << setprecision(8) << i * dth  << " ";
			ofile << setw(15) << setprecision(8) << j * dphi << endl;
		}
	}
}
